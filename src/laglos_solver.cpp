#include "laglos_solver.hpp"
#include <cassert>

extern "C" {
   void __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      double *in_rhol, double *in_ul, double *in_el, double *in_pl,
      double *in_rhor, double *in_ur, double *in_er, double *in_pr,
      double *in_tol, bool *no_iter,double *lambda_maxl_out,
      double *lambda_maxr_out, double *pstar, double *vstar, int *k, double *b_covolume);
}

namespace mfem
{

double fRand(double fMin, double fMax);
void random_mesh_v(const Vector &x, Vector &res);

namespace hydrodynamics
{

void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    ParGridFunction &gf, const char *title,
                    int x, int y, int w, int h, bool vec)
{
   gf.HostRead();
   ParMesh &pmesh = *gf.ParFESpace()->GetParMesh();
   MPI_Comm comm = pmesh.GetComm();

   int num_procs, myid;
   MPI_Comm_size(comm, &num_procs);
   MPI_Comm_rank(comm, &myid);

   bool newly_opened = false;
   int connection_failed;

   do
   {
      if (myid == 0)
      {
         if (!sock.is_open() || !sock)
         {
            sock.open(vishost, visport);
            sock.precision(8);
            newly_opened = true;
         }
         sock << "solution\n";
      }

      pmesh.PrintAsOne(sock);
      gf.SaveAsOne(sock);

      if (myid == 0 && newly_opened)
      {
         const char* keys = (gf.FESpace()->GetMesh()->Dimension() == 2)
                            ? "mAcRjl" : "mmaaAcl";

         sock << "window_title '" << title << "'\n"
              << "window_geometry "
              << x << " " << y << " " << w << " " << h << "\n"
              << "keys " << keys;
         if ( vec ) { sock << "vvv"; }
         sock << std::endl;
      }

      if (myid == 0)
      {
         connection_failed = !sock && !newly_opened;
      }
      MPI_Bcast(&connection_failed, 1, MPI_INT, 0, comm);
   }
   while (connection_failed);
}

template<int dim, int problem>
LagrangianLOOperator<dim, problem>::LagrangianLOOperator(ParFiniteElementSpace &h1,
                                                         ParFiniteElementSpace &l2,
                                                         ParFiniteElementSpace &l2v,
                                                         ParLinearForm *m,
                                                         bool use_viscosity,
                                                         bool mm, 
                                                         double CFL) :
   H1(h1),
   L2(l2),
   L2V(l2v),
   pmesh(H1.GetParMesh()),
   m_lf(m),
   Vsize_H1(H1.GetVSize()),
   TVSize_H1(H1.TrueVSize()),
   GTVSize_H1(H1.GlobalTrueVSize()),
   NDofs_H1(H1.GetNDofs()),
   NVDofs_H1(H1.GetNVDofs()),
   Vsize_L2(L2.GetVSize()),
   TVSize_L2(L2.TrueVSize()),
   GTVSize_L2(L2.GlobalTrueVSize()),
   NDofs_L2(L2.GetNDofs()),
   Vsize_L2V(L2V.GetVSize()),
   TVSize_L2V(L2V.TrueVSize()),
   GTVSize_L2V(L2V.GlobalTrueVSize()),
   NDofs_L2V(L2V.GetNDofs()),
   block_offsets(6),
   BdrElementIndexingArray(pmesh->GetNumFaces()),
   BdrVertexIndexingArray(pmesh->GetNV()),
   vstar_arr(pmesh->GetNumFaces()),
   // ess_tdofs(ess_tdofs),
   NE(pmesh->GetNE()),
   l2dofs_cnt(L2.GetFE(0)->GetDof()),
   h1dofs_cnt(H1.GetFE(0)->GetDof()),
   l2vdofs_cnt(L2V.GetFE(0)->GetDof()),
   num_faces(L2.GetNF()),
   num_vertices(pmesh->GetNV()),
   // Options
   use_viscosity(use_viscosity),
   mm(mm),
   CFL(CFL)

{
   block_offsets[0] = 0;
   block_offsets[1] = block_offsets[0] + Vsize_H1;
   block_offsets[2] = block_offsets[1] + Vsize_H1;
   block_offsets[3] = block_offsets[2] + Vsize_L2;
   block_offsets[4] = block_offsets[3] + Vsize_L2V;
   block_offsets[5] = block_offsets[4] + Vsize_L2;

   // This way all face indices that do not correspond to bdr elements will have val -1.
   // I.e. if global face index k corresponds to an interior face, then BdrElementIndexingArray[k] = -1.
   BdrElementIndexingArray = -1.;
   CreateBdrElementIndexingArray();

   // Similar to BdrElementIndexingArray, set default to -1.
   BdrVertexIndexingArray = -1.;
   CreateBdrVertexIndexingArray();

   // Build mass vector
   m_hpv = m_lf->ParallelAssemble();

   // Set size of intermediate face velocities vector
   v_face_intermediate.SetSize(num_faces * dim);
   v_face_intermediate = 0.;

   // Print some Dimension information
   // cout << "Vsize_H1: " << Vsize_H1 << endl;
   // cout << "TVSize_H1: " << TVSize_H1 << endl;
   // cout << "GTVSize_H1: " << GTVSize_H1 << endl;
   // cout << "NDofs_H1: " << NDofs_H1 << endl;
}

template<int dim, int problem>
LagrangianLOOperator<dim, problem>::~LagrangianLOOperator()
{
   // delete pmesh, m_lf, m_hpv;
}

template<int dim, int problem>
double LagrangianLOOperator<dim, problem>::GetCFL()
{
   return this->CFL;
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::SetCFL(const double &_CFL)
{
   this->CFL = _CFL;
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::CalculateTimestep(const Vector &S)
{
   /* TODO: Fix, issue with both timestep dependence on d and d dependence on timestep */
   // int n = m_hpv->Size(); // NDofs_L2
   double t_min = 1.;
   double t_temp = 0;
   double mi = 0;
   int ci = 0, cj = 0;

   Array<int> fids, oris;
   Vector c(dim), n(dim);
   Vector U_i(dim+2), U_j(dim+2);
   double c_norm = 0., d=0., temp_sum = 0.;

   mfem::Mesh::FaceInformation FI;

   for (ci; ci < L2.GetNE(); ci++) // Cell iterator
   {
      temp_sum = 0.;
      mi = m_hpv->Elem(ci);
      assert(mi > 0); // Assumption, equation (3.6)

      GetCellStateVector(S, ci, U_i);

      H1.ExchangeFaceNbrData();

      switch (dim)
      {
         case 1: 
         {
            pmesh->GetElementVertices(ci, fids);
            break;
         }
         case 2:
         {
            pmesh->GetElementEdges(ci, fids, oris);
            break;
         }
         case 3:
         {
            pmesh->GetElementFaces(ci, fids, oris);
         }
      }

      for (int j=0; j < fids.Size(); j++) // Face iterator
      {
         CalcOutwardNormalInt(S, ci, fids[j], c);

         c_norm = c.Norml2();
         n = c;
         n /= c_norm;

         FI = pmesh->GetFaceInformation(fids[j]);

         if (FI.IsInterior())
         {
            // Get index information/state vector for second cell
            if (ci == FI.element[0].index) { 
               cj = FI.element[1].index; 
            }
            else { 
               cj = FI.element[0].index; 
            }

            GetCellStateVector(S, cj, U_j);

            // viscosity contribution
            d = compute_lambda_max(U_i, U_j, n) * c_norm; 

            temp_sum += d/mi;
         }
      }
      
      t_temp = CFL / temp_sum;

      if (t_temp < t_min && t_temp > 1e-12) { t_min = t_temp; }
   } // End cell iterator

   this->timestep = t_min;
}

template<int dim, int problem>
double LagrangianLOOperator<dim, problem>::GetTimestep()
{
   return timestep;
}

template<int dim, int problem>
double LagrangianLOOperator<dim, problem>::GetTimeStepEstimate(const Vector &S)
{
   // Todo: implement parallelized time step estimator here.
   return 0.001;
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::IterateOverCells()
{
   cout << "-----Iterating over cells-----\n";
   Array<int> fids, fids2, oris, oris2;
   Array<int> verts;

   int ci = 0, cj = 0;
   double d = 0., c_norm = 0.;

   mfem::Mesh::FaceInformation FI;

   /* Preliminary info */
   cout << "GetNE: " << pmesh->GetNE() << endl;
   cout << "GetNBE: " << pmesh->GetNBE() << endl;
   cout << "GetNEdges: " << pmesh->GetNEdges() << endl;
   cout << "GetNumFaces: " << pmesh->GetNumFaces() << endl;
   
   for (ci; ci < L2.GetNE(); ci++) // cell iterator
   {
      cout << "We are on cell: " << ci << endl;

      pmesh->GetElementEdges(ci, fids, oris);

      // iterate over first fids
      for (int j=0; j < fids.Size(); j++)
      {
         cout << "############ fids(" << j << "): " << fids[j] << endl;
         pmesh->GetEdgeVertices(fids[j], verts);
         for (int k=0; k < verts.Size(); k++)
         {
            cout << "Verts(" << k << "): " << verts[k] << endl;
         }
      }
      // for (int j=0; j < fids2.Size(); j++)
      // {
      //    cout << "############ fids2(" << j << "): " << fids2[j] << endl;
      // }
   }
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::CreateBdrElementIndexingArray()
{
   // cout << "Constructing BdrElementIndexingArray:\n";
   for (int i = 0; i < pmesh->GetNBE(); i++)
   {
      int index = pmesh->GetBdrElementEdgeIndex(i);
      BdrElementIndexingArray[index] = i;
      // cout << "bdr el index: " << i << ", corresponding edge index: " << index << endl;
   }
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::CreateBdrVertexIndexingArray()
{
   // cout << "Creating boundary vertex indexing array.\n";
   // Iterate over boundary elements, grab edges, and fill in the corresponding
   // boundary attribute for the adjacent vertices
   // 3DTODO: Will need to modify this for faces instead of edges
   Array<int> fids, oris;
   Array<int> verts;

   for (int i = 0; i < pmesh->GetNBE(); i++)
   {
      int bdr_attr = pmesh->GetBdrAttribute(i);
      pmesh->GetBdrElementEdges(i, fids, oris);

      // iterate over boundary faces
      for (int j = 0; j < fids.Size(); j++)
      {
         // cout << "face: " << fids[j] << endl;
         pmesh->GetEdgeVertices(fids[j], verts);
         for (int k = 0; k < verts.Size(); k++)
         {
            int index = verts[k];
            // cout << "vert index: " << index << endl;
            switch (problem)
            {  
               case 7:
               {
                  // Replace the bdr attribute in the array as long as it is not
                  // the dirichlet condition (For Saltzman Problem)
                  if (BdrVertexIndexingArray[index] != 1)
                  {
                     BdrVertexIndexingArray[index] = bdr_attr;
                     // cout << "node " << index << " has bdr attr: " << bdr_attr << endl;
                  }
                  break;
               }
               default:
               {
                  BdrVertexIndexingArray[index] = bdr_attr;
               }
            } // End switch
         } // end vertex iterator
      } // end boundary faces
   } // end boundary elements
}

template<int dim, int problem>
bool LagrangianLOOperator<dim, problem>::IsBdrVertex(const int & node)
{
   if (BdrVertexIndexingArray[node] == -1)
   {
      return false;
   }
   return true;
}


template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::CreateVStarArr(const Vector &S)
{
   // Function is called at the beginning of compute_node_velocity_RP
   mfem::Mesh::FaceInformation FI;
   int c, cp;
   Vector Uc(dim+2), Ucp(dim+2), c_vec(dim), n_vec(dim);
   double in_rhol, in_ul, in_el, in_pl, in_rhor, in_ur, in_er, in_pr;

   for (int face = 0; face < pmesh->GetNumFaces(); face++)
   {
      FI = pmesh->GetFaceInformation(face);

      if (FI.IsInterior())
      {
         c = FI.element[0].index;
         cp = FI.element[1].index;  
         GetCellStateVector(S,c,Uc);
         GetCellStateVector(S,cp,Ucp);
         CalcOutwardNormalInt(S, c, face, c_vec);
         c_vec *= 2.; // Need to accomodate for the 1/2 in 4.7 that isn't present in 5.7a
         n_vec = c_vec;
         n_vec /= c_vec.Norml2();

         /* Compute vstar here. */
         in_rhol = 1. / Uc[0];
         in_ul = velocity(Uc) * n_vec; 
         in_el = specific_internal_energy(Uc);
         in_pl = pressure(Uc);

         in_rhor = 1. / Ucp[0]; 
         in_ur = velocity(Ucp) * n_vec; 
         in_er = specific_internal_energy(Ucp);
         in_pr = pressure(Ucp);

         // double b_covolume = 0.1/max(in_rhol,in_rhor);
         double b_covolume = 0.;

         double in_tol = 0.0000000000001,
                  lambda_maxl_out = 0.,
                  lambda_maxr_out = 0.,
                  pstar = 0.,
                  vstar = 0.;
         bool no_iter = true; 
         int k = 0; // Tells you how many iterations were needed for convergence

         __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
            &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,
            &no_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&vstar,&k,&b_covolume);
         
         vstar_arr[face] = vstar;
      }

   } // End face loop
}


template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::MakeTimeStep(Vector &S, double & t, double dt)
{
   Vector S_new = S; // We need a place to store updated cell values.

   // Retrieve data from Monolithic block vector S
   Vector* sptr = const_cast<Vector*>(&S_new);

   ParGridFunction x_gf, mv_gf_new;
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);
   mv_gf_new.MakeRef(&H1, *sptr, block_offsets[1]);
   
   // Variables needed for iteration
   Vector val(dim+2), c(dim), n(dim), U_i(dim+2), U_j(dim+2), sums(dim+2);
   Array<int> fids, fids2, oris, oris2, verts; // TODO: Remove fids2, oris2, These are for testing.

   int ci = 0, cj = 0;
   double d, c_norm;

   mfem::Mesh::FaceInformation FI;
   H1.ExchangeFaceNbrData();

   for (ci; ci < L2.GetNE(); ci++) // Cell iterator
   {
      GetCellStateVector(S, ci, U_i);
      val = U_i;

      switch (dim)
      {
         case 1: 
         {
            pmesh->GetElementVertices(ci, fids);
            break;
         }
         case 2:
         {
            pmesh->GetElementEdges(ci, fids, oris);
            break;
         }
         case 3:
         {
            pmesh->GetElementFaces(ci, fids, oris);
         }
      }

      DenseMatrix F_i = flux(U_i);
      sums = 0.;
      for (int j=0; j < fids.Size(); j++) // Face iterator
      {
         CalcOutwardNormalInt(S, ci, fids[j], c);

         c_norm = c.Norml2();
         n = c;
         n /= c_norm;
         FI = pmesh->GetFaceInformation(fids[j]);

         if (FI.IsInterior())
         {
            // Get index information/state vector for second cell
            if (ci == FI.element[0].index) { 
               cj = FI.element[1].index; 
            }
            else { 
               cj = FI.element[0].index; 
            }
            GetCellStateVector(S, cj, U_j); 

            // flux contribution
            DenseMatrix dm = flux(U_j);
            dm += F_i; 
            Vector y(dim+2);
            dm.Mult(c, y);
            sums -= y;

            /* viscosity contribution */
            if (use_viscosity)
            {
               d = compute_lambda_max(U_i, U_j, n) * c_norm; 
               Vector z = U_j;
               z -= U_i;
               sums.Add(d, z);
            }
         }
         else
         {
            assert(FI.IsBoundary());
            Vector y_temp(dim+2), y_temp_bdry(dim+2), U_i_bdry(dim+2);

            F_i.Mult(c, y_temp);
            U_i_bdry = U_i;
            y_temp_bdry = 0.;

            /* Enforce Boundary Conditions */
            switch (problem)
            {
               case 7: // Saltzman 
               {
                  // Check bdry flag
                  int BdrElIndex = BdrElementIndexingArray[fids[j]];
                  int bdr_attribute = pmesh->GetBdrAttribute(BdrElIndex);

                  if (bdr_attribute == 1) // Dirichlet (star states "ghost")
                  {
                     double _rho = 3.9992502342988532;
                     double _p = 1.3334833281256551;
                     Vector _v(dim);
                     _v = 0.;
                     _v[0] = 1.; //e_x
                     Vector _v_neg = _v, _vp = _v;
                     _v_neg *= -1.;
                     _vp *= _p;

                     DenseMatrix F_i_bdry(dim+2, dim);
                     F_i_bdry.SetRow(0, _v_neg);
                     for (int i = 0; i < dim; i++)
                     {
                        F_i_bdry(i+1, i) = _p;
                     }
                     F_i_bdry.SetRow(dim+1, _vp);

                     F_i_bdry.Mult(c, y_temp_bdry);
                  }

                  else if (bdr_attribute == 2) // slip
                  {
                     // Negate velocity
                     for (int _it = 0; _it < dim; _it++)
                     {
                        U_i_bdry[_it + 1] = U_i_bdry[_it + 1] * -1;
                     }
                     DenseMatrix F_i_slip = flux(U_i_bdry);
                     F_i_slip.Mult(c, y_temp_bdry);
                  }
                  else
                  {
                     cout << "invalid boundary attribute: " << bdr_attribute << endl;
                     cout << "cell : " << ci << endl;
                     cout << "face: " << fids[j] << endl;  
                     y_temp *= 2; 
                  }

                  y_temp += y_temp_bdry;
                  break;
               }

               default:
               {
                  y_temp *= 2;
                  break;
               }
            } // End switch case

            // Add in boundary contribution
            sums -= y_temp;
         } // End boundary face        

      } // End Face iterator

      sums *= dt;
      sums /= m_hpv->Elem(ci);
      val += sums;

      SetCellStateVector(S_new, ci, val);
      
   } // End cell iterator

   if (mm)
   {
      MoveMesh(S, x_gf, mv_gf_new, t, dt);
   }

   S = S_new;
   t += dt;

   CheckMassConservation(S);
} 


template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::GetCellStateVector(const Vector &S, 
                                              const int cell, 
                                              Vector &U)
{
   U.SetSize(dim + 2);

   // This happens when a face is a boundary face, one of the cells index is -1
   if (cell == -1) {
      U = 0.;
      return;
   }

   // Retrieve information from BlockVector S
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction sv_gf, v_gf, ste_gf;
   sv_gf.MakeRef(&L2, *sptr, block_offsets[2]);
   v_gf.MakeRef(&L2V, *sptr, block_offsets[3]);
   ste_gf.MakeRef(&L2, *sptr, block_offsets[4]);

   // Retrieve specific volume on cell
   U[0] = sv_gf.Elem(cell);

   // Retrieve cell velocity 
   // Here the corresponding gridfunction is a stacked vector:
   // Ex: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
   for (int i = 0; i < dim; i++)
   {
      int index = cell + i*NDofs_L2V;
      U[i+1] = v_gf.Elem(index);
   }

   // Retrieve cell specific total energy
   U[dim+1] = ste_gf.Elem(cell);
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::SetCellStateVector(Vector &S_new,
                                              const int cell,
                                              const Vector &U)
{
   // cout << "Setting cell state vector\n";
   Array<int> dofs;
   Array<int> sub_dofs;
   dofs.Append(cell);
   sub_dofs.Append(1);

   // Get current state grid functions
   Vector* sptr = const_cast<Vector*>(&S_new);
   ParGridFunction sv_gf, v_gf, ste_gf;
   sv_gf.MakeRef(&L2, *sptr, block_offsets[2]);
   v_gf.MakeRef(&L2, *sptr, block_offsets[3]);
   ste_gf.MakeRef(&L2, *sptr, block_offsets[4]);

   // Set subvectors according to U
   sv_gf.SetSubVector(dofs, U[0]);
   ste_gf.SetSubVector(dofs, U[dim+1]);

   for (int i = 1; i < dim; i++)
   {
      dofs.Append(cell + i*NDofs_L2V);
      sub_dofs.Append(i+1);
   }

   Vector vel;
   // cout << "sub_dofs:\n";
   // sub_dofs.Print(cout);
   // cout << "dofs:\n";
   // dofs.Print(cout);
   U.GetSubVector(sub_dofs, vel);
   v_gf.SetSubVector(dofs, vel);

   // Sync update gridfunctions
   sv_gf.SyncAliasMemory(S_new);
   v_gf.SyncAliasMemory(S_new);
   ste_gf.SyncAliasMemory(S_new);
}

/*
*
* cij computation
*
*/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::CalcOutwardNormalInt(const Vector &S, const int cell, const int face, Vector & res)
{
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction x_gf, mv_gf;
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);
   // Fill nodes with mesh nodes for retrieval
   const int nNodes = x_gf.Size() / dim; // TODO: Implement class variables for this

   Array<int> face_dofs;
   Vector a(dim);
   Vector face_node(dim); // Temp
   res = 0.;

   // Get int der shape functions
   Vector shapes = GetIntDerRefShapeFunctions();
   assert (dim == 2); // "This code only works for dim==2."

   H1.GetFaceDofs(face, face_dofs);

   // Compute C_{cc'} = 1/2 \Sum_{i\in {1:3}} [a_i^n \int_0^1 \theta_i(\xi) d\xi]
   for (int d = 0; d < face_dofs.Size(); d++)
   {
      a = 0.;
      for (int _dim = 0; _dim < dim; _dim++)
      {
         int index = _dim * nNodes + face_dofs[d];
         a[_dim] = x_gf(index);
      }
      if (d == face_dofs.Size() - 1)
      {
         face_node = a;
      }

      // Now that vectors are retrieved, sum
      a *= shapes[d];
      add(res, a, res);
   }

   Orthogonal(res);
   res *= 0.5;

   /* Ensure orientation of normal */
   const auto FI = pmesh->GetFaceInformation(face);

   if (FI.IsInterior()) // Interior face
   {
      // Orientation of the normal vector depends on the indices of
      // the cells that share that face.  The normal points towards
      // the cell with the lower global index.  This normal must be
      // flipped if the cell we are working with is the one of the
      // lower index.
      if (cell == min(FI.element[0].index, FI.element[1].index))
      {
         res *= -1.;
      }
   }
   else // Boundary face
   {
      // By default, boundary normals point toward cell.  This always
      // needs to be reversed to compute outward normal.
      assert(FI.IsBoundary());
      res *= -1.;
   }

   // Output normal
   // {
   //    cout << "n_vecide cell: " << cell << endl;
   //    cout << "face: " << face << endl;
   //    cout << "Face node coords: ";
   //    face_node.Print(cout);
   //    cout << "Normal vector: ";
   //    res.Print(cout);
   //    cout << endl;
   // }
}

template<int dim, int problem>
double LagrangianLOOperator<dim, problem>::cosineSimilarity(const Vector &v1, const Vector &v2)
{
   // cout << "computing cosine similarity\n";
   // cout << "v1 * v2: " << v1 * v2 << endl;
   double val = 0.;

   // Only compute cosine similarity if we have nonzero vectors
   if (v1.Norml2() != 0 && v2.Norml2() != 0)
   {
      val = (v1 * v2) / (v1.Norml2() * v2.Norml2());
   }
   // cout << "cosine similarity is: " << val << endl;
   return val;
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::Orthogonal(Vector &v)
{
   if (dim == 2)
   {
      double x = v(0), y = v(1);
      v(0) = -1 * y;
      v(1) = x;
   }
}

template<int dim, int problem>
Vector LagrangianLOOperator<dim, problem>::GetIntDerRefShapeFunctions() 
{
   /* 
   Test stuff 
   Can use this to evaluate shape functions on the reference face.
   */
   // const FiniteElement * test_fe_p = H1FEC.FiniteElementForGeometry(Geometry::SEGMENT); 
   // cout << "FE derive type: " << test_fe_p->GetDerivType() << endl;
   // const IntegrationRule ir = test_fe_p->GetNodes();
   // cout << "num integration points: " << ir.GetNPoints() << endl;

   // for (int p = 0; p < ir.GetNPoints(); p++)
   // {
   //    const IntegrationPoint &ip = ir.IntPoint(p);
   //    cout << "ip x: " << ip.x << ", ip y: " << ip.y << ", ip weight: " << ip.weight << endl;
   //    DenseMatrix dm(3,2);
   //    test_fe_p->CalcDShape(ip, dm);
   //    dm.Print(cout);
   // }
   assert (dim == 2); // "This code only works for dim==2."

   Vector r(3);
   r[0] = -1.;
   r[1] = 1.;
   r[2] = 0.;

   return r;
}

/*
*
* Enforce Boundary Conditions
*  Enforces on v_gf
*
*/
// template<int dim, int problem>
// void LagrangianLOOperator<dim, problem>::EnforceBoundaryConditions(Vector &S)
// {
//    cout << "enforcing BCs\n";
//    switch (problem)
//    {
//       case 7: // Saltzman
//       {
//          // v = e_x on left boundary (bdr_marker = 1)
//          // v.n = 0 elsewhere (bdr_marker = 2)

//          Vector* sptr = const_cast<Vector*>(&S);
//          ParGridFunction v_gf;
//          v_gf.MakeRef(&L2, *sptr, block_offsets[3]);

//          Array<int> left_cells;
 
//          /* iterate over boundary elements */
//          for (int bdr_el = 0; bdr_el < pmesh->GetNBE(); bdr_el +=1)
//          {
//             int el = -1;
//             int _info = -1;
//             int bdr_attribute = pmesh->GetBdrAttribute(bdr_el);
//             pmesh->GetBdrElementAdjacentElement(bdr_el, el, _info);

//             // cout << "bdr_el: " << bdr_el << endl;
//             // cout << "bdr attr: " << pmesh->GetBdrAttribute(bdr_el) << endl;
//             // cout << "Adjacent element: " << el << endl;
//             if (bdr_attribute == 2)
//             {
//                // Get outward normal vector
//                Vector n(dim);
//                CalcOutwardNormalInt(S, el, bdr_el, n);
//                n /= n.Norml2();

//                // Get velocity for cell el
//                Vector vel(dim);
//                for (int i = 0; i < dim; i++)
//                {
//                   int index = el + i*NDofs_L2V;
//                   vel[i] = v_gf.Elem(index);
//                }
//                // cout << "Printing old velocity for cell: " << el << endl;
//                // vel.Print(cout);

//                // cout << "Printing normal for cell: " << el << endl;
//                // n.Print(cout);

//                // Compute tangentional velocity
//                double mag = vel * n;
//                Vector norm_comp = n;
//                norm_comp *= mag;
//                Vector tan_comp(2);
//                subtract(vel, norm_comp, tan_comp);

//                double _result = tan_comp * n;
//                if (_result > pow(10.,-12))
//                {
//                   cout << "Printing new velocity for cell: " << el << endl;
//                   vel.Print(cout);
//                   cout << "Should be 0: " << _result << endl;
//                }

//                // Set velocity to tangentional velocity
//                // for (int i = 0; i < dim; i++)
//                // {
//                //    int index = el + i*NDofs_L2V;
//                //    v_gf.Elem(index) = tan_comp[i];
//                // }
//             }
//             else
//             {
//                assert(bdr_attribute == 1);
//                cout << "Element " << el << " has boundary marker 1.\n";
//                left_cells.Append(el);
//             }
//          }

//          // Handle leftmost boundary last to avoid issues with top left 
//          // and bottom left corner cells which have both a boundary marked 
//          // 1 and boundary marked 2
//          for (int i = 0; i < left_cells.Size(); i++)
//          {
//             int cell = left_cells[i];
//             cout << "Setting leftmost BC for element: " << cell << endl;
//             for (int j = 0; j < dim; j++)
//             {
//                int index = cell + j*NDofs_L2V;
//                if (j == 0)
//                {
//                   v_gf.Elem(index) = 1.;
//                }
//                else
//                {
//                   v_gf.Elem(index) = 0.;
//                }
//             }
//          }
//       }
//       default:
//       {
//          return;
//       }
//    }
// }

/* 
*
* MESH MOTION
*
*/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::MoveMesh(Vector &S, GridFunction & x_gf, GridFunction & mv_gf_new,const double & t, const double & dt)
{
   cout << "Moving mesh\n";
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction mv_gf_old;
   mv_gf_old.MakeRef(&H1, *sptr, block_offsets[1]);

   switch (problem)
   {
      case 4: // Noh problem
      {
         VectorFunctionCoefficient velocity_coeff(dim, InitialValues<problem, dim>::v0);
         // VectorFunctionCoefficient velocity(dim, random_mesh_v);
         velocity_coeff.SetTime(t);
         mv_gf_old.ProjectCoefficient(velocity_coeff);
         break;
      }
      default: // Move mesh according to paper
      {
         // Compute intermediate face velocities according to (5.7)
         mfem::Mesh::FaceInformation FI;
         int c, cp;
         Vector Uc(dim+2), Ucp(dim+2), c_vec(dim), n_vec(dim), Vf(dim); 
         double d, c_norm, F;

         for (int face = 0; face < num_faces; face++) // face iterator
         {
            Vf = 0.;
            FI = pmesh->GetFaceInformation(face);
            c = FI.element[0].index;
            cp = FI.element[1].index;

            GetCellStateVector(S, c, Uc);

            if (FI.IsInterior())
            {
               // this vector is only needed for interior faces
               GetCellStateVector(S, cp, Ucp);
               // Get normal, d, and |F|
               CalcOutwardNormalInt(S, c, face, c_vec);
               double c_norm = c_vec.Norml2();
               c_vec *= 2.; // Needed to accomodate for the 1/2 in 4.7 that isn't present in 5.7a
               F = c_vec.Norml2();
               n_vec = c_vec;
               n_vec /= F;

               d = compute_lambda_max(Uc, Ucp, n_vec) * c_norm;

               Vf = velocity(Uc);
               Vf += velocity(Ucp);
               Vf *= 0.5;
               double coeff = d * (Ucp[0] - Uc[0]) / F;
               Vf.Add(coeff, n_vec);
            }
            else 
            {
               assert(FI.IsBoundary());
               Vf = velocity(Uc);             
            }

            // Put face velocity into object
            for (int i = 0; i < dim; i++)
            {
               int index = face + i * num_faces;
               // cout << "face: " << face << ", index: " << index << endl;
               v_face_intermediate[index] = Vf[i];
            }

         } // End face iterator

         // Construct ti, face normals
         // ti = [0. 0.]^T in 2D
         // 3DTODO: Modify this according to seciton 5.2

         compute_node_velocity_LS(S, t, dt);
         Array<int> row;

         cout << "Computing velocity on faces now\n";

         // TODO: Compute velocities on faces now
         compute_face_velocity(S, dt);

         // Since we cannot use Serendipity elements, we must update cell center velocities
         Vector Vc(dim);
         Vector node_v(dim);
         Array<int> verts;
         for (int ci = 0; ci < L2.GetNE(); ci++)
         {
            GetCellStateVector(S, ci, Uc);
            Vc = velocity(Uc);
            // Vc = 0.;
            // pmesh->GetElementVertices(ci, verts);
            // for (int j = 0; j < verts.Size(); j++)
            // {
            //    get_node_velocity(S, verts[j], node_v);
            //    Vc += node_v;
            // }
            // Vc /= 4.;
            int cell_vdof = NVDofs_H1 + num_faces + ci;
            update_node_velocity(S, cell_vdof, Vc);
         }

      }
      break;
   }

   cout << "Finally, move x_gf according to mv_gf.\n";
   // Finally, move x_gf according to mv_gf
   add(x_gf, dt, mv_gf_old, x_gf);
   mv_gf_new = mv_gf_old;
   pmesh->NewNodes(x_gf, false);
   // cout << "Printing MV_gf at time " << t << endl;
   // mv_gf_new.Print(cout);
   cout << "Done moving mesh\n";
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::CheckMassConservation(const Vector &S)
{
   Vector U_i(dim + 2);
   
   for (int ci = 0; ci < L2.GetNE(); ci++)
   {
      const double m = m_hpv->Elem(ci);
      const double k = pmesh->GetElementVolume(ci);
      GetCellStateVector(S, ci, U_i);

      if (abs(k / U_i[0] - m) > pow(10, -13))
      {
         cout << "cell: " << ci << endl;
         cout << "MASS CONSERVATION BROKEN!!\n";
         cout << "k / U_i[0] - m = " << k / U_i[0] - m << endl;
         cout << "K: " << k << endl;
         cout << "T: " << U_i[0] << endl;
         cout << "m: " << m << endl;
         cout << "K/T: " << k / U_i[0] << endl;
         cout << endl;
      }
   }
}

template<int dim, int problem>
double LagrangianLOOperator<dim, problem>::CheckMassLoss(const Vector &S)
{
   // cout << "Checking Mass Conservation\n";
   Vector U_i(dim + 2);

   double num = 0., denom = 0.;
   
   for (int ci = 0; ci < L2.GetNE(); ci++)
   {
      const double m = m_hpv->Elem(ci);
      const double k = pmesh->GetElementVolume(ci);
      GetCellStateVector(S, ci, U_i);
      num += abs(k / U_i[0] - m);
      denom += abs(m);
   }
   // cout << "Total relative mass loss: " << num / denom << endl;
   return num / denom;
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::get_intermediate_face_velocity(const int & face, Vector & vel)
{
      // Retrieve face velocity from object
      for (int i = 0; i < dim; i++)
      {
         int index = face + i * num_faces;
         vel[i] = v_face_intermediate[index];
      }
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::tensor(const Vector & v1, const Vector & v2, DenseMatrix & dm)
{
   const int v1_len = v1.Size(), v2_len = v2.Size();
   for (int i = 0; i < v1_len; i++)
   {
      for (int j = 0; j < v2_len; j++)
      {
         dm.Elem(i,j) = v1[i] * v2[j];
      }
   }
   // cout << "First vector:\n";
   // v1.Print(cout);
   // cout << "Second vector:\n";
   // v2.Print(cout);
   // cout << "Resulting dm:\n";
   // dm.Print(cout); 
}

// cell weighted average
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::compute_node_velocity_cwa(Vector &S, const double & t, const double & dt)
{
   // Still need face information to implement slip BCs
   mfem::Mesh::FaceInformation FI;
   Table * edge_vertex = pmesh->GetEdgeVertexTable(); // How to iterate over faces attached to nodes
   Table vertex_edge;
   // Ref: https://mfem.org/howto/nav-mesh-connectivity/
   Transpose(*edge_vertex, vertex_edge);
   Array<int> vertex_edge_row;
   int vertex_edge_row_length;
   Vector n_vec(dim);

   Table * vertex_element = pmesh->GetVertexToElementTable();

   Array<int> row;
   int row_length;
   double denom; // change var name

   double mi = 0., Ti = 0., Ki=0.;
   Vector Uc(dim+2), Vi(dim);
   int vertex_bdr_attribute = 0;

   for (int vertex = 0; vertex < num_vertices; vertex++)
   {
      Vi = 0.; // Reset vertex velocity 
      vertex_element->GetRow(vertex, row);
      row_length = row.Size();
      denom = 0.;

      // Iterate over cells that share this vertex and average their contribution to vertex velocity
      for (int element_it = 0; element_it < row_length; element_it++)
      {
         int el_index = row[element_it];
         // Get mass of cell
         mi = m_hpv->Elem(el_index);

         // Get cell state vector
         GetCellStateVector(S, el_index, Uc);
         Ti = Uc[0];

         // Get measure of cell
         Ki = mi * Ti;
         assert(Ki >= 0.);
         Vi.Add(Ki, velocity(Uc));
         denom += Ki;
      }
      Vi /= denom;

      /* If we have a boundary node, adjust the velocity accordingly */
      // Check if the node is a boundary node
      vertex_bdr_attribute = BdrVertexIndexingArray[vertex];
      switch (problem)
      {
         case 7: // Saltzman
         {
            switch (vertex_bdr_attribute)
            {
               case 1: // Dirichlet BCs (for piston)
               {
                  Vi = 0.;

                  if (timestep_first == 0.)
                  {
                     timestep_first = timestep;
                  }
                  double _xi = t / (2*timestep_first);

                  double _psi = (4 - (_xi + 1) * (_xi - 2) * ((_xi - 2) - (abs(_xi-2) + (_xi-2)) / 2)) / 4.;

                  Vi[0] = 1. * _psi;

                  break;
               }
               case 2: // Slip BCs
               {
                  vertex_edge.GetRow(vertex, vertex_edge_row);
                  vertex_edge_row_length = vertex_edge_row.Size();
                  for (int face_it = 0; face_it < vertex_edge_row_length; face_it++)
                  {
                     int face = vertex_edge_row[face_it];
                     FI = pmesh->GetFaceInformation(face);
                     if (FI.IsBoundary())
                     {
                        CalcOutwardNormalInt(S, FI.element[0].index, face, n_vec);
                        n_vec /= n_vec.Norml2();

                        Vector Vi_normal_component = n_vec;
                        double _scale = Vi * n_vec;
                        Vi_normal_component *= _scale;

                        Vi -= Vi_normal_component;
                     }
                  }

                  break;
               }
               default:
               {
                  // do nothing
               }
            } // End Vertex attribute

            break;
         } // End case Saltzman
         default:
         {
            // no change
         }
      } // End: boundary node enforcement

      // finally, update the vertex velocity in the data structure
      update_node_velocity(S, vertex, Vi);
   } // end vertex iterator
   
}

// Riemann Problem node movement
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::compute_node_velocity_RP(Vector &S, const double & t, const double & dt)
{
   // Values needed
   mfem::Mesh::FaceInformation FI;
   Table * edge_vertex = pmesh->GetEdgeVertexTable(); // How to iterate over faces attached to nodes
   Table vertex_edge;
   // Ref: https://mfem.org/howto/nav-mesh-connectivity/
   Transpose(*edge_vertex, vertex_edge);
   Array<int> vertex_edge_row;
   int vertex_edge_row_length;
   Vector n_vec(dim);

   Table * vertex_element = pmesh->GetVertexToElementTable();

   Array<int> row;
   int row_length;

   double mi = 0., Ti = 0.;
   Vector Uc(dim+2), Vi(dim), Vi_perp(dim);
   int vertex_bdr_attribute = 0;

   // This method uses vstar_arr which have been formed using CreateVStarArr
   CreateVStarArr(S);

   // Iterate over nodes
   for (int vertex = 0; vertex < num_vertices; vertex++)
   {
      // Reset values and grab adjacent elements
      Vi = 0.; // Reset vertex velocity 
      vertex_bdr_attribute = 0; // reset boundary attribute

      vertex_element->GetRow(vertex, row);
      row_length = row.Size();
      Vector vertex_coords(2);
      get_node_position(S, vertex, vertex_coords);

      // Average velocity of adjacent cells and normalize to get a direction
      for (int element_it = 0; element_it < row_length; element_it++)
      {
         int el_index = row[element_it];

         // Get cell state vector
         GetCellStateVector(S, el_index, Uc);

         Vi += velocity(Uc);
      } // End cell iterator

      Vi /= row_length; 

      // If Vi is zero, no more computation is needed since the node will not move.
      // Hence, only compute something if the directional vector is nonzero
      if (Vi.Norml2() != 0)
      {
         Vi /= Vi.Norml2();

         // Compute orthogonal vector to normal direction
         Vi_perp = Vi;
         Orthogonal(Vi_perp);
         
         double tempCosineSimilarity = 0.;
         int closest_face1 = -1, closest_face2 = -1;
         double closest_face1_cs = -1., closest_face2_cs = -1.; // value to temporary store cosineSimilarity

         // Iterate over adjacent faces
         vertex_edge.GetRow(vertex, vertex_edge_row);
         vertex_edge_row_length = vertex_edge_row.Size();
         for (int face_it = 0; face_it < vertex_edge_row_length; face_it++)
         {
            // Get Face information
            int face = vertex_edge_row[face_it];
            FI = pmesh->GetFaceInformation(face);

            // Compute coordinates for face node
            Array<int> face_dof_row;
            H1.GetFaceDofs(face, face_dof_row);
            int face_dof = face_dof_row[2];
            Vector face_x(2);
            get_node_position(S, face_dof, face_x);

            // Compute the secant vector to each of the adjact face midpoint
            Vector secant(2);
            subtract(face_x, vertex_coords, secant);

            // Compute the cosine similarity from the reference line
            tempCosineSimilarity = abs(cosineSimilarity(secant, Vi_perp));

            // Determine if this face is one of the closest two
            // by maximizing cosine similarity
            if (tempCosineSimilarity >= closest_face2_cs)
            {
               if (tempCosineSimilarity >= closest_face1_cs)
               {
                  // We have the new closest face
                  closest_face2_cs = closest_face1_cs;
                  closest_face2 = closest_face1;

                  closest_face1_cs = tempCosineSimilarity;
                  closest_face1 = face;
               }
               else {
                  // We have the new second closest face
                  closest_face2_cs = tempCosineSimilarity;
                  closest_face2 = face;
               }
            }
         } // End face iterator

         // Average the magnitude of those vstar values  
         double calc_vel = (abs(vstar_arr[closest_face1]) + abs(vstar_arr[closest_face2])) / 2.;

         Vi *= calc_vel;

         /* If we have a boundary node, adjust the velocity accordingly */
         // Check if the node is a boundary node
         vertex_bdr_attribute = BdrVertexIndexingArray[vertex];
         switch (problem)
         {
            case 7: // Saltzman
            {
               switch (vertex_bdr_attribute)
               {
                  case 1: // Dirichlet BCs (for piston)
                  {
                     Vi = 0.;

                     if (timestep_first == 0.)
                     {
                        timestep_first = timestep;
                     }
                     double _xi = t / (2*timestep_first);

                     double _psi = (4 - (_xi + 1) * (_xi - 2) * ((_xi - 2) - (abs(_xi-2) + (_xi-2)) / 2)) / 4.;

                     Vi[0] = 1. * _psi;

                     break;
                  }
                  case 2: // Slip BCs
                  {
                     Vi *= 2.;
                     vertex_edge.GetRow(vertex, vertex_edge_row);
                     vertex_edge_row_length = vertex_edge_row.Size();
                     for (int face_it = 0; face_it < vertex_edge_row_length; face_it++)
                     {
                        // Vi *= 2.;
                        int face = vertex_edge_row[face_it];
                        FI = pmesh->GetFaceInformation(face);
                        if (FI.IsBoundary())
                        {
                           CalcOutwardNormalInt(S, FI.element[0].index, face, n_vec);
                           n_vec /= n_vec.Norml2();

                           Vector Vi_normal_component = n_vec;
                           double _scale = Vi * n_vec;
                           Vi_normal_component *= _scale;

                           Vi -= Vi_normal_component;
                        }
                     }

                     break;
                  }
                  default:
                  {
                     // do nothing
                  }
               } // End Vertex attribute

               break;
            } // End case Saltzman
            default:
            {
               // no change
            }
         } // End: boundary node enforcement

      } // End if Vi.Norml2() != 0

      // Update the vertex velocity
      update_node_velocity(S, vertex, Vi);  
   } // End vertex iteration
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::
   compute_node_velocity_LS(Vector &S, const double & t, 
                            const double & dt, 
                            const string flag, // Default NA
                            void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
   mfem::Mesh::FaceInformation FI;
   Table * edge_vertex = pmesh->GetEdgeVertexTable(); // How to iterate over faces attached to nodes
   Table vertex_edge;
   // Ref: https://mfem.org/howto/nav-mesh-connectivity/
   Transpose(*edge_vertex, vertex_edge);
   Array<int> row;
   int row_length;

   Vector Vf(dim); // Intermediate Face Velocity or exact velocity if testing
   Vector n_vec(dim); // For enforcing BCs

   bool is_boundary_node = false;

   // Iterate over vertices
   for (int vertex = 0; vertex < num_vertices; vertex++) // Vertex iterator
   {
      // Reset flag indicating boundary node
      is_boundary_node = false;

      // Construct A, B_x, B_y
      vertex_edge.GetRow(vertex, row);
      row_length = row.Size();

      // Instantiate A
      DenseMatrix A(row_length, dim + 1), At(dim+1, row_length), LHS(dim+1, dim+1);
      Vector Bx(row_length), By(row_length), vertex_x(dim), vertex_v(dim);
      Vector Bx_mod(dim+1), By_mod(dim+1); // To solve LS problem
      get_node_position(S, vertex, vertex_x);
      
      // 3DTODO: Will need to modify this for faces instead of edges in 3D
      for (int face_it = 0; face_it < row_length; face_it++) // Adjacent face iterator
      {
         // Reset face velocity object
         Vf = 0.;

         int face = row[face_it];
         FI = pmesh->GetFaceInformation(face);
         if (FI.IsBoundary()) { is_boundary_node = true; }

         // Retrieve coordinates for face node
         Array<int> face_dof_row;
         H1.GetFaceDofs(face, face_dof_row);
         int face_dof = face_dof_row[2]; // Dof corresponding to the center of the face
         Vector face_x(dim);
         get_node_position(S, face_dof, face_x);

         // Fill corresponding row of A, Bx, By
         for (int i = 0; i < dim; i++)
         {
            A(face_it, i) = face_x(i) - vertex_x(i); // Fill with nodal locations
         }
         A(face_it, dim) = 1.;

         if (flag == "testing")
         {
            // Used for make test
            test_vel(face_x, t, Vf);
         }
         else
         {
            get_intermediate_face_velocity(face, Vf);
         }

         // In either case, we set the corresponding row to the face velocity,
         // whether that is given by the intermediate face velocity or the exact
         Bx(face_it) = Vf(0);
         By(face_it) = Vf(1);
      } // End face iterator

      // We should have an overdetermined system 
      // assert(row_length > dim + 1);
      At.Transpose(A);
      Mult(At, A, LHS);
      At.Mult(Bx, Bx_mod);
      At.Mult(By, By_mod);

      // Solve linear systems A x_A = B_x, A x_B = B_y
      Vector res_x(dim+1), res_y(dim+1);

      if (is_boundary_node)
      {
         CGSolver cg;
         cg.SetOperator(LHS);
         cg.SetRelTol(1e-12);
         cg.SetAbsTol(1e-12);
         cg.SetMaxIter(1000);
         cg.SetPrintLevel(-1);
         cg.Mult(Bx_mod, res_x);
         cg.Mult(By_mod, res_y);
      } // End Boundary Node
      else
      {
         // Solve directly
         LHS.Invert();
         LHS.Mult(Bx_mod, res_x);
         LHS.Mult(By_mod, res_y);
      }
      
      // Solve for nodal velocity
      vertex_v(0) = res_x(dim); // Set to cx
      vertex_v(1) = res_y(dim); // Set to cy

      // Enforce boundary conditions
      if (is_boundary_node)
      {
         /* Enforce BCs (Copy/Paste from previous nodal computation) */
         switch (problem)
         {
            case 7: // Saltzman
            {
               // cout << "Working on boundary node: " << vertex << endl;
               is_boundary_node = false;
               for (int face_it = 0; face_it < row_length; face_it++)
               {
                  int face = row[face_it];
                  FI = pmesh->GetFaceInformation(face);
                  if (FI.IsBoundary())
                  {
                     int BdrElIndex = BdrElementIndexingArray[face];
                     int bdr_attribute = pmesh->GetBdrAttribute(BdrElIndex);

                     if (bdr_attribute == 1)
                     {
                        vertex_v = 0.;
                        
                        if (timestep_first == 0.)
                        {
                           timestep_first = timestep;
                        }
                        double _xi = t / (2*timestep_first);

                        double _psi = (4 - (_xi + 1) * (_xi - 2) * ((_xi - 2) - (abs(_xi-2) + (_xi-2)) / 2)) / 4.;

                        vertex_v[0] = 1. * _psi;

                        break;
                     }
                     else 
                     {
                        assert (bdr_attribute == 2);
                        CalcOutwardNormalInt(S, FI.element[0].index, face, n_vec);
                        n_vec /= n_vec.Norml2();

                        Vector vertex_v_normal_component(dim);
                        vertex_v_normal_component = n_vec;
                        double _scale = vertex_v * n_vec;
                        vertex_v_normal_component *= _scale;

                        vertex_v -= vertex_v_normal_component;
                     }
                  } // End boundary face
               } // End: Face iterator
            } // End: Saltzman
            default:
            {
               // No change
            }
         } // Switch case end
      }

      
      update_node_velocity(S, vertex, vertex_v);
   } // End Vertex iterator
}


template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::compute_node_velocity(Vector &S, const double & t, const double & dt)
{
   // Compute node velocities according to (5.10)
   mfem::Mesh::FaceInformation FI;
   Table * edge_vertex = pmesh->GetEdgeVertexTable(); // How to iterate over faces attached to nodes
   Table vertex_edge;
   // Ref: https://mfem.org/howto/nav-mesh-connectivity/
   Transpose(*edge_vertex, vertex_edge);
   // cout << "Printing transpose table (vertex to edge)\n";
   // vertex_edge.Print(cout);
   Array<int> row;
   int row_length;
   DenseMatrix LHS(dim), dm_tmp(dim);
   Vector RHS(dim), y(dim), Vi(dim), Vf(dim), n_vec(dim);
   bool is_boundary_node = false;
   double c_norm = 0.;

   // Vector Vf(dim);
   for (int vertex = 0; vertex < num_vertices; vertex++) // Vertex iterator
   {
      // cout << "Computing vertex velocity at vertex: " << vertex << endl;
      LHS = 0.;
      RHS = 0.;
      // cout << "vertex: " << vertex << endl;
      vertex_edge.GetRow(vertex, row);
      row_length = row.Size();
      // 3DTODO: Will need to modify this for faces instead of edges in 3D
      for (int face_it = 0; face_it < row_length; face_it++) // Adjacent face iterator
      {
         Vf = 0.;
         y = 0.;
         int face = row[face_it];
         // cout << "Using face: " << face << endl;
         // cout << "Adjacent face: " << face << endl;
         get_intermediate_face_velocity(face, Vf);
         FI = pmesh->GetFaceInformation(face);
         if (FI.IsBoundary()) { is_boundary_node = true; }

         // cout << "Face elements: 1: " << FI.element[0].index << ", 2: " << FI.element[1].index << endl;

         // Retrieve corresponding normal (orientation doesn't matter)
         CalcOutwardNormalInt(S, FI.element[0].index, face, n_vec);
         c_norm = n_vec.Norml2();
         n_vec /= c_norm;

         tensor(n_vec, n_vec, dm_tmp);
         LHS += dm_tmp;
         dm_tmp.Mult(Vf, y);
         RHS += y;
      }

      LHS.Invert();
      LHS.Mult(RHS, Vi);

      /* If we have a boundary node, adjust the velocity accordingly */
      if (is_boundary_node)
      {
         switch (problem)
         {
            case 7: // Saltzman
            {
               // cout << "Working on boundary node: " << vertex << endl;
               is_boundary_node = false;
               for (int face_it = 0; face_it < row_length; face_it++)
               {
                  int face = row[face_it];
                  FI = pmesh->GetFaceInformation(face);
                  if (FI.IsBoundary())
                  {
                     int BdrElIndex = BdrElementIndexingArray[face];
                     int bdr_attribute = pmesh->GetBdrAttribute(BdrElIndex);

                     if (bdr_attribute == 1)
                     {
                        Vi = 0.;
                        
                        if (timestep_first == 0.)
                        {
                           timestep_first = timestep;
                        }
                        double _xi = t / (2*timestep_first);

                        double _psi = (4 - (_xi + 1) * (_xi - 2) * ((_xi - 2) - (abs(_xi-2) + (_xi-2)) / 2)) / 4.;

                        Vi[0] = 1. * _psi;

                        break;
                     }
                     else 
                     {
                        assert (bdr_attribute == 2);
                        CalcOutwardNormalInt(S, FI.element[0].index, face, n_vec);
                        n_vec /= n_vec.Norml2();

                        Vector Vi_normal_component(dim);
                        Vi_normal_component = n_vec;
                        double _scale = Vi * n_vec;
                        Vi_normal_component *= _scale;

                        Vi -= Vi_normal_component;
                     }
                  } // End boundary face
               } // End: Face iterator
            } // End: Saltzman
            default:
            {
               // No change
            }
         } // Switch case end
         
      } // End: Boundary node

      update_node_velocity(S, vertex, Vi);
   } // end: vertex iterator 
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::compute_face_velocity(Vector &S, const double & dt)
{
   /* Parameters needed for face velocity calculations */
   mfem::Mesh::FaceInformation FI;
   Vector Vf(dim), c_vec(dim), n_vec(dim);
   Array<int> row;
   int face_dof = 0, face_vdof2 = 0, c = 0;

   // Iterate over faces
   for (int face = 0; face < num_faces; face++) // face iterator
   {  
      // Get intermediate face velocity, face information, and face normal
      Vector face_velocity(2);
      Vector n_vec_temp(2);
      Vf = 0.;
      get_intermediate_face_velocity(face, Vf);
      FI = pmesh->GetFaceInformation(face);

      /* adjacent corner indices */
      H1.GetFaceDofs(face, row);
      int face_vdof1 = row[1], face_vdof2 = row[0], face_dof = row[2]; // preserve node orientation discussed in appendix A

      // retrieve old face and corner locations
      Vector face_x(2), vdof1_x(2), vdof2_x(2), vdof1_v(2), vdof2_v(2);
      get_node_position(S, face_dof, face_x);
      get_node_position(S, face_vdof1, vdof1_x);
      get_node_position(S, face_vdof2, vdof2_x);

      // retrieve corner velocities
      get_node_velocity(S, face_vdof1, vdof1_v);
      get_node_velocity(S, face_vdof2, vdof2_v);

      if (FI.IsInterior())
      {
         c = FI.element[0].index;
         CalcOutwardNormalInt(S, c, face, c_vec);

         // Face area
         // https://github.com/mfem/mfem/issues/951
         Vector nor(dim);
         ElementTransformation *T = pmesh->GetFaceTransformation(face);
         T->SetIntPoint(&Geometries.GetCenter(pmesh->GetFaceBaseGeometry(face)));
         CalcOrtho(T->Jacobian(), nor);
         double face_length = nor.Norml2();

         n_vec = c_vec;
         n_vec /= n_vec.Norml2();

         // Calculate new corner locations and half  step locations
         Vector vdof1_x_new(2), vdof2_x_new(2), vdof1_x_half(2), vdof2_x_half(2);
         vdof1_x_new = vdof1_x;
         vdof1_x_new.Add(dt, vdof1_v);
         vdof1_x_half = vdof1_x;
         vdof1_x_half.Add(dt/2, vdof1_v);

         vdof2_x_new = vdof2_x;
         vdof2_x_new.Add(dt, vdof1_v);
         vdof2_x_half = vdof2_x;
         vdof2_x_half.Add(dt/2, vdof2_v);

         // calculate a_{12}^{n+1} (new tangent midpoint)
         Vector vdof12_x_new(2);
         vdof12_x_new = vdof1_x_new;
         vdof12_x_new += vdof2_x_new;
         vdof12_x_new /= 2.;

         // Compute D (A.4c)
         Vector n_vec_R(dim), n_vec_perp(dim), temp_vec(dim), temp_vec_2(dim);
         n_vec_R = n_vec;
         Orthogonal(n_vec_R);

         n_vec_perp = n_vec_R;
         n_vec_perp.Neg();

         subtract(vdof1_v, vdof2_v, temp_vec); // V1 - V2 = temp_vec

         subtract(vdof2_x_half, vdof1_x_half, temp_vec_2); // A2-A1

         Orthogonal(temp_vec_2);

         double D = dt * (temp_vec * n_vec_R) + 2 * (n_vec * temp_vec_2);

         // Compute c1 (A.4a)
         subtract(vdof2_v, vdof1_v, temp_vec); // only change temp_vec, since temp_vec_2 is same from D calculation (half step representation)
         double c1 = ( dt * (temp_vec * n_vec) - 2 * (temp_vec_2 * n_vec_perp) ) / D;

         // Compute c0 (A.4b)
         double bmn = Vf * n_vec; // calculate flux across face (5.8)
         bmn *= face_length;
         temp_vec = vdof1_x_half;
         Orthogonal(temp_vec);
         temp_vec_2 = vdof2_x_half;
         Orthogonal(temp_vec_2);
         double const1 = vdof1_v * temp_vec - vdof2_v * temp_vec_2; // V1*A1R - V2*A2R
         double const2 = vdof1_v * temp_vec_2 - vdof2_v * temp_vec; // V1*A2R - V2*A1R
         temp_vec = face_x;
         Orthogonal(temp_vec);
         subtract(vdof2_v, vdof1_v, temp_vec_2);
         double const3 = temp_vec_2 * temp_vec; // (V2 - V1) * a3n
         double c0 = (3. / D) * (bmn + const1 / 2. + const2 / 6. + 2. * const3 / 3.);
         
         // Compute V3n_perp
         subtract(vdof2_x_new, vdof1_x_new, temp_vec);
         subtract(face_x, vdof12_x_new, temp_vec_2);
         temp_vec_2.Add(c0*dt, n_vec);
         const1 = temp_vec * temp_vec_2; // numerator
         temp_vec_2 = n_vec_perp;
         temp_vec_2.Add(c1, n_vec);
         const2 = temp_vec * temp_vec_2;
         const2 *= dt; // denominator
         double V3nperp = -1. * const1 / const2;

         // Compute V3n (5.11)
         double V3n = c1 * V3nperp + c0;

         // Compute face velocity (5.11)
         // const1 = Vf * n_vec;
         // temp_vec = n_vec;
         // temp_vec *= const1;
         // subtract(Vf, temp_vec, face_velocity);
         // face_velocity.Add(V3n, n_vec);

         // Compute face velocity (Appendix A)
         face_velocity = 0.;
         face_velocity.Add(V3n, n_vec);
         face_velocity.Add(V3nperp, n_vec_perp);

         Vector face_x_new(2);
         face_x_new = face_x;
         face_x_new.Add(dt, face_velocity);

         // Check perpendicular
         subtract(face_x_new, vdof12_x_new, temp_vec);
         subtract(vdof2_x_new, vdof1_x_new, temp_vec_2);
         if (abs(temp_vec * temp_vec_2) > pow(10, -12))
         {
            cout << "temp_vec:\n";
            temp_vec.Print(cout);
            cout << "temp_vec_2:\n";
            temp_vec_2.Print(cout);
            cout << "################## Vectors are not orthogonal!";
            MFEM_ABORT("vectors are not orthogonal!\n");
         }
      }
      else
      {
         assert(FI.IsBoundary());
         // TODO: Add in boundary conditions similar to corner vertex bcs
         for (int j = 0; j < dim; j++)
         {
            face_velocity[j] = (vdof1_v[j] + vdof2_v[j]) / 2;
         }
      }

      // Lastly, put face velocity into gridfunction object
      update_node_velocity(S, face_dof, face_velocity);
   }
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::update_node_velocity(Vector &S, const int & node, const Vector & vel)
{
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction mv_gf;
   mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);

   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      mv_gf[index] = vel[i];
   }
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::get_node_velocity(const Vector &S, const int & node, Vector & vel)
{
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction mv_gf;
   mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);

   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      vel[i] = mv_gf[index];
   }
}  

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::get_node_position(const Vector &S, const int & node, Vector & x)
{
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction x_gf;
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);

   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      x[i] = x_gf[index];
   }
}

/* ProblemDescription Functions */
template<int dim, int problem>
double LagrangianLOOperator<dim, problem>::internal_energy(const Vector & U)
{
   const double &rho = 1./U[0];
   const double &e = specific_internal_energy(U);
   return rho * e;
}

template<int dim, int problem>
double LagrangianLOOperator<dim, problem>::specific_internal_energy(const Vector & U)
{
   const Vector  v   = velocity(U);
   const double &E   = U[dim + 1]; // specific total energy
   return E - 0.5 * pow(v.Norml2(), 2);
}

template<int dim, int problem>
double LagrangianLOOperator<dim, problem>::pressure(const Vector & U)
{
   switch (problem)
   {
      case 0:
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      case 7:
      default:
      {
         double gamma = InitialValues<problem,dim>::gamma_func();
         return (gamma - 1.) * internal_energy(U);
      }
   }
}

template<int dim, int problem>
double LagrangianLOOperator<dim, problem>::compute_lambda_max(const Vector & U_i,
                                                     const Vector & U_j,
                                                     const Vector & n_ij,
                                                     const string flag) // default is 'NA'
{
   double in_rhol, in_ul, in_el, in_pl, in_rhor, in_ur, in_er, in_pr;
   if (flag == "testing")
   {
      in_rhol = U_i[0];
      in_ul = U_i[1];
      in_pl = U_i[2];
      in_el = U_i[3];

      in_rhor = U_j[0]; 
      in_ur = U_j[1];
      in_pr = U_j[2];
      in_er = U_j[3];
   }
   else 
   {
      assert(flag == "NA");

      in_rhol = 1. / U_i[0];
      in_ul = velocity(U_i) * n_ij; 
      in_el = specific_internal_energy(U_i);
      in_pl = pressure(U_i);

      in_rhor = 1. / U_j[0]; 
      in_ur = velocity(U_j) * n_ij; 
      in_er = specific_internal_energy(U_j);
      in_pr = pressure(U_j);
   }

   // double b_covolume = 0.1/max(in_rhol,in_rhor);
   double b_covolume = 0.;

   double in_tol = 0.0000000000001,
          lambda_maxl_out = 0.,
          lambda_maxr_out = 0.,
          pstar = 0.,
          vstar = 0.;
   bool no_iter = false;
   int k = 0; // Tells you how many iterations were needed for convergence

   // cout << "CLM pre fortran function.\n";
   __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,
      &no_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&vstar,&k,&b_covolume);

   return std::max(std::abs(lambda_maxl_out), std::abs(lambda_maxr_out));
   // return 1.;
}

template<int dim, int problem>
Vector LagrangianLOOperator<dim, problem>::velocity(const Vector & U)
{
   Vector v;
   v.SetSize(dim);
   Array<int> dofs;
   for (int i = 0; i < dim; i++)
   {
      dofs.Append(i + 1);
   }
   U.GetSubVector(dofs, v);

   return v;
}

template<int dim, int problem>
DenseMatrix LagrangianLOOperator<dim, problem>::flux(const Vector &U)
{
   DenseMatrix result(dim+2, dim);

   const Vector v = velocity(U);
   const double p = pressure(U);

   // * is not overridden for Vector class, but *= is
   Vector v_neg = v, vp = v;
   v_neg *= -1.;
   vp *= p;

   // Set f(U) according to (2.1c)
   result.SetRow(0,v_neg);

   for (int i = 0; i < dim; i++)
   {
      result(i+1, i) = p;
   }
   result.SetRow(dim+1, vp);

   return result;
}

/* Explicit n_vectantiation */
template class LagrangianLOOperator<1, 2>;
template class LagrangianLOOperator<2, 2>;
template class LagrangianLOOperator<3, 2>;

template class LagrangianLOOperator<1, 1>;
template class LagrangianLOOperator<2, 1>;
template class LagrangianLOOperator<3, 1>;

template class LagrangianLOOperator<1, 0>;
template class LagrangianLOOperator<2, 0>;
template class LagrangianLOOperator<3, 0>;

template class LagrangianLOOperator<1, 3>;
template class LagrangianLOOperator<2, 3>;
template class LagrangianLOOperator<3, 3>;

template class LagrangianLOOperator<1, 4>;
template class LagrangianLOOperator<2, 4>;
template class LagrangianLOOperator<3, 4>;

template class LagrangianLOOperator<1, 5>;
template class LagrangianLOOperator<2, 5>;
template class LagrangianLOOperator<3, 5>;

template class LagrangianLOOperator<1, 6>;
template class LagrangianLOOperator<2, 6>;
template class LagrangianLOOperator<3, 6>;

template class LagrangianLOOperator<1, 7>;
template class LagrangianLOOperator<2, 7>;
template class LagrangianLOOperator<3, 7>;

template class LagrangianLOOperator<1, 8>;
template class LagrangianLOOperator<2, 8>;
template class LagrangianLOOperator<3, 8>;

} // end ns hydrodynamics


double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


/*
*
* Initial Conditions
*
*/
void random_mesh_v(const Vector &x, Vector &res)
{
   const double r_max = 1.;
   res[0] = fRand(0, r_max);
   res[1] = fRand(0, r_max);
   return;
}

} // end ns mfem
