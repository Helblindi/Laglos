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

// TODO: Create basic object instantiation
// template<int dim, int problem>
// LagrangianLOOperator<dim, problem>::LagrangianLOOperator() :
//    H1(new FiniteElementSpace()),
//    use_viscosity(true),
//    mm(true),
//    CFL(1.)
// {

// }

template<int dim, int problem>
LagrangianLOOperator<dim, problem>::LagrangianLOOperator(ParFiniteElementSpace &h1,
                                                         ParFiniteElementSpace &l2,
                                                         ParFiniteElementSpace &l2v,
                                                         ParFiniteElementSpace &cr,
                                                         ParLinearForm *m,
                                                         bool use_viscosity,
                                                         bool mm, 
                                                         double CFL) :
   H1(h1),
   L2(l2),
   L2V(l2v),
   CR(cr),
   CRc(CR.GetParMesh(), CR.FEColl(), 1),
   v_CR_gf(&CR),
   pmesh(H1.GetParMesh()),
   m_lf(m),
   Vsize_H1(H1.GetVSize()),
   TVSize_H1(H1.TrueVSize()),
   GTVSize_H1(H1.GlobalTrueVSize()),
   NDofs_H1(H1.GetNDofs()),
   NVDofs_H1(H1.GetNVDofs()), // Scalar vertex dofs
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
   num_faces(L2.GetNF()),
   num_vertices(pmesh->GetNV()),
   num_edges(pmesh->GetNEdges()),
   vertex_element(pmesh->GetVertexToElementTable()),
   face_element(pmesh->GetFaceToElementTable()),
   // element_face(pmesh->GetElementToFaceTable()),
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
   cout << "Initial mass vector:\n";
   m_hpv->Vector::Print(std::cout);

   // resize v_CR_gf to correspond to the number of faces
   if (dim == 1)
   {
      v_CR_gf.SetSize(num_faces);
   }
   
   // Initialize values of intermediate face velocities
   v_CR_gf = 0.;

   assert(v_CR_gf.Size() == dim * num_faces);

   // Set integration rule for Rannacher-Turek space
   // TODO: Modify this to be set by a parameter rather than hard-coded
   IntegrationRules IntRulesLo(0, Quadrature1D::GaussLobatto);
   RT_ir = IntRules.Get(CR.GetFE(0)->GetGeomType(), 2);

   // Print some Dimension information
   cout << "Vsize_H1: " << Vsize_H1 << endl;
   cout << "TVSize_H1: " << TVSize_H1 << endl;
   cout << "GTVSize_H1: " << GTVSize_H1 << endl;
   cout << "NDofs_H1: " << NDofs_H1 << endl;
   cout << "NVDofs_H1: " << NVDofs_H1 << endl;
   cout << "NDofs_L2: " << NDofs_L2 << endl;
   cout << "NDofs_L2V: " << NDofs_L2V << endl;
   cout << "Vsize_L2V: " << Vsize_L2V << endl;
   cout << "CR.GetNDofs(): " << CR.GetNDofs() << endl;
   cout << "pmesh->GetNFaces(): " << pmesh->GetNFaces() << endl;
   cout << "num_faces: " << num_faces << endl;
   cout << "num_vertices: " << num_vertices << endl;
   cout << "num_edges: " << num_edges << endl;
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
   Vector c(dim), n(dim), n_int(dim);;
   Vector U_i(dim+2), U_j(dim+2);
   double c_norm = 0., d=0., temp_sum = 0.;

   mfem::Mesh::FaceInformation FI;

   for (; ci < NDofs_L2; ci++) // Cell iterator
   {
      temp_sum = 0.;
      mi = m_hpv->Elem(ci);

      GetCellStateVector(S, ci, U_i);
      // cout << "cell state vector: \n";
      // U_i.Print(cout);
      assert(mi > 0); // Assumption, equation (3.6)

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
         CalcOutwardNormalInt(S, ci, fids[j], n_int);
         c = n_int;
         c /= 2.;
         c_norm = c.Norml2();
         n = n_int;
         double n_norm = n.Norml2();
         n /= n_norm;

         assert(1. - n.Norml2() < 1e-12);

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
            // cout << "\t---n:\n";
            // n.Print(cout);
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

/****************************************************************************************************
* Function: GetEntityDof
* Parameters:
*     GDof   - Global degree of freedom with a max of H1.GetNDofs()
*     Entity - enum representing what type of node GDof corresponds to
*     EDof   - Entity index of the GDof
*
* Purpose:
*     This function is a helper function to identify which entity a node corresponds to
*     and to convert from the global numbering to the entity's numbering.
****************************************************************************************************/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::GetEntityDof(const int GDof, DofEntity & entity, int & EDof)
{
   // cout << "================== \n GetEntityDof \n================== \n";
   
   switch (dim)
   {
      case 1:
      {
         if (GDof < num_vertices)
         {
            // In 1d a face is a vertex, however our method to compute RT velocities
            // should return the intermediate face velocity.  Thus we shall treat
            // them as faces.
            entity = DofEntity::face;
            EDof = GDof;
         }
         else if (GDof < Vsize_H1) 
         {
            entity = DofEntity::cell;
            EDof = GDof - num_vertices;
         }
         break;
      }
      case 2:
      {
         if (GDof < num_vertices)
         {
            entity = DofEntity::corner;
            EDof = GDof;
         }
         else if (GDof < num_vertices + num_faces)
         {
            entity = DofEntity::face;
            EDof = GDof - num_vertices;
         }
         // TODO: Add edge implementation
         else
         {
            entity = DofEntity::cell;
            EDof = GDof - num_vertices - num_faces;
         }
         break;
      }
      case 3:
      {
         MFEM_ABORT("3D not implemented yet.");
         break;
      }
      default:
      {
         MFEM_ABORT("Incorrect dim provided.");
      }
   }
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
   
   for (; ci < NDofs_L2; ci++) // cell iterator
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


/****************************************************************************************************
* Function: MakeTimeStep
* Parameters:
*     S  - BlockVector that stores mesh information, mesh velocity, and state variables.
*     t  - Current time
*     dt - Current timestep
*
* Purpose:
*     This function first computes the new state variables at the (n+1)th timestep, then
*     computes the mesh motion and moves the mesh accordingly.  Finally, we check mass
*     conservation to ensure this is preserved locally.
*     
****************************************************************************************************/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::MakeTimeStep(Vector &S, double & t, double dt)
{
   cout << "Making timestep\n";
   Vector S_new = S; // We need a place to store updated cell values.

   // Retrieve data from Monolithic block vector S
   Vector* sptr = const_cast<Vector*>(&S_new);
   ParGridFunction x_gf, mv_gf_new;
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);
   mv_gf_new.MakeRef(&H1, *sptr, block_offsets[1]);

   // Update state variables contained in S_new
   ComputeStateUpdate(S_new, t, dt);

   // Move the mesh
   if (mm)
   {
      MoveMesh(S, x_gf, mv_gf_new, t, dt);
   }

   // Update the paramters passed by reference
   S = S_new;
   t += dt;

   // Verify that our algorithm is locally mass conservative
   CheckMassConservation(S);
} 


/****************************************************************************************************
* Function: ComputeStateUpdate
* Parameters:
*     S_new - BlockVector that stores mesh information, mesh velocity, and state variables.
*     t     - Current time
*     dt    - Current timestep
*
* Purpose:
*     This function takes in a BlockVector by reference and updates the state variables. 
*     If the option 'use_viscosity' is used at runtime, then the invariant domain preserving
*     method as described in (Guermond et al, eqn 4.9) is used.
****************************************************************************************************/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::ComputeStateUpdate(Vector &S_new, const double &t, const double dt)
{
   // cout << "========================================\n";
   // cout << "Computing State Update\n";
   // cout << "========================================\n";

   Vector val(dim+2), c(dim), n(dim), n_int(dim), U_i(dim+2), U_j(dim+2), sums(dim+2);
   Array<int> fids, fids2, oris, oris2, verts; // TODO: Remove fids2, oris2, These are for testing.

   int ci = 0, cj = 0;
   double d, c_norm;

   mfem::Mesh::FaceInformation FI;
   H1.ExchangeFaceNbrData();

   for (; ci < NDofs_L2; ci++) // Cell iterator
   {
      GetCellStateVector(S_new, ci, U_i);
      
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
         CalcOutwardNormalInt(S_new, ci, fids[j], n_int);
         c = n_int;
         c /= 2.;
         c_norm = c.Norml2();
         n = n_int;
         double F = n.Norml2();
         n /= F;
         assert(1. - n.Norml2() < 1e-12);
         FI = pmesh->GetFaceInformation(fids[j]);

         if (FI.IsInterior())
         {
            // Get index information/state vector for second cell
            if (ci == FI.element[0].index) { 
               cj = FI.element[1].index; 
            }
            else { 
               assert(ci == FI.element[1].index);
               cj = FI.element[0].index; 
            }
            GetCellStateVector(S_new, cj, U_j); 

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
}


template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::GetCellStateVector(const Vector &S, 
                                              const int cell, 
                                              Vector &U)
{
   U.SetSize(dim + 2);

   // This happens when a face is a boundary face, one of the cells index is -1q
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
   Array<int> dofs;
   Array<int> sub_dofs;
   dofs.Append(cell);
   sub_dofs.Append(1);

   // Get current state grid functions
   Vector* sptr = const_cast<Vector*>(&S_new);
   ParGridFunction sv_gf, v_gf, ste_gf;
   sv_gf.MakeRef(&L2, *sptr, block_offsets[2]);
   v_gf.MakeRef(&L2V, *sptr, block_offsets[3]);
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
/****************************************************************************************************
* Function: CalcOutwardNormalInt
* Parameters:
*     S    - BlockVector corresponding to nth timestep that stores mesh information, 
*            mesh velocity, and state variables.
*     cell - Inside integer cell index
*     face - Face index to compute normal on
*     res  - Resulting vector passed by reference and modified.
*
* Purpose:
*     This function returns the integral over the given face of the normal pointing away from
*     the cell c and storing this value into the vector res.
*
*        /\                                 /\
*        |                                  |
*        | n ds   = [orthogonal(s) / ||s||] | 1 ds   = orthogonal(s).
*        |                                  |
*       \/ (face)                          \/ (face)
****************************************************************************************************/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::CalcOutwardNormalInt(const Vector &S, const int cell, const int face, Vector & res)
{
   res = 0.;

   mfem::Mesh::FaceInformation FI;
   Array<int> row;
   FI = pmesh->GetFaceInformation(face);

   int c, cp, face_dof;
   c = FI.element[0].index;
   cp = FI.element[1].index;

   switch(dim)
   {
      case 1:
      {
         Vector cell_center_x(dim), face_x(dim);
         face_dof = face;

         vertex_element->GetRow(face, row);

         int cell_gdof = cell + num_faces;
         get_node_position(S, cell_gdof, cell_center_x);
         get_node_position(S, face, face_x);

         subtract(face_x, cell_center_x, res);
         res *= 1. / res.Norml2();
         break;
      }
      case 2:
      {
         H1.GetFaceDofs(face, row);
         int face_vdof1 = row[1], face_vdof2 = row[0];
         face_dof = row[2];

         Vector face_x(dim), vdof1_x(dim), vdof2_x(dim), vdof1_v(dim), vdof2_v(dim);
         get_node_position(S, face_dof, face_x);
         get_node_position(S, face_vdof1, vdof1_x);
         get_node_position(S, face_vdof2, vdof2_x);
         cout << "caloutwardnormalint, cell: " << cell << ", face: " << face << endl;
         cout << "face_x: ";
         face_x.Print(cout);
         cout << "vdof1_x: ";
         vdof1_x.Print(cout);
         cout << "vdof2_x: ";
         vdof2_x.Print(cout);

         Vector secant(dim);
         subtract(vdof2_x, vdof1_x, secant);

         res = secant;
         Orthogonal(res);

         /* Ensure orientation of normal */
         if (FI.IsInterior()) // Interior face
         {
            // Orientation of the normal vector depends on the indices of
            // the cells that share that face.  The normal points towards
            // the cell with the lower global index.  This normal must be
            // flipped if the cell we are working with is the one of the
            // lower index.
            if (cell == max(c, cp))
            {
               // cout << "flipping normal\n";
               res *= -1.;
            }
         }
         cout << "int n: ";
         res.Print(cout);

         break;
      }
      case 3:
      {
         MFEM_ABORT("3D not implemented.\n");
      }
      default:
      {
         MFEM_ABORT("Invalid dim provided.\n");
      }
   }
}


template<int dim, int problem>
double LagrangianLOOperator<dim, problem>::cosineSimilarity(const Vector &v1, const Vector &v2)
{
   double val = 0.;

   // Only compute cosine similarity if we have nonzero vectors
   if (v1.Norml2() != 0 && v2.Norml2() != 0)
   {
      val = (v1 * v2) / (v1.Norml2() * v2.Norml2());
   }

   return val;
}

/*
* Function: Orthogonal
*
* Parameters:
*     v - Vector to be rotated
*
* Purpose:
*     Rotate a vector counterclockwise of angle pi/2.
*
* Example:
*  (0,-1) ---> (1,0)
*/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::Orthogonal(Vector &v)
{
   if (dim == 2)
   {
      double x = v(0), y = v(1);
      v(0) = -1. * y;
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
//          v_gf.MakeRef(&L2V, *sptr, block_offsets[3]);

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

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::
   compute_intermediate_face_velocities(const Vector &S,
                                        const double t,
                                        const string flag, // Default NA
                                        void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
   /***
    * Compute intermediate face velocities according to (5.7) 
    ***/
   mfem::Mesh::FaceInformation FI;
   int c, cp;
   Vector Uc(dim+2), Ucp(dim+2), n_int(dim), c_vec(dim), n_vec(dim), Vf(dim); 
   double d, c_norm, F;

   for (int face = 0; face < num_faces; face++) // face iterator
   {
      Vf = 0.;
      FI = pmesh->GetFaceInformation(face);
      c = FI.element[0].index;
      cp = FI.element[1].index;

      GetCellStateVector(S, c, Uc);

      if (flag == "NA")
      {
         if (FI.IsInterior())
         {
            // this vector is only needed for interior faces
            GetCellStateVector(S, cp, Ucp);

            // Get normal, d, and |F|
            CalcOutwardNormalInt(S, c, face, n_int);
            n_vec = n_int;

            double F = n_vec.Norml2();
            n_vec /= F;
            assert(1. - n_vec.Norml2() < 1e-12);
            c_vec = n_int;
            c_vec /= 2.;
            double c_norm = c_vec.Norml2();

            // cout << "(mm)\tn:\n";
            // n_vec.Print(cout);
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
      }
      else 
      {
         assert(flag=="testing");
         Vector face_x(dim);
         int face_dof;

         // dof ordering depends on dimension
         switch (dim)
         {
            case 1:
            {
               face_dof = face;
               break;
            }
            case 2:
            {
               Array<int> row;
               H1.GetFaceDofs(face, row);
               face_dof = row[2];               
               break;
            }
            case 3:
            {
               MFEM_ABORT("3D not implemented.\n");
            }
            default:
            {
               MFEM_ABORT("Invalid dimension value.\n");
            }
         }

         get_node_position(S, face_dof, face_x);
         test_vel(face_x, t, Vf);
      }

      // Put face velocity into object
      for (int i = 0; i < dim; i++)
      {
         int index = face + i * num_faces;
         v_CR_gf[index] = Vf[i];
      }

   } // End face iterator

}


template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::get_intermediate_face_velocity(const int & face, Vector & vel)
{
   assert(face < num_faces);
   // Retrieve face velocity from object
   for (int i = 0; i < dim; i++)
   {
      int index = face + i * num_faces;
      vel[i] = v_CR_gf[index];
   }
}

/****************************************************************************************************
* Function: MoveMesh
* Parameters:
*     S         - BlockVector corresponding to nth timestep that stores mesh information, 
*                 mesh velocity, and state variables.
*     x_gf      - ParGridFunction representing the mesh that references BlockVector at 
                  (n+1)th timestep
*     mv_gf_new - ParGridFunction representing the mesh velocity that references 
                  BlockVector at (n+1)th timestep
*     t         - Current time
*     dt        - Current timestep
*
* Purpose:
*     This function accomplishes the following steps:
*        1) Constructs intermediate face velocity objects for all faces at the nth timestep.
*        2) Computes corrected node velocities on the corner vertices.
*        3) Fills unnecessary node at cell center with the average velocity of the corner nodes.
*        4) Computes the corrective face velocities.
*        5) Modify x_gf to represent the moved nodes (Computes the (n+1)th mesh location).
****************************************************************************************************/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::MoveMesh(Vector &S, GridFunction &x_gf, GridFunction & mv_gf_new,const double & t, const double & dt)
{
   // cout << "Moving mesh\n";
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction mv_gf_old;
   mv_gf_old.MakeRef(&H1, *sptr, block_offsets[1]);

   switch (problem)
   {
      case 4: // Noh problem
      {
         VectorFunctionCoefficient velocity_coeff(dim, InitialValues<problem, dim>::v0);
         velocity_coeff.SetTime(t);
         mv_gf_old.ProjectCoefficient(velocity_coeff);

         break;
      }
      default: // Move mesh according to paper
      {
         compute_intermediate_face_velocities(S, t);
         // Construct ti, face normals
         // ti = [0. 0.]^T in 2D
         // 3DTODO: Modify this according to seciton 5.2

         compute_node_velocities(S, t, dt);
         if (dim > 1)
         {
            // Don't need to fill face velocities with average in dim=1
            //3DTODO: Modify for 3d case
            fill_face_velocities_with_average(S);
            // compute_corrective_face_velocities(S, t, dt);
         }
         fill_center_velocities_with_average(S);
      }
      break;
   }

   // Finally, move x_gf according to mv_gf_old
   add(x_gf, dt, mv_gf_old, x_gf);
   pmesh->NewNodes(x_gf, false);

   // Set mv_gf that corresponds to S_new.
   mv_gf_new = mv_gf_old;
}


template<int dim, int problem>
double LagrangianLOOperator<dim, problem>::CalcMassLoss(const Vector &S)
{
   // cout << "Checking Mass Conservation\n";
   Vector U_i(dim + 2);

   double num = 0., denom = 0.;

   for (int ci = 0; ci < NDofs_L2; ci++)
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
void LagrangianLOOperator<dim, problem>::CheckMassConservation(const Vector &S)
{
   Vector U_i(dim + 2);
   int counter = 0;
   
   for (int ci = 0; ci < NDofs_L2; ci++)
   {
      const double m = m_hpv->Elem(ci);
      const double k = pmesh->GetElementVolume(ci);
      GetCellStateVector(S, ci, U_i);

      if (abs(k / U_i[0] - m) > pow(10, -8))
      {
         counter++;
         // cout << "cell: " << ci << endl;
         // cout << "MASS CONSERVATION BROKEN!!\n";
         // cout << "k / U_i[0] - m = " << k / U_i[0] - m << endl;
         // cout << "K: " << k << endl;
         // cout << "T: " << U_i[0] << endl;
         // cout << "m: " << m << endl;
         // cout << "K/T: " << k / U_i[0] << endl;
         // cout << endl;
      }
   }
   double cell_ratio = counter / (double)NDofs_L2;
   if (cell_ratio > 0.)
   {
      cout << "Percentage of cells where mass conservation was broken: " << cell_ratio << endl;
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
}


/* * * * * Corrected Node velocity calculations * * * * */
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::compute_A(const DenseMatrix & C, const double d, const double &dt, DenseMatrix &A)
{

   DenseMatrix _mat(C);
   _mat *= -dt;

   _mat(0,0) += d;
   _mat(1,1) += d;

   _mat.Invert();

   Mult(_mat, C, A);
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::compute_B(const DenseMatrix &C, const Vector & D, const double d, const double &dt, Vector &B)
{
   DenseMatrix _mat(C);
   _mat *= -dt;
   _mat(0,0) += d;
   _mat(1,1) += d;
   _mat.Invert();
   _mat.Mult(D, B);
}

// Assumes dt, not dt/2 is passed in
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::compute_determinant(const DenseMatrix &C, const double &dt, double & d)
{
   double a = 1.;
   double b = -1. * (1. + (dt / 2.) * (C(0,0) + C(1,1)));
   double c = (pow(dt, 2) / 4.) * (C(0,1) * C(1,0) + C(0,0) * C(1,1));

   cout << "a: " << a << ", b: " << b << ", c: " << c << endl;

   cout << "discriminant: " << pow(b,2) - 4. * a * c << endl;
   double pm = sqrt(pow(b,2) - 4. * a * c);
   cout << "pm: " << pm << endl;

   double d1 = -1 * b + pm;
   d1 /= (2. * a);

   double d2 = -1 * b - pm;
   d2 /= (2. * a);

   cout << "d1: " << d1 << endl;
   cout << "d2: " << d2 << endl;

   d = std::max(d1, d2);
}

// template<int dim, int problem>
// void LagrangianLOOperator<dim, problem>::compute_determinant(const DenseMatrix &C, const double &dt, double & d)
// {
//    double trace = C.Trace();
//    double det = C.Det();
//    // cout << "C:\n";
//    // C.Print(cout);
//    // cout << "trace: " << trace << ", det: " << det << ", dt: " << dt << endl;

//    double num = 0., denom = 0., pm = 0.;
//    double d1 = 0., d2 = 0.;

//    num = dt * trace + 1.;
//    denom = 2.;
//    pm = sqrt(pow(dt, 2) * pow(trace, 2) + 2. * dt * trace + 1. - 4. * pow(dt,2) * det);
//    d1 = (num + pm) / denom;
//    d2 = (num - pm) / denom;

//    if (isnan(d1) || isnan(d2))
//    {
//       MFEM_ABORT("NaN values returned in compute_determinant.\n");
//    }

//    // Will need to change
//    // How to pick which d?
//    // Check invertability of (d I - dt C)?
//    // One of the values is 0. Hmmm.
//    if (d1 > 0.)
//    {
//       d = d1;
//    }
//    else
//    {
//       assert(d2 > 0.);
//       d = d2;
//    }
// }

/*
Function: compute_corrected_node_velocity
Parameters:
   C, D     - Matrix and vector representing linear velocity field at node
   dt       - timestep
   vertex_x - Vector coordinates for node
   vertex_v - Computed velocity
   flag     - flag used to indicate testing
   test_vel - velocity used for testing
Purpose:
   This function takes in the velocity field that is computed at a node 
   using Least Squares or the Rannacher-Turek FE, and returns the corrected
   velocity that preserves linear motion (i.e. face nodes don't have to correct
   when the velocity is linear).
*/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::
   compute_corrected_node_velocity(const DenseMatrix &C, 
                                   const Vector &D, 
                                   const double & dt, 
                                   const Vector &vertex_x,
                                   Vector &vertex_v, 
                                   const string flag, 
                                   void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
   cout << "========================================\n";
   cout << "Computing corrected node_velocity\n";
   cout << "========================================\n";
   // V_exact = Cx + D
   // Will need access to what C and D are from the computed velocity from LS
   double d = 0.;
   DenseMatrix A(dim, dim);
   Vector B(dim);

   // First, compute d
   compute_determinant(C, dt/2., d);
   cout << "d: " << d << endl;

   // Second, compute A
   compute_A(C, d, dt/2., A);
   // cout << "A:\n";
   // A.Print(cout);

   // Third, compute B
   compute_B(C, D, d, dt/2., B);
   // cout << "B:\n";
   // B.Print(cout);

   /* Print values */
   cout << "- LS Velocity -\n";
   cout << "C:\n";
   C.Print(cout);
   cout << "D:\n";
   D.Print(cout);

   cout << "- Corrected Velocity -\n";
   cout << "A:\n";
   A.Print(cout);
   cout << "B:\n";
   B.Print(cout);

   cout << "vertex_x:\n";
   vertex_x.Print(cout);
   cout << "Pre corrected velocity:\n";
   vertex_v.Print(cout);
   // cout << "C:\n";
   // C.Print(cout);
   // cout << "A:\n";
   // A.Print(cout);
   // cout <<"D:\n";
   // D.Print(cout);
   // cout << "B:\n";
   // B.Print(cout);

   // Finally, evaluate velocity at vertex_x and put resulting value into vertex_v
   A.Mult(vertex_x, vertex_v);
   vertex_v += B;
   cout << "Post corrected velocity:\n";
   vertex_v.Print(cout);
}

/*
Function: compute_node_velocity_LS
Parameters:
   S        - BlockVector representing FiniteElement information
   t        - Current time
   dt       - Timestep
   flag     - Flag used to indicate testing
   test_vel - Velocity used for testing
Purpose:
   This function computes the nodal velocity of all nodes in the mesh.  It does
   this using the intermediate face velocities from the adjacent faces, forming 
   an overdetermined system for the local velocity field under the assumption that
   it is a linear field, and solves the overdetermined system using Least Squares.
*/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::
   compute_node_velocity_LS(const Vector &S,
                            const Table &vertex_edge, 
                            const int & vertex,
                            const double &t,
                            const double &dt,
                            Vector &vertex_v,
                            DenseMatrix &C,
                            Vector &D,
                            const string flag, // Default NA
                            void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
   cout << "========================================\n";
   cout << "Computing node velocity using Least Squares.\n";
   cout << "========================================\n";

   /* Code from previous LS */
   bool is_boundary_node = false;
   mfem::Mesh::FaceInformation FI;
   Array<int> row;
   Array<int> face_dof_row;
   int row_length, face, face_dof;

   /* Get row information */
   vertex_edge.GetRow(vertex, row);
   row_length = row.Size();

   DenseMatrix A(row_length, dim + 1), At(dim+1, row_length), LHS(dim+1, dim+1);
   Vector Bx(row_length), By(row_length), vertex_x(dim);
   Vector Bx_mod(dim+1), By_mod(dim+1); // To solve LS problem 
   Vector Vf(dim); // Intermediate Face Velocity or exact velocity if testing
   Vector n_vec(dim); // For enforcing BCs 
   Vector face_x(dim);
   
   // Get position of node
   get_node_position(S, vertex, vertex_x);
   
   // 3DTODO: Will need to modify this for faces instead of edges in 3D
   for (int face_it = 0; face_it < row_length; face_it++) // Adjacent face iterator
   {
      // Reset face velocity object
      Vf = 0.;

      face = row[face_it];
      FI = pmesh->GetFaceInformation(face);
      if (FI.IsBoundary()) { is_boundary_node = true; }

      // Retrieve coordinates for face node
      H1.GetFaceDofs(face, face_dof_row);
      face_dof = face_dof_row[2]; // Dof corresponding to the center of the face
      get_node_position(S, face_dof, face_x);

      // Fill corresponding row of A, Bx, By
      for (int i = 0; i < dim; i++)
      {
         A(face_it, i) = face_x(i) - vertex_x(i); // Fill with nodal locations
      }
      A(face_it, dim) = 1.;

      get_intermediate_face_velocity(face, Vf);

      // In either case, we set the corresponding row to the face velocity,
      // whether that is given by the intermediate face velocity or the exact
      Bx(face_it) = Vf(0);
      By(face_it) = Vf(1);
   } // End face iterator

   // Form system
   At.Transpose(A);
   Mult(At, A, LHS);
   At.Mult(Bx, Bx_mod);
   At.Mult(By, By_mod);

   // Solve linear systems A x_A = B_x, A x_B = B_y
   Vector res_x(dim+1), res_y(dim+1);
   
   /*
   ** New stuff
   **
   ** */

   // Must handle underdetermined system
   if (row_length < dim + 1)
   {
      if (flag=="testing") { 
         test_vel(vertex_x, t, vertex_v); 

         // Get C and D from the test_vel function
         Vector temp_vec(dim), temp_res(dim);
         temp_vec = 0.;
         test_vel(temp_vec, t, temp_res);
         D(0) = temp_res(0);
         D(1) = temp_res(1);
         temp_vec(0) = 1.;
         test_vel(temp_vec, t, temp_res);
         C(0,0) = temp_res(0) - D(0); // a
         C(1,0) = temp_res(1) - D(1); // d

         temp_vec = 0.;
         temp_vec(1) = 1.;
         test_vel(temp_vec, t, temp_res);
         C(0,1) = temp_res(0) - D(0); // a
         C(1,1) = temp_res(1) - D(1); // d

         cout << "testing C:\n";
         C.Print(cout);
         cout << "testing D:\n";
         D.Print(cout);
      }
      else 
      {
         CGSolver cg;
         cg.SetOperator(LHS);
         cg.SetRelTol(1e-12);
         cg.SetAbsTol(1e-12);
         cg.SetMaxIter(1000);
         cg.SetPrintLevel(-1);
         cg.Mult(Bx_mod, res_x);
         cg.Mult(By_mod, res_y);

         // Solve for nodal velocity
         vertex_v(0) = res_x(dim); // Set to cx
         vertex_v(1) = res_y(dim); // Set to cy

         // Form C, D
         C(0, 0) = res_x(0);
         C(0, 1) = res_x(1);
         C(1, 0) = res_y(0);
         C(1, 1) = res_y(1);

         D(0) = res_x(dim) - (res_x(0) * vertex_x(0) + res_x(1) * vertex_x(1));
         D(1) = res_y(dim) - (res_y(0) * vertex_x(0) + res_y(1) * vertex_x(1));
      } 
   }
   else // determined or overdetermined system
   { 
      LHS.Invert();
      LHS.Mult(Bx_mod, res_x);
      LHS.Mult(By_mod, res_y);

      // Solve for nodal velocity
      vertex_v(0) = res_x(dim); // Set to cx
      vertex_v(1) = res_y(dim); // Set to cy

      // TODO: We will need to ensure the timestep restriction
      // Form C, D
      C(0, 0) = res_x(0);
      C(0, 1) = res_x(1);
      C(1, 0) = res_y(0);
      C(1, 1) = res_y(1);

      D(0) = res_x(dim) - (res_x(0) * vertex_x(0) + res_x(1) * vertex_x(1));
      D(1) = res_y(dim) - (res_y(0) * vertex_x(0) + res_y(1) * vertex_x(1));

      // cout << "determined or overdetermined system:\n";
      // cout << "testing C:\n";
      // C.Print(cout);
      // cout << "testing D:\n";
      // D.Print(cout);
   }
   
   // Enforce boundary conditions
   if (is_boundary_node)
   {
      /* Enforce BCs (Copy/Paste from previous nodal computation) */
      switch (problem)
      {
         // case 2:
         // case 3:
         // {
         //    vertex_v = 0.;
         //    break;
         // }
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
                     double _xi = t / (2.*timestep_first);

                     double _psi = (4. - (_xi + 1.) * (_xi - 2.) * ((_xi - 2.) - (abs(_xi-2.) + (_xi-2.)) / 2.)) / 4.;

                     vertex_v[0] = 1. * _psi;

                     break;
                  }
                  else 
                  {
                     assert (bdr_attribute == 2);
                     CalcOutwardNormalInt(S, FI.element[0].index, face, n_vec);
                     n_vec /= n_vec.Norml2();
                     assert(1. - n_vec.Norml2() < 1e-12);

                     Vector vertex_v_normal_component(dim);
                     vertex_v_normal_component = n_vec;
                     double _scale = vertex_v * n_vec;
                     vertex_v_normal_component *= _scale;

                     vertex_v -= vertex_v_normal_component;
                  }
               } // End boundary face
            } // End: Face iterator
            break;
         } // End: Saltzman
         default:
         {
            // No change
         }
      } // Switch case end
   }
}

/*
Function: compute_node_velocity_RT
Parameters:
   dt       - Timestep
Purpose:

*/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::
   compute_node_velocity_RT(const int & node, const double & dt, Vector &node_v)
{
   switch (dim)
   {
      case 1:
      {
         // Remark 5.3 in the paper states that since each geometric node belongs to one
         // face only, the RT velocity is equal to the IFV previously computed.
         assert(node < num_faces);
         get_intermediate_face_velocity(node, node_v);
         break;
      }
      default:
      {
         DenseMatrix Ci(dim);
         Vector Vgeo(dim);
         double d = 0.;
         
         Ci = 0., Vgeo = 0., node_v = 0.;
         compute_geo_C(node, Ci);
         compute_geo_V(node, Vgeo);

         compute_determinant(Ci, dt, d);

         // Compute V_i^n
         DenseMatrix _mat(Ci);
         _mat *= - dt / 2.;
         for (int i = 0; i < dim; i++)
         {
            _mat(i,i) += d;
         }
         _mat.Invert();
         _mat.Mult(Vgeo, node_v);
         break;
      }
   }
}


/***********************************************************************************************************
Function: RT_nodal_velocity
Parameters:
   cell - index corrseponding to the cell (K_c)
   node - global index of the node to calculate the velocity on (node < H1.GetNDofs)
   vel  - returned velocity v_h^v_K(x_i)
Purpose:
   The purpose of this function is to compute the Rannacher-Turek constructed velocity of a given cell
   at a given node in the mesh.  If that node is not contained in the cell, this function will throw
   an error.

   NOTE: This function assumes dim > 1.
   NOTE: This function assumes that the function LagrangianLOOperator:compute_intermediate_face_velocities()
   has already been called.  If this function has not been called, then the returned velocity will be 0.
***********************************************************************************************************/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::RT_nodal_velocity(const int & cell, const int & node, Vector &vel)
{
   assert(node < NDofs_H1); // "Invalid nodal index"
   assert(cell < NDofs_L2); // "Invalid cell index"

   DofEntity entity;
   int EDof;
   GetEntityDof(node, entity, EDof);

   // Get DoFs corresponding to cell
   Array<int> element_face_row, element_face_oris;
   switch (dim)
   {
      case 1:
      {
         pmesh->GetElementVertices(cell, element_face_row);
         // In 1D, "face" refers to vertex
         break;
      }
      case 2: 
      {
         pmesh->GetElementEdges(cell, element_face_row, element_face_oris);
         // In 2D, "face" refers to edge
         break;
      }
      case 3:
      {
         MFEM_ABORT("3D not implemented.\n");
         // In 3D, "face" refers to face
         // pmesh->GetElementFaces(cell, element_face_row, element_face_oris);
         break;
      }
      default:
      {
         MFEM_ABORT("Improper value provided.\n");
      }
   }
   // element_face.GetRow(cell, element_face_row);
   
   int row_length = element_face_row.Size();

   // cout << "element faces for cell " << cell << ": ";
   // element_face_row.Print(cout);
   // cout << "looking for EDof: " << EDof << endl;

   // Reset velocity
   vel = 0.;


   switch (entity)
   {
      case 0: // corner
      {
         // Get integration point corresponding to node location
         Array<int> verts;
         pmesh->GetElementVertices(cell, verts);

         // cout << "cell: " << cell << ", verts: \n";
         // verts.Print(cout);

         IntegrationPoint test_ip;
         if (node == verts[0]) {
            test_ip.Set2(0., 0.);
         } else if (node == verts[1]) {
            test_ip.Set2(1., 0.);
         } else if (node == verts[2]) {
            test_ip.Set2(1., 1.);
         } else if (node == verts[3]) {
            test_ip.Set2(0., 1.);
         } else { 
            // Invalid if node is not contained in cell
            MFEM_ABORT("Incorrect node provided.\n"); 
         }

         // Evaluate reference shape functions at integration point
         Vector shapes(4);
         const FiniteElement * fe = CR.GetFE(cell);
         fe->CalcShape(test_ip, shapes);

         // Sum over faces evaluating local velocity at vertex
         Vector Ci(dim);
         Ci = 0.;
         Vector face_vel(dim);

         for (int face_it = 0; face_it < row_length; face_it++)
         {
            int face = element_face_row[face_it];

            // Get intermediate face velocities corresponding to faces
            get_intermediate_face_velocity(face, face_vel);
            vel.Add(shapes[face_it], face_vel);
         }
         break;
      }
      case 1: // face
      {
         bool face_is_part_of_cell = false;
         // Check is face is part of the provided cell
         for (int face_it = 0; face_it < row_length; face_it++)
         {
            if (EDof == element_face_row[face_it])
            {
               face_is_part_of_cell = true;
            }
         }

         if (!face_is_part_of_cell)
         {
            cout << "cell: " << cell << ", node " << node << endl;
            MFEM_ABORT("Provided face is not adjacent to cell");
         }

         // Simply return the intermediate face velocity
         get_intermediate_face_velocity(EDof, vel);
         break;
      }
      case 2: // cell
      {
         break;
      }
      default:
      {
         MFEM_ABORT("Invalid entity value.\n");
      }
   }
}

/***********************************************************************************************************
Function: RT_int_grad
Parameters:
   cell - index corrseponding to the cell (K_c)
   res  - (dim x dim) DenseMatrix representing the outer product of Vfm with 
          the gradient of its corresponding scalar RT shape function.
Purpose:
   This function calculates the integral of the gradient of the RT velocity function
   on a given cell.

   NOTE: This function assumes that dim > 1.
   NOTE: This function assumes that the function LagrangianLOOperator:compute_intermediate_face_velocities()
   has already been called.  If this function has not been called, then the returned velocity will be 0.
***********************************************************************************************************/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::RT_int_grad(const int cell, DenseMatrix & res)
{
   // cout << "RT_int_grad funcall.\n";
   ParGridFunction CRc_gf(&CRc);
   const int size = CRc.GetVSize();

   ElementTransformation * trans = pmesh->GetElementTransformation(cell);
   Vector grad(dim), row(dim);
   res = 0., row = 0.;
   res.SetSize(dim);

   // Iterate over quadrature
   for (int i = 0; i < RT_ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = RT_ir.IntPoint(i);
      trans->SetIntPoint(&ip);
      // cout << "ip.x: " << ip.x << ", ip.y: " << ip.y << ", weight: " << ip.weight << endl;

      // cout << "el " << cell << " at integration point " << i << endl;
      for (int j = 0; j < dim; j++)
      {
         // cout << "dim: " << j << endl;
         CRc_gf.MakeRef(&CRc, v_CR_gf, j*size);
         CRc_gf.GetGradient(*trans, grad);
         // cout << "grad: ";
         // grad.Print(cout);

         // Put information into Dense Matrix
         res.GetRow(j, row);
         row.Add(ip.weight, grad);
         res.SetRow(j, row);
      }
   }
}

/***********************************************************************************************************
Function: compute_geo_V
Parameters:
   node - Index corresponding to the global node index, to include face nodes and cell center nodes.  
          Note that we should not be computing any sort of averaged velocity at the cell centers.
   res  - Averaged Vector resulting from the average value of the cell wise velocities of all
          adjacent cells at the node
Purpose:
   Part of the process to reconstructing the continuous geometric velocity field from the discontinuous
   RT representation.  The exact equation to be solved is given by equation (5.11).

   NOTE: This function assumes dim > 1.

   NOTE: This function assumes that the function LagrangianLOOperator:compute_intermediate_face_velocities()
   has already been called.  If this function has not been called, then the returned velocity will be 0.
***********************************************************************************************************/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::compute_geo_V(const int &node, Vector & res)
{
   assert(node < NDofs_H1); // "Invalid nodal index"
   res.SetSize(dim);
   res = 0.;

   Vector temp(dim);
   mfem::Array<int> row;

   DofEntity entity;
   int EDof;
   GetEntityDof(node, entity, EDof);

   switch (entity)
   {
      case 0: // corner
      {
         vertex_element->GetRow(EDof, row);
         break;
      }
      case 1: // face
      {
         face_element->GetRow(EDof, row);
         break;
      }
      case 2: // cell
      {
         MFEM_ABORT("No need to compute V_i at cell centers");
      }
      default:
      {
         MFEM_ABORT("Invalid entity\n");
      }
   }

   int row_length = row.Size();
   for (int row_it = 0; row_it < row_length; row_it++)
   {
      int row_el = row[row_it];

      temp = 0.;
      RT_nodal_velocity(row_el, node, temp);
      // cout << "velocity computed on cell " << row_el << " for node: " << node << endl;
      // temp.Print(cout);
      res.Add(1., temp);
   }

   // cout << "denom: " << denom << endl;
   res *= 1./row_length;

   // cout << "resulting velocity: ";
   // res.Print(cout);
}

/***********************************************************************************************************
Function: compute_geo_C
Parameters:
   node - Index corresponding to the global node index, to include face nodes and cell center nodes.  
          Note that we should not be computing any sort of averaged velocity at the cell centers.
   res  - Averaged DenseMatrix resulting from the average gradient of the cell wise velocities of all
          adjacent cells.
Purpose:
   Part of the process to reconstructing the continuous geometric velocity field from the discontinuous
   RT representation.  The exact equation to be solved is given by equation (5.11)

   NOTE: This function assumes dim > 1.
   NOTE: This function assumes that the function LagrangianLOOperator:compute_intermediate_face_velocities()
   has already been called.  If this function has not been called, then the returned velocity will be 0.
***********************************************************************************************************/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::compute_geo_C(const int &node, DenseMatrix & res)
{
   assert(node < NDofs_H1); // "Invalid nodal index"
   res.SetSize(dim);
   res = 0.;

   // Check if node corresponds to a face or to a vertex
   // and get adjacent cell indices
   double denom = 0.;
   DenseMatrix dm_temp(dim);
   mfem::Array<int> row;

   DofEntity entity;
   int EDof;
   GetEntityDof(node, entity, EDof);
   // cout << "entity: " << entity << ", Edof: " << EDof << endl;

   switch (entity)
   {
      case 0: // corner
      {
         // cout << "corner node\n";
         vertex_element->GetRow(node, row);
         break;
      }
      case 1: // face
      {
         // cout << "face node\n";
         face_element->GetRow(EDof, row);
         break;
      }
      case 2: // cell
      {
         MFEM_ABORT("No need to compute C_i at cell centers");
      }
      default:
      {
         MFEM_ABORT("Invalid entity\n");
      }
   }
   // cout << "vertex_element row: " << endl;
   // row.Print(cout);

   int row_length = row.Size();
   for (int row_it = 0; row_it < row_length; row_it++)
   {
      int row_el = row[row_it];
      denom += abs(pmesh->GetElementVolume(row_el));

      dm_temp = 0.;
      RT_int_grad(row_el, dm_temp);
      res.Add(1., dm_temp);
      // cout << "int_grad for row el: " << row_el << endl;
      // dm_temp.Print(cout);
   }

   // cout << "denom: " << denom << endl;
   res *= 1./denom;

   // cout << "Ci for node: " << node << endl;
   // res.Print(cout);
}


template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::
   compute_node_velocities(Vector &S, 
                           const double & t, 
                           const double & dt, 
                           const string flag, // Default NA
                           void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
   Vector vertex_v(dim);

   // Iterate over vertices
   for (int vertex = 0; vertex < num_vertices; vertex++) // Vertex iterator
   {
      compute_node_velocity_RT(vertex, dt, vertex_v);

      update_node_velocity(S, vertex, vertex_v);
   } // End Vertex iterator
}


/*
Function: compute_corrective_face_velocities
Parameters:
   S        - BlockVector representing FiniteElement information
   t        - Current time
   dt       - Timestep
   flag     - Flag used to indicate testing
   test_vel - Velocity used for testing
Purpose:
   This function computes the node velocities on the faces which
   are designed to bubble in the direction of the normal vector 
   to conserve mass locally.
*/
template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::
   compute_corrective_face_velocities(Vector &S, const double & t, const double & dt,
                                    const string flag, // Default NA
                                    void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
   assert(dim > 1); // No need to correct face velocities in dim=1
   // cout << "========================================\n";
   // cout << "Computing corrective interior face velocities.\n";
   // cout << "========================================\n";
   /* Parameters needed for face velocity calculations */
   mfem::Mesh::FaceInformation FI;
   Vector Vf(dim), n_int(dim), n_vec(dim), face_velocity(dim);
   Vector Uc(dim+2), Ucp(dim+2);
   Vector cell_c_v(dim), cell_cp_v(dim), cell_center_v(dim); // For computation of bmn
   Array<int> row;
   int face_dof = 0, face_vdof2 = 0, c = 0, cp = 0;

   // Iterate over faces
   for (int face = 0; face < num_faces; face++) // face iterator
   {  
      // Get intermediate face velocity, face information, and face normal
      face_velocity = 0.;
      Vf = 0.;

      FI = pmesh->GetFaceInformation(face);

      /* adjacent corner indices */
      H1.GetFaceDofs(face, row);
      int face_vdof1 = row[1], face_vdof2 = row[0], face_dof = row[2]; // preserve node orientation discussed in appendix A

      // retrieve old face and corner locations
      Vector face_x(dim), vdof1_x(dim), vdof2_x(dim), vdof1_v(dim), vdof2_v(dim);
      get_node_position(S, face_dof, face_x);
      get_node_position(S, face_vdof1, vdof1_x);
      get_node_position(S, face_vdof2, vdof2_x);

      // cout << "------------------------------------\n";
      // cout << "dt: " << dt << endl;
      // cout << "We are iterating on face: " << face << endl;
      // cout << "This face is located at: \n";
      // face_x.Print(cout);
      // cout << "The adjacent vertices are vertex 1: " << face_vdof1 << ", located at: \n";
      // vdof1_x.Print(cout);
      // cout << "and vertex 2: " << face_vdof2 << ", located at: \n";
      // vdof2_x.Print(cout);

      // retrieve corner velocities
      get_node_velocity(S, face_vdof1, vdof1_v);
      get_node_velocity(S, face_vdof2, vdof2_v);
      
      // cout << "The computed velocity at vertex: " << face_vdof1 << " is: \n";
      // vdof1_v.Print(cout);
      // cout << "The computed velocity at vertex: " << face_vdof2 << " is: \n";
      // vdof2_v.Print(cout);

      if (FI.IsInterior())
      {
         // Get adjacent cell information
         c = FI.element[0].index;
         cp = FI.element[1].index;
         GetCellStateVector(S, c, Uc);
         GetCellStateVector(S, cp, Ucp);
         cell_c_v = velocity(Uc);
         cell_cp_v = velocity(Ucp);
         
         // cout << "We have an interior face: " << face << endl;
         // cout << "Cell c: " << c << ", cell cp: " << cp << endl;
         // cout << "Velocity on cell c:\n";
         // cell_c_v.Print(cout);
         // cout << "Velocity on cell cp:\n";
         // cell_cp_v.Print(cout);

         // Calculate outer normal
         CalcOutwardNormalInt(S, c, face, n_int);
         n_vec = n_int;
         double F = n_vec.Norml2();
         n_vec /= F;

         assert(1. - n_vec.Norml2() < 1e-12);

         // cout << "outward normal integral: ";
         // n_int.Print(cout);
         // cout << "Face Length: " << F << endl;
         // cout << "The outward normal vector is: \n";
         // n_vec.Print(cout);

         // Calculate new corner locations and half  step locations
         Vector vdof1_x_new(dim), vdof2_x_new(dim), vdof1_x_half(dim), vdof2_x_half(dim);
         vdof1_x_new = vdof1_x;
         vdof1_x_new.Add(dt, vdof1_v);
         vdof1_x_half = vdof1_x;
         vdof1_x_half.Add(dt/2., vdof1_v);

         vdof2_x_new = vdof2_x;
         vdof2_x_new.Add(dt, vdof2_v);
         vdof2_x_half = vdof2_x;
         vdof2_x_half.Add(dt/2., vdof2_v);

         // calculate a_{12}^{n+1} (new tangent midpoint)
         Vector vdof12_x_new(dim);
         vdof12_x_new = vdof1_x_new;
         vdof12_x_new += vdof2_x_new;
         vdof12_x_new /= 2.;

         // Output new nodal location information
         // cout << "New locations for vdof1. full:\n";
         // vdof1_x_new.Print(cout);
         // cout << "half: \n";
         // vdof1_x_half.Print(cout);

         // cout << "New locations for vdof2. full:\n";
         // vdof2_x_new.Print(cout);
         // cout << "half: \n";
         // vdof2_x_half.Print(cout);

         // cout << "new tangent midpoint, vdof12_x_new: \n";
         // vdof12_x_new.Print(cout);

         // Compute D (A.4c)
         Vector n_vec_R(dim), temp_vec(dim), temp_vec_2(dim);
         n_vec_R = n_vec;
         Orthogonal(n_vec_R);
         // cout << "n_vec_R: \n";
         // n_vec_R.Print(cout);

         subtract(vdof1_v, vdof2_v, temp_vec); // V1 - V2 = temp_vec
         // cout << "V1-V2:\n";
         // temp_vec.Print(cout);

         subtract(vdof2_x_half, vdof1_x_half, temp_vec_2); // A2-A1
         // cout << "A2-A1:\n";
         // temp_vec_2.Print(cout);

         Orthogonal(temp_vec_2);
         // cout << "(A2-A1)^R:\n";
         // temp_vec_2.Print(cout);

         double D = dt * (temp_vec * n_vec_R) + 2. * (n_vec * temp_vec_2);
         // cout << "D: " << D << endl;

         // Compute c1 (A.4a)
         subtract(vdof2_v, vdof1_v, temp_vec); // only change temp_vec, since temp_vec_2 is same from D calculation (half step representation)
         // cout << "V2-V1:\n";
         // temp_vec.Print(cout);
         // cout << "c1 computation n_vec:\n";
         // n_vec.Print(cout);
         double c1 = ( dt * (temp_vec * n_vec) + 2. * (temp_vec_2 * n_vec_R) ) / D; // TRYING SOMETHING HERE. WILL NEED CORRECTED
         // cout << "c1 computation n_vec_R:\n";
         // n_vec_R.Print(cout);
         // cout << "c1: " << c1 << endl;

         // Compute c0 (A.4b)
         // cout << "!!!Computing c0\n";

         // cout << "\n === First we compute bmn ===\n";
         Vector n_vec_half(dim);
         subtract(vdof2_x_half, vdof1_x_half, n_vec_half);
         Orthogonal(n_vec_half);
         // cout << "n_vec_half:\n";
         // n_vec_half.Print(cout);
         // Vector v_exact_at_face(dim);
         // test_vel(face_x, 0., v_exact_at_face);
         // double bmn = v_exact_at_face * n_vec_half;
         
         get_intermediate_face_velocity(face, Vf);
         
         double bmn = Vf * n_vec;
         bmn *= F;

         // cout << "face length: " << face_length << endl;
         // bmn *= face_length;
         // cout << "bmn: " << bmn << endl;

         // cout << "\n\n==========================\n";
         // cout << "corrective face\nface: " << face << endl;
         // cout << "'inside' cell: " << c << endl;
         // cout << "Vf:\n";
         // Vf.Print(cout);
         // cout << "n_vec:\n";
         // n_vec.Print(cout);
         // cout << "=========================\n";

         temp_vec = vdof1_x_half;
         Orthogonal(temp_vec); // A1R
         // cout << "A1^R:\n";
         // temp_vec.Print(cout);
         temp_vec_2 = vdof2_x_half;
         Orthogonal(temp_vec_2); // A2R
         // cout << "A2^R:\n";
         // temp_vec_2.Print(cout);
         double const1 = vdof1_v * temp_vec - vdof2_v * temp_vec_2; // V1*A1R - V2*A2R
         double const2 = vdof1_v * temp_vec_2 - vdof2_v * temp_vec; // V1*A2R - V2*A1R
         
         temp_vec = face_x;
         Orthogonal(temp_vec);
         // cout << "a3^R:\n";
         // temp_vec.Print(cout);
         subtract(vdof2_v, vdof1_v, temp_vec_2);
         // cout << "V2 - V1:\n";
         // temp_vec_2.Print(cout);
         double const3 = temp_vec_2 * temp_vec; // (V2 - V1) * a3nR
         // cout << "const 1: " << const1 << ", const 2: " << const2 << ", const3: " << const3 << endl;
         double c2 = (3. / D) * (const1 / 2. + const2 / 6. + 2. * const3 / 3.);
         // cout << "c0: " << c0 << endl;
         
         // Compute V3n_perp
         // cout << "Computing V3n_perp\n";
         subtract(vdof2_x_new, vdof1_x_new, temp_vec);
         // cout << "a2n+1 - a1n+1:\n";
         // temp_vec.Print(cout);
         subtract(face_x, vdof12_x_new, temp_vec_2);
         temp_vec_2.Add((c2+ (3. * bmn / D))*dt, n_vec);
         // cout << "a3n - a12n+1 + c0dtn_vec:\n";
         // temp_vec_2.Print(cout);
         const1 = temp_vec * temp_vec_2; // numerator
         temp_vec_2 = n_vec_R;
         temp_vec_2 *= -1.;
         temp_vec_2.Add(c1, n_vec);
         // cout << "c1n_vec + nperp:\n";
         // temp_vec_2.Print(cout);
         const2 = temp_vec * temp_vec_2;
         const2 *= dt; // denominator
         double V3nperp = -1. * const1 / const2;
         // cout << "V3nperp: " << V3nperp << endl;

         // Compute V3n (5.11)
         double V3n = c1 * V3nperp + c2 + 3. * bmn / D;
         // cout << "V3n: " << V3n << endl;

         // Compute face velocity (5.11)
         // const1 = Vf * n_vec;
         // temp_vec = n_vec;
         // temp_vec *= const1;
         // subtract(Vf, temp_vec, face_velocity);
         // face_velocity.Add(V3n, n_vec);

         // Compute face velocity (Appendix A)
         face_velocity = 0.;
         face_velocity.Add(V3n, n_vec);
         Vector n_vec_perp(dim);
         n_vec_perp = n_vec_R;
         n_vec_perp *= -1.;
         // cout << "n_vec_perp:\n";
         // n_vec_perp.Print(cout);
         // cout << "n_vec_R:\n";
         // n_vec_R.Print(cout);
         face_velocity.Add(V3nperp, n_vec_perp);

         Vector face_x_new(dim);
         face_x_new = face_x;
         face_x_new.Add(dt, face_velocity);

         // Check perpendicular
         subtract(face_x_new, vdof12_x_new, temp_vec);
         subtract(vdof2_x_new, vdof1_x_new, temp_vec_2);
         if (abs(temp_vec * temp_vec_2) > pow(10, -12))
         {
            // cout << "temp_vec:\n";
            // temp_vec.Print(cout);
            // cout << "temp_vec_2:\n";
            // temp_vec_2.Print(cout);
            // cout << "################## Vectors are not orthogonal!";
            MFEM_ABORT("vectors are not orthogonal!\n");
         }

         // if (flag == "testing")
         // {
            // cout << "The coordinate location of this face is: \n";
            // face_x.Print(cout);
            // cout << "The computed velocity on this face is: \n";
            // face_velocity.Print(cout);
            // cout << "Since we shouldn't have any correction on the faces, the exact face velocity should be: \n";
            // Vector face_v_exact(dim);
            // test_vel(face_x, 0., face_v_exact);
            // face_v_exact.Print(cout);
            // Vector temp_vel(dim);
            // subtract(face_v_exact, face_velocity, temp_vel);

            // if (temp_vel.Norml2() > 10e-6)
            // {
            //    MFEM_ABORT("Incorrect face velocity.\n");
            // }
         // }
      } // End interior face

      // On the boundary, we just need the average of the adjacent 
      // vertices to account for the corrective corner velocities
      // NOTE that even in the test case, the corner velocities will
      // be modified and so prescribing the exact velocities on the
      // boundary faces will not preserve mass.
      else
      {
         // cout << "Boundary face~~\n";
         assert(FI.IsBoundary());
         // TODO: Add in boundary conditions similar to corner vertex bcs
         // TODO: Add testing validation step here
         for (int j = 0; j < dim; j++)
         {
            face_velocity[j] = (vdof1_v[j] + vdof2_v[j]) / 2.;
         }
        
      } // End boundary face

      // Lastly, put face velocity into gridfunction object
      update_node_velocity(S, face_dof, face_velocity);
   }

   // assert(false);
}


template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::
   fill_center_velocities_with_average(Vector &S,
                                       const string flag, 
                                       void (*test_vel)(const Vector&, const double&, Vector&))
{
   // Since we cannot use Serendipity elements, we must update cell center velocities
   Vector Vc(dim), Uc(dim+2), node_v(dim);

   Array<int> verts;
   Vector face_x(dim);
   for (int ci = 0; ci < NDofs_L2; ci++)
   {
      // Get center node dof
      int cell_vdof;
      switch (dim)
      {
         case 1:
         {
            cell_vdof = num_faces + ci;
            break;
         }
         case 2:
         {
            cell_vdof = NVDofs_H1 + num_faces + ci;
            break;
         }
         case 3:
         {
            MFEM_ABORT("3D not implemented.\n");
         }
         default:
         {
            MFEM_ABORT("Incorrect dim value provided.\n");
         }
      }

      get_node_position(S, cell_vdof, face_x);
      
      GetCellStateVector(S, ci, Uc);
      Vc = velocity(Uc);
      Vc = 0.;
      pmesh->GetElementVertices(ci, verts);
      for (int j = 0; j < verts.Size(); j++)
      {
         get_node_velocity(S, verts[j], node_v);
         Vc += node_v;
      }
      Vc /= verts.Size();

      // Update the velocity for the center node
      update_node_velocity(S, cell_vdof, Vc);
   }
}

template<int dim, int problem>
void LagrangianLOOperator<dim, problem>::fill_face_velocities_with_average(Vector &S, const string flag, void (*test_vel)(const Vector&, const double&, Vector&))
{
   /* Parameters needed for face velocity calculations */
   mfem::Mesh::FaceInformation FI;
   Vector Uc(dim+2), Ucp(dim+2);
   Vector cell_c_v(dim), cell_cp_v(dim); // For computation of bmn
   Array<int> row;
   int face_dof = 0, face_vdof1 = 0, face_vdof2 = 0;

   // Iterate over faces
   for (int face = 0; face < num_faces; face++) // face iterator
   {  
      // Get face information
      Vector face_velocity(dim);
      face_velocity = 0.;

      FI = pmesh->GetFaceInformation(face);

      /* adjacent corner indices */
      H1.GetFaceDofs(face, row);
      int face_vdof1 = row[1], face_vdof2 = row[0], face_dof = row[2]; // preserve node orientation discussed in appendix A

      // retrieve old face and corner locations
      Vector face_x(dim), vdof1_v(dim), vdof2_v(dim);

      // Get face location
      get_node_position(S, face_dof, face_x);

      // retrieve corner velocities
      get_node_velocity(S, face_vdof1, vdof1_v);
      get_node_velocity(S, face_vdof2, vdof2_v);
      

      // Average nodal velocities
      for (int j = 0; j < dim; j++)
      {
         face_velocity[j] = (vdof1_v[j] + vdof2_v[j]) / 2;
      }

      // test_vel(face_x, 0, face_velocity);

      // Lastly, put face velocity into gridfunction object
      update_node_velocity(S, face_dof, face_velocity);
   }
}

/*
Function: update_node_velocity
Parameters:
   S        - BlockVector representing FiniteElement information
   node     - 
   test_vel - Velocity used for testing
Purpose:
   This function computes the node velocities on the faces which
   are designed to bubble in the direction of the normal vector 
   to conserve mass locally.
*/
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
      {
         assert(dim==1);
         return 1.;
      }
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
   bool no_iter = true; // Had to change for test case to run
   int k = 0; // Tells you how many iterations were needed for convergence

   // cout << "CLM pre fortran function.\n";
   __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,
      &no_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&vstar,&k,&b_covolume);

   double d = std::max(std::abs(lambda_maxl_out), std::abs(lambda_maxr_out));

   if (isnan(d))
   {
      cout << "nij:\n";
      n_ij.Print(cout);
      cout << "Ui:\n";
      U_i.Print(cout);
      cout << "Uj:\n";
      U_j.Print(cout);
      cout << "in_rhol: " << in_rhol << ", ul: " << in_ul << ", el: " << in_el << ", pl: " << in_pl << endl;
      cout << "in_rhor: " << in_rhor << ", ur: " << in_ur << ", er: " << in_er << ", pr: " << in_pr << endl;
      MFEM_ABORT("NaN values returned by lambda max computation!\n");
   }

   return d;
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
