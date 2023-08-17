#include "laglos_solver.hpp"
#include <cassert>


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


/****************************************************************************************************
* Function: Constructor
* Parameters:
*  h1            - ParFiniteElementSpace to be used for the mesh motion.  The node corresponding to
*                  the cell center is an averaged valued computed in fill_center_velocities_with_average
*  l2            - ParFiniteElementSpace to be used for the hydrodynamics variables.  We use a DG0
*                  approximation in this code.
*  l2v           - ParFiniteElementSpace for the velocity component.
*  cr            - ParFiniteElementSpace to assist with the mesh velocity reconstruction.  This object
*                  holds the intermediate face velocity object computed in equation (5.7) of the paper.
*  m             - ParLinearForm containing the initial density of each cell to be used to compute the
*                  initial cell mass, which is guaranteed to be conserved locally.
*  _pb           - ProblemBase class containing problem specific functions and quantities.
* Options:
*  use_viscosity - Boolean for the optional addition of viscosity to ensure our method is IDP.
*  mm            - Boolean to allow for mesh motion to guarantee local conservation of mass.
*  CFL           - Double representing our time step restriction.
*
* Purpose: Instantiate LagrangianLOOperator class.
****************************************************************************************************/
template<int dim>
LagrangianLOOperator<dim>::LagrangianLOOperator(ParFiniteElementSpace &h1,
                                                ParFiniteElementSpace &l2,
                                                ParFiniteElementSpace &l2v,
                                                ParFiniteElementSpace &cr,
                                                ParLinearForm *m,
                                                ProblemBase<dim> *_pb,
                                                bool use_viscosity,
                                                bool mm, 
                                                double CFL) :
   H1(h1),
   L2(l2),
   L2V(l2v),
   CR(cr),
   CRc(CR.GetParMesh(), CR.FEColl(), 1),
   v_CR_gf(&CR),
   v_geo_gf(&h1),
   pmesh(H1.GetParMesh()),
   m_lf(m),
   pb(_pb),
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
   cout << "Instantiating hydro op\n";
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

template<int dim>
LagrangianLOOperator<dim>::~LagrangianLOOperator()
{
   // delete pmesh, m_lf, m_hpv;
}


/****************************************************************************************************
* Function: CalculateTimestep
* Parameters:
*  S - BlockVector that stores mesh information, mesh velocity, and state variables.
*
* Purpose:
*  Computes the maximum time step that can be taken under the CFL condition given in section 4.4:
*           dt * Sum_{c' \in J^*(c)} \frac{d_cc'}{mcrho} \leq \frac{CFL}{2} \forall c.
*
*  This is accomplished by iterative over each cell and computing the corresponding max dt, then 
*  a minimum of the current dt, and the computed restriction.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::CalculateTimestep(const Vector &S)
{
   cout << "CalculateTimestep\n";
   // int n = m_hpv->Size(); // NDofs_L2
   double t_min = 1.;
   double t_temp = 0;
   double mi = 0;
   int cj = 0;

   Array<int> fids, oris;
   Vector c(dim), n(dim), n_int(dim);;
   Vector U_i(dim+2), U_j(dim+2);
   double c_norm = 0., d=0., temp_sum = 0.;

   mfem::Mesh::FaceInformation FI;

   for (int ci = 0; ci < NDofs_L2; ci++) // Cell iterator
   {
      cout << "cell: " << ci << endl;
      temp_sum = 0.;
      mi = m_hpv->Elem(ci); // TODO

      GetCellStateVector(S, ci, U_i);

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
         cout << "cell: " << ci << ", face: " << fids[j] << endl;
         CalcOutwardNormalInt(S, ci, fids[j], n_int);
         c = n_int;
         c /= 2.;
         c_norm = c.Norml2();
         n = n_int;
         double n_norm = n.Norml2();
         n /= n_norm;

         
         if (1. - n.Norml2() > 1.e-12)
         {
            cout << "1. - n.Norml2(): " << 1. - n.Norml2() << endl;
            cout << "n.Norml2(): " << n.Norml2() << endl;
            cout << "n vec causing trouble: ";
            n.Print(cout);
            MFEM_ABORT("Normal vec not normal.\n");
         }

         FI = pmesh->GetFaceInformation(fids[j]);

         if (FI.IsInterior())
         {
            cout << "interior face: " << fids[j] << endl;
            // Get index information/state vector for second cell
            if (ci == FI.element[0].index) { 
               cj = FI.element[1].index; 
            }
            else { 
               cj = FI.element[0].index; 
            }

            GetCellStateVector(S, cj, U_j);

            // viscosity contribution
            cout << "getting pressure\n";
            double pl = pb->pressure(U_i);
            double pr = pb->pressure(U_j);
            cout << "pl: " << pl << ", pr: " << pr << endl;
            cout << "Ul: ";
            U_i.Print(cout);
            cout << "Ur: ";
            U_j.Print(cout);
            cout << "computing lambda max\n";
            d = pb->compute_lambda_max(U_i, U_j, n, pl, pr, pb->get_b()) * c_norm; 

            cout << "d for ci " << ci << " and cj " << cj << ": " << d << endl;

            temp_sum += d;
         }
      }
      
      t_temp = 0.5 * ((CFL * mi) / temp_sum );

      if (t_temp < t_min && t_temp > 1e-12) { 
         cout << "t_min: " << t_min << ", t_temp: " << t_temp << endl;
         t_min = t_temp; }
   } // End cell iterator

   this->timestep = t_min;
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
template<int dim>
void LagrangianLOOperator<dim>::GetEntityDof(const int GDof, DofEntity & entity, int & EDof)
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


/****************************************************************************************************
* Function: CreateBdrElementIndexingArray
*  
* Purpose: 
*  To fill the Vector BdrElementIndexingArray of size NumFaces. If a face is a boundary face, then
*  BdrElementIndexingArray[face] = 1, and if the face is an interior face, then we will have
*  BdrElementIndexingArray[face] = -1.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::CreateBdrElementIndexingArray()
{
   // cout << "Constructing BdrElementIndexingArray:\n";
   for (int i = 0; i < pmesh->GetNBE(); i++)
   {
      int index = pmesh->GetBdrElementEdgeIndex(i);
      BdrElementIndexingArray[index] = i;
      // cout << "bdr el index: " << i << ", corresponding edge index: " << index << endl;
   }
}


/****************************************************************************************************
* Function: CreateBdrVertexIndexingArray
*
* Purpose: 
*  To fill the Vector CreateBdrVertexIndexingArray of size Numvertices. If a vertex is a boundary vertex, 
*  then CreateBdrVertexIndexingArray[vertex] = 1, and if the vertex is an interior vertex, then we will 
*  have CreateBdrVertexIndexingArray[vertex] = -1. This is done by iterating over boundary elements, 
*  grab edges, and fill in the corresponding boundary attribute for the adjacent vertices
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::CreateBdrVertexIndexingArray()
{
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
            if (pb->indicator == "saltzman")
            {
               // Replace the bdr attribute in the array as long as it is not
               // the dirichlet condition (For Saltzman Problem)
               if (BdrVertexIndexingArray[index] != 1)
               {
                  BdrVertexIndexingArray[index] = bdr_attr;
                  // cout << "node " << index << " has bdr attr: " << bdr_attr << endl;
               }

            }
            else
            {
               BdrVertexIndexingArray[index] = bdr_attr;
            }

         } // end vertex iterator
      } // end boundary faces
   } // end boundary elements
}


/****************************************************************************************************
* Function: IsBdrVertex
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
bool LagrangianLOOperator<dim>::IsBdrVertex(const int & node)
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
template<int dim>
void LagrangianLOOperator<dim>::MakeTimeStep(Vector &S, const double & t, double & dt)
{
   // cout << "Making timestep\n";
   if (mm)
   {
      // Compute mesh velocities
      // This function may change the timestep
      ComputeMeshVelocities(S, t, dt);
   }

   // Update state variables contained in S_new
   ComputeStateUpdate(S, t, dt);

   // Move the mesh
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction x_gf, mv_gf;
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);
   mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);

   add(x_gf, dt, mv_gf, x_gf);
   pmesh->NodesUpdated();
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
template<int dim>
void LagrangianLOOperator<dim>::ComputeStateUpdate(Vector &S, const double &t, const double dt)
{
   // cout << "========================================\n";
   // cout << "Computing State Update\n";
   // cout << "========================================\n";

   // We need a place to store the new state variables
   Vector S_new = S;

   Vector val(dim+2), c(dim), n(dim), n_int(dim), U_i(dim+2), U_j(dim+2), sums(dim+2);
   Array<int> fids, fids2, oris, oris2, verts; // TODO: Remove fids2, oris2, These are for testing.

   int cj = 0;
   double d, c_norm;

   mfem::Mesh::FaceInformation FI;
   H1.ExchangeFaceNbrData();

   bool is_boundary_cell = false;

   Vector sum_validation(dim);
   for (int ci = 0; ci < NDofs_L2; ci++) // Cell iterator
   {
      is_boundary_cell = false;

      sum_validation = 0.;
      // cout << "CSU::ci: " << ci << endl;
      GetCellStateVector(S, ci, U_i);
      // cout << "cell state vector:\n";
      // U_i.Print(cout);
      
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

      DenseMatrix F_i = pb->flux(U_i);
      sums = 0.;
      // cout << "Flux at U_i:\n";
      // F_i.Print(cout);

      for (int j=0; j < fids.Size(); j++) // Face iterator
      {
         CalcOutwardNormalInt(S, ci, fids[j], n_int);
         c = n_int;
         c /= 2.;
         sum_validation.Add(1., c);
         c_norm = c.Norml2();
         n = n_int;
         double F = n.Norml2();
         n /= F;
         assert(1. - n.Norml2() < 1e-12);
         FI = pmesh->GetFaceInformation(fids[j]);
         // cout << "\tcell: " << ci << ", face: " << fids[j] << endl;

         if (FI.IsInterior())
         {
            // cout << "\t\tCSU::Interior face\n";
            // Get index information/state vector for second cell
            if (ci == FI.element[0].index) { 
               cj = FI.element[1].index; 
            }
            else { 
               assert(ci == FI.element[1].index);
               cj = FI.element[0].index; 
            }
            GetCellStateVector(S, cj, U_j); 

            // flux contribution
            DenseMatrix dm = pb->flux(U_j);
            dm += F_i; 
            Vector y(dim+2);
            dm.Mult(c, y);

            sums -= y;

            /* viscosity contribution */
            if (use_viscosity)
            {
               double pl = pb->pressure(U_i);
               double pr = pb->pressure(U_j);
               d = pb->compute_lambda_max(U_i, U_j, n, pl, pr, pb->get_b()) * c_norm; 
               // cout << "viscosity: " << d << endl;
               Vector z = U_j;
               z -= U_i;
               sums.Add(d, z);
            }
            
         }
         else
         {
            // cout << "\t\tCSU::boundary\n";
            assert(FI.IsBoundary());
            is_boundary_cell = true;
            Vector y_temp(dim+2), y_temp_bdry(dim+2), U_i_bdry(dim+2);

            F_i.Mult(c, y_temp);
            U_i_bdry = U_i;
            y_temp_bdry = 0.;

            /* Enforce Boundary Conditions */
            // if (pb->indicator == "noh")
            // {
            //    // Since we will enforce Dirichlet conditions on all boundary cells, 
            //    // this enforcement will be done after the face iterator
            // }
            // else if (pb->indicator == "saltzman")
            // {
            //    // Check bdry flag
            //    int BdrElIndex = BdrElementIndexingArray[fids[j]];
            //    int bdr_attribute = pmesh->GetBdrAttribute(BdrElIndex);

            //    if (bdr_attribute == 1) // Dirichlet (star states "ghost")
            //    {
            //       double _rho = 3.9992502342988532;
            //       double _p = 1.3334833281256551;
            //       Vector _v(dim);
            //       _v = 0.;
            //       _v[0] = 1.; //e_x
            //       Vector _v_neg = _v, _vp = _v;
            //       _v_neg *= -1.;
            //       _vp *= _p;

            //       DenseMatrix F_i_bdry(dim+2, dim);
            //       F_i_bdry.SetRow(0, _v_neg);
            //       for (int i = 0; i < dim; i++)
            //       {
            //          F_i_bdry(i+1, i) = _p;
            //       }
            //       F_i_bdry.SetRow(dim+1, _vp);

            //       F_i_bdry.Mult(c, y_temp_bdry);
            //    }

            //    else if (bdr_attribute == 2) // slip
            //    {
            //       // Negate velocity
            //       for (int _it = 0; _it < dim; _it++)
            //       {
            //          U_i_bdry[_it + 1] = U_i_bdry[_it + 1] * -1;
            //       }
            //       DenseMatrix F_i_slip = pb->flux(U_i_bdry);
            //       F_i_slip.Mult(c, y_temp_bdry);
            //    }
            //    else
            //    {
            //       cout << "invalid boundary attribute: " << bdr_attribute << endl;
            //       cout << "cell : " << ci << endl;
            //       cout << "face: " << fids[j] << endl;  
            //       y_temp *= 2.; 
            //    }

            //    y_temp += y_temp_bdry;
            // }
            
            y_temp *= 2.;
            

            // Add in boundary contribution
            sums -= y_temp;
         } // End boundary face        

      } // End Face iterator

      // Enforce exact condition on boundary for Noh
      // if (pb->indicator == "noh" && is_boundary_cell)
      // {
      //    // Compute cell center and corresponding velocity at this location
      //    // Get center node dof
      //    int cell_vdof;
      //    Vector cell_x(dim);
      
      //    switch (dim)
      //    {
      //       case 1:
      //       {
      //          cell_vdof = num_faces + ci;
      //          break;
      //       }
      //       case 2:
      //       {
      //          cell_vdof = NVDofs_H1 + num_faces + ci;
      //          break;
      //       }
      //       case 3:
      //       {
      //          MFEM_ABORT("3D not implemented.\n");
      //       }
      //       default:
      //       {
      //          MFEM_ABORT("Incorrect dim value provided.\n");
      //       }
      //    }

      //    get_node_position(S, cell_vdof, cell_x);

      //    // TODO: Fill in val with exact solution
      //    val[0] = pb->rho0(cell_x, t+dt); // Is this t or t + dt?
      //    Vector v_exact(dim);
      //    pb->v0(cell_x, t+dt, v_exact);
      //    for (int i = 0; i < dim; i++)
      //    {
      //       val[1 + i] = v_exact[i];
      //    }
      //    val[dim + 1] = pb->ste0(cell_x, t+dt);
      // }
      // else
      // {
      // Do the normal thing
      assert(sum_validation.Norml2() < 1e-12);

      sums *= dt;
      double k = pmesh->GetElementVolume(ci);
      double _mass = k / U_i[0];
      sums /= _mass;
      val += sums;
      // }

      // In either case, update the cell state vector
      SetCellStateVector(S_new, ci, val);
      
   } // End cell iterator
   S = S_new;
}


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::GetCellStateVector(const Vector &S, 
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
      int index = cell + i*NDofs_L2;
      assert(index < v_gf.Size());
      U[i+1] = v_gf.Elem(index);
   }

   // Retrieve cell specific total energy
   U[dim+1] = ste_gf.Elem(cell);
}


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::SetCellStateVector(Vector &S_new,
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
      dofs.Append(cell + i*NDofs_L2);
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
template<int dim>
void LagrangianLOOperator<dim>::CalcOutwardNormalInt(const Vector &S, const int cell, const int face, Vector & res)
{
   cout << "=======================================\n"
        << "  Calculating outward normal integral  \n"
        << "=======================================\n";
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
         // cout << "face dofs row for face: " << face << endl;
         // row.Print(cout);
         int face_vdof1 = row[1], face_vdof2 = row[0];
         face_dof = row[2];

         Vector face_x(dim), vdof1_x(dim), vdof2_x(dim), vdof1_v(dim), vdof2_v(dim);
         get_node_position(S, face_dof, face_x);
         get_node_position(S, face_vdof1, vdof1_x);
         get_node_position(S, face_vdof2, vdof2_x);
         // cout << "caloutwardnormalint, cell: " << cell << ", face: " << face << endl;
         // cout << "face_dof: " << row[2] << endl;
         // cout << "face_x: ";
         // face_x.Print(cout);
         // cout << "vdof1_x: ";
         // vdof1_x.Print(cout);
         // cout << "vdof2_x: ";
         // vdof2_x.Print(cout);

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

         cout << "cell: " << cell << ", face: " << face << ", outwardnormint: ";
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


/****************************************************************************************************
* Function: Orthogonal
* Parameters:
*     v - Vector to be rotated
*
* Purpose:
*     Rotate a vector counterclockwise of angle pi/2.
*
* Example:
*  (0,-1) ---> (1,0)
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::Orthogonal(Vector &v)
{
   if (dim == 2)
   {
      double x = v(0), y = v(1);
      v(0) = -1. * y;
      v(1) = x;
   }
}

/*
*
* Enforce Boundary Conditions
*  Enforces on v_gf
*
*/
// template<int dim>
// void LagrangianLOOperator<dim>::EnforceBoundaryConditions(Vector &S)
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


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::
   compute_intermediate_face_velocities(const Vector &S,
                                        const double t,
                                        const string flag, // Default NA
                                        void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
   /***
    * Compute intermediate face velocities according to (5.7) 
    ***/
   cout << "Compute intermediate face velocities\n";
   mfem::Mesh::FaceInformation FI;
   int c, cp;
   Vector Uc(dim+2), Ucp(dim+2), n_int(dim), c_vec(dim), n_vec(dim), Vf(dim); 
   double d, c_norm, F;

   for (int face = 0; face < num_faces; face++) // face iterator
   {
      cout << "ifv for face: " << face << endl;
      Vf = 0.;
      FI = pmesh->GetFaceInformation(face);
      c = FI.element[0].index;
      cp = FI.element[1].index;

      GetCellStateVector(S, c, Uc);

      if (flag == "NA")
      {
         if (FI.IsInterior())
         {
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

            // Compute max wave speed
            double pl = pb->pressure(Uc);
            double pr = pb->pressure(Ucp);
            double val = pb->compute_lambda_max(Uc, Ucp, n_vec, pl, pr, pb->get_b());
            d = val * c_norm; 

            // Compute intermediate face velocity
            Vf = pb->velocity(Uc);
            Vf += pb->velocity(Ucp);
            Vf *= 0.5;

            double coeff = d * (Ucp[0] - Uc[0]) / F;
            Vf.Add(coeff, n_vec);

            double coeff2 = val / 2.;
            coeff2 *= (Ucp[0] - Uc[0]);
            cout << "coeff: " << coeff << ", coeff2: " << coeff2 << endl;
            
            // Debugging purposes
            if (coeff != 0.)
            {
               cout << "ifv::coeff: " << coeff << endl;
               cout << "cell 1: " << c << ", pressure: " << pb->pressure(Uc) << ", density: " <<  1./Uc[0] << ", sv: " << Uc[0] << endl;
               cout << "cell 2: " << cp << ", pressure: " << pb->pressure(Ucp) << ", density: " <<  1./Ucp[0] << ", sv: " << Ucp[0] << endl;
               cout << "mws: " << val << ", d: " << d << endl;

               cout << "F: " << F << endl;
               cout << "n_vec:\n";
               n_vec.Print(cout);
               cout << "final face velocity:\n";
               Vf.Print(cout);
            } 
         }
         else 
         {
            assert(FI.IsBoundary());
            Vf = pb->velocity(Uc);             
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


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::get_intermediate_face_velocity(const int & face, Vector & vel)
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
* Function: ComputeMeshVelocities
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
template<int dim>
void LagrangianLOOperator<dim>::ComputeMeshVelocities(Vector &S, 
                                                      const double & t, 
                                                      double & dt,
                                                      const string flag, // Default NA
                                                      void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL)
{
   cout << "Compute mesh velocities\n";
   // Vector* sptr = const_cast<Vector*>(&S);
   // ParGridFunction x_gf, mv_gf;
   // x_gf.MakeRef(&H1, *sptr, block_offsets[0]);
   // mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);

   compute_intermediate_face_velocities(S, t, flag, test_vel);  

   // Compute v_geo_gf
   compute_geo_V();

   compute_node_velocities(S, t, dt);
   if (dim > 1)
   {
      // Don't need to fill face velocities with average in dim=1
      // fill_face_velocities_with_average(S);
      compute_corrective_face_velocities(S, t, dt, flag, test_vel);
      // S.Print(cout);

      // Construct ti, face normals
      // ti = [0. 0.]^T in 2D
      // 3DTODO: Modify this according to seciton 5.2
   }
   else
   {
      // Validate conditions for mass conservation in 1D given in Theorem 5.4
      Array<int> fids;
      Vector n_int(dim), n(dim), Vf(dim);
      double val = 0., k=0., temp_sum = 0.;
      for (int ci = 0; ci < NDofs_L2; ci++)
      {  
         // cout << "cell: " << ci << endl;
         val = 1.;
         k = pmesh->GetElementVolume(ci);
         temp_sum = 0.;
         pmesh->GetElementVertices(ci, fids);
         for (int j=0; j < fids.Size(); j++) // Face iterator
         {
            CalcOutwardNormalInt(S, ci, fids[j], n);
            get_intermediate_face_velocity(fids[j], Vf);
            temp_sum += n[0] * Vf[0];
         }
         val += (dt / k) * temp_sum;
         assert(val > 0.);
      }
   }
   fill_center_velocities_with_average(S);
}


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
double LagrangianLOOperator<dim>::CalcMassLoss(const Vector &S)
{
   // cout << "Checking Mass Loss\n";
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


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::CheckMassConservation(const Vector &S)
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
   double cell_ratio = (double)counter / (double)NDofs_L2;

   cout << "Percentage of cells where mass conservation was broken: " << cell_ratio << endl;
}


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::tensor(const Vector & v1, const Vector & v2, DenseMatrix & dm)
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


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
// Assumes dt, not dt/2 is passed in
template<int dim>
void LagrangianLOOperator<dim>::compute_determinant(const DenseMatrix &C, const double &dt, double & d)
{
   cout << "=====================\n";
   cout << "Computing determinant\n";
   cout << "=====================\n";
   double trace = C.Trace();
   double det = C.Det();
   cout << "Matrix C:\n";
   C.Print(cout);
   cout << "trace: " << trace << ", det: " << det << endl;

   double a = 1.;
   double b = -1. * (1. + (dt / 2.) * trace);
   double c = (pow(dt, 2) / 4.) * det;

   cout << "tau: " << dt / 2. << ", a: " << a << ", b: " << b << ", c: " << c << endl;

   cout << "discriminant: " << pow(b,2) - 4. * a * c << endl;
   double pm = sqrt(pow(b,2) - 4. * a * c);
   cout << "pm: " << pm << endl;

   double d1 = -1. * b + pm;
   d1 /= (2. * a);

   double d2 = -1. * b - pm;
   d2 /= (2. * a);

   d = std::max(d1, d2);

   if (d <= 0.)
   {
      cout << "d1: " << d1 << ", d2: " << d2 << endl;
      MFEM_ABORT("Alpha_i should be positive.\n");
   }

   cout << "d1: " << d1 << endl;
   cout << "d2: " << d2 << endl;
   cout << "d: " << d << endl;
}


/****************************************************************************************************
* Function: compute_node_velocity_RT
* Parameters:
*  dt - timestep
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::
   compute_node_velocity_RT(const int & node, double & dt, Vector &node_v, bool &is_dt_changed)
{
   cout << "-----\nCompute node velocity RT\n-----\n";
   // cout << "Node: " << node << endl;
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

         // Enforce time restrition imposed by calculation of alpha_i
         if (dim == 2)
         {
            // Verify determinant > 0
            double trace = Ci.Trace();
            double det = Ci.Det();
            if (det <= 1.e-12)
            {
               cout << "Negative determinant.\n";
               double _c = abs(Ci(0,0));
               _c += 1.;
               if (dt > 2. / _c)
               {
                  dt = 2. / _c;
                  is_dt_changed = true;
                  cout << "averting alpha_i <= 0 abort call.\n";
               }
            }
            
            // Enforce timestep restriction
            double val = 2. * sqrt(det);
            double zero1 = 1. / (-trace - val);
            double zero2 = 1. / (-trace + val);

            cout << "trace: " << trace << ", det: " << det << ", z1: " << zero1 << ", z2: " << zero2 << endl;

            /* old timestep restriction */
            if (zero1 <= 0.)
            {
               if (zero2 > 0.) // z1 0 z2
               {
                  // We must enforce dt <= 2*z2
                  if (dt > 2.*zero2) { 
                     cout << "restricting timestep, z1,0,z2\n";
                     cout << "old timestep: " << dt << ", new timestep: " << zero2 << endl;
                     dt = 2.*zero2; 
                     is_dt_changed = true;
                  }
               }
               // else // z1 z2 0 - no enforcement on timestep
            }
            else // 0 z1 z2
            {
               if (dt > 2.*zero1) // 0 z1 dt z2    --- or --- 0 z1 z2 dt
               {
                  cout << "restricting timestep, 0,z1,z2\n";
                  cout << "old timestep: " << dt << ", new timestep: " << zero1 << endl;
                  dt = 2.*zero1;
                  is_dt_changed = true;
               }
            }

            /* new timestep restriction */
            // No time restriction if b < 0
            // a = trace, b = determinant
            // if (det >= 0.)
            // {
            //    if (pow(trace,2) - 4 * det > 0.)
            //    {
            //       // if (trace > 0) -> no timestep restriction
            //       if (trace <= 0.)
            //       {
            //          if (dt > 2. * zero1)
            //          {
            //             cout << "timestep restricted (1)\n";
            //             is_dt_changed = true;
            //             dt = 2. * zero1;
            //          }
            //       }
            //    }
            //    else if (pow(trace,2) - 4 * det < 0.)
            //    {
            //       if (dt > 2. * zero2)
            //       {
            //          cout << "timestep restricted (2)\n";
            //          is_dt_changed = true;
            //          dt = 2. * zero2;
            //       }
            //    }
            //    else
            //    {
            //       assert (pow(trace,2) - 4 * det == 0.);
            //       if (trace < 0)
            //       {
            //          if (dt > -1. / trace)
            //          {
            //             cout << "timestep restricted (2)\n";
            //             is_dt_changed = true;
            //             dt = -1. / trace;
            //          }
            //       }
            //    }
            // }
         }
         get_vi_geo(node, Vgeo);
         cout << "cnv_RT Vgeo: ";
         Vgeo.Print(cout);

         compute_determinant(Ci, dt, d);

         // Compute V_i^n
         DenseMatrix _mat(Ci);
         _mat *= - dt / 2.;
         for (int i = 0; i < dim; i++)
         {
            _mat(i,i) += d;
         }

         cout << "det: " << d << endl;
         cout << "Ci: ";
         Ci.Print(cout);
         cout << "Pre inverse: \n";
         _mat.Print(cout);

         _mat.Invert();
         cout << "Inverse: \n";
         _mat.Print(cout);

         _mat.Mult(Vgeo, node_v);
         cout << "cnv_RT post mult V geo: ";
         Vgeo.Print(cout);
         cout << "node_v for node " << node << ": ";
         node_v.Print(cout);
         break;
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
template<int dim>
void LagrangianLOOperator<dim>::RT_int_grad(const int cell, DenseMatrix & res)
{
   cout << "RT_int_grad funcall on cell: " << cell << endl;
   ParGridFunction CRc_gf(&CRc);
   const int _size = CRc.GetVSize();

   ElementTransformation * trans = pmesh->GetElementTransformation(cell);
   Vector grad(dim), row(dim);
   res.SetSize(dim);
   res = 0., row = 0.;

   // Iterate over quadrature
   cout << "num Q points: " << RT_ir.GetNPoints() << endl;
   for (int i = 0; i < RT_ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = RT_ir.IntPoint(i);
      trans->SetIntPoint(&ip);
      cout << "ip.x: " << ip.x << ", ip.y: " << ip.y << ", weight: " << ip.weight << endl;

      cout << "el " << cell << " at integration point " << i << endl;
      for (int j = 0; j < dim; j++)
      {
         cout << "dim: " << j << endl;
         CRc_gf.MakeRef(&CRc, v_CR_gf, j*_size);
         CRc_gf.GetGradient(*trans, grad);
         // if (grad.Norml2() > 0.)
         // {
         //    cout << "grad: ";
         //    grad.Print(cout);
         //    cout << "at cell: " << cell << endl;
         // }

         // cout << "v_CR_gf: " << endl;
         // v_CR_gf.Print(cout);
         // cout << "CRc_gf: " << endl;
         // CRc_gf.Print(cout);

         cout << "grad: ";
         grad.Print(cout);

         // Put information into Dense Matrix
         res.GetRow(j, row);
         row.Add(ip.weight, grad);
         res.SetRow(j, row);
      }
   }
   // if (res.FNorm() > 0.)
   // {
   //    cout << "Resulting nonzero matrix:\n";
   //    res.Print(cout);
   //    assert(false);
   // }  
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
template<int dim>
void LagrangianLOOperator<dim>::compute_geo_C(const int &node, DenseMatrix & res)
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
         cout << "compute_geo_C::corner node\n";
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
         MFEM_ABORT("No need to compute C_i at cell centers.\n");
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
      double cell_vol = pmesh->GetElementVolume(row_el);
      denom += cell_vol;

      dm_temp = 0.;
      RT_int_grad(row_el, dm_temp);
      // res.Add(1., dm_temp);
      res.Add(cell_vol, dm_temp);
      cout << "cell_vol: " << cell_vol << endl;
      cout << "int_grad for row el: " << row_el << endl;
      dm_temp.Print(cout);
   }

   cout << "denom: " << denom << endl;
   cout << "res before dividing: \n";
   res.Print(cout);
   res *= 1./denom;
   cout << "Ci for node " << node << ": \n";
   res.Print(cout);
}

/***********************************************************************************************************
Function: compute_geo_V
Purpose:
   Reconstruct and average RT velocity field.  This projects the discontinuous GF from the 
   Crouzeix-Raviart FE space onto the H1 space, filling the valies into v_geo_gf. (5.11)

   NOTE: This function assumes that the function LagrangianLOOperator:compute_intermediate_face_velocities()
   has already been called.  If this function has not been called, then the returned velocity will be 0.
***********************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::compute_geo_V()
{
   VectorGridFunctionCoefficient vcr_coeff(&v_CR_gf);
   v_geo_gf.ProjectDiscCoefficient(vcr_coeff, mfem::ParGridFunction::AvgType::ARITHMETIC);  
}


template<int dim>
void LagrangianLOOperator<dim>::
   compute_node_velocities(Vector &S, 
                           const double & t, 
                           double & dt, 
                           const string flag, // Default NA
                           void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
   cout << "Compute node velocities\n";
   Vector vertex_v(dim);

   bool is_dt_changed = false;
   // Iterate over vertices
   for (int vertex = 0; vertex < num_vertices; vertex++) // Vertex iterator
   {
      compute_node_velocity_RT(vertex, dt, vertex_v, is_dt_changed);

      if (vertex_v[0] != vertex_v[0] || vertex_v[1] != vertex_v[1])
      {
         cout << "NaN velocity encountered in compute_node_velocities at vertex: " << vertex << endl;
         MFEM_ABORT("Aborting due to NaNs.\n");
      }

      update_node_velocity(S, vertex, vertex_v);

      // If we restricted the timestep, we must recompute the vertex velocities that were computed previously
      if (is_dt_changed)
      {
         vertex = -1;
         is_dt_changed = false;
         cout << "Restarting vertex iterator\n";
      }
   } // End Vertex iterator
   // assert(false);
}


/***********************************************************************************************************
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
***********************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::
   compute_corrective_face_velocities(Vector &S, const double & t, const double & dt,
                                    const string flag, // Default NA
                                    void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
   assert(dim > 1); // No need to correct face velocities in dim=1
   cout << "========================================\n";
   cout << "Computing corrective interior face velocities.\n";
   cout << "========================================\n";
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
         cell_c_v = pb->velocity(Uc);
         cell_cp_v = pb->velocity(Ucp);
         
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

         if (flag == "testing")
         {
            cout << "The coordinate location of face " << face << " is: ";
            face_x.Print(cout);
            cout << "The computed velocity on this face is: ";
            face_velocity.Print(cout);
            cout << "The exact face velocity is: ";
            Vector face_v_exact(dim);
            test_vel(face_x, 0., face_v_exact);
            face_v_exact.Print(cout);
            Vector temp_vel(dim);
            // subtract(face_v_exact, face_velocity, temp_vel);

            // if (temp_vel.Norml2() > 10e-6)
            // {
            //    MFEM_ABORT("Incorrect face velocity.\n");
            // }
         }
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
      if (face_velocity[0] != face_velocity[0] || face_velocity[1] != face_velocity[1])
      {
         cout << "Houston we have a problem.\n";
         MFEM_ABORT("Exit");
      }
      update_node_velocity(S, face_dof, face_velocity);
   }

   // assert(false);
}


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::
   fill_center_velocities_with_average(Vector &S,
                                       const string flag, 
                                       void (*test_vel)(const Vector&, const double&, Vector&))
{
   // Since we cannot use Serendipity elements, we must update cell center velocities
   Vector Vc(dim), Uc(dim+2), node_v(dim);

   Array<int> verts;
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
      
      // GetCellStateVector(S, ci, Uc);
      // Vc = pb->velocity(Uc);
      Vc = 0.;
      if (dim == 1)
      {
         pmesh->GetElementVertices(ci, verts);
         for (int j = 0; j < verts.Size(); j++)
         {
            get_node_velocity(S, verts[j], node_v);
            Vc += node_v;
         }
         Vc /= verts.Size();
      }
      else if (dim == 2)
      {
         Array<int> element_face_row, element_face_oris;
         pmesh->GetElementEdges(ci, element_face_row, element_face_oris);
         for (int face_it = 0; face_it < element_face_row.Size(); face_it++)
         {
            int face = element_face_row[face_it];

            // Get intermediate face velocities corresponding to faces
            // get_intermediate_face_velocity(face, node_v);
            get_node_velocity(S, face + NVDofs_H1, node_v);
            if (face_it%2 == 0)
            {
               Vc[1] = Vc[1] + node_v[1];
            }
            else
            {
               Vc[0] = Vc[0] + node_v[0];
            }
         }

         Vc /= 2.;
      }

      // Update the velocity for the center node
      if (Vc[0] != Vc[0] || Vc[1] != Vc[1])
      {
         MFEM_ABORT("Center velocity average resulted in NaN.\n");
      }
      update_node_velocity(S, cell_vdof, Vc);
   }
}


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::fill_face_velocities_with_average(Vector &S, const string flag, void (*test_vel)(const Vector&, const double&, Vector&))
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
template<int dim>
void LagrangianLOOperator<dim>::update_node_velocity(Vector &S, const int & node, const Vector & vel)
{

   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction mv_gf;
   mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);

   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      // if (vel[i] != vel[i]) { 
      //    cout << "Node: " << node << endl;
      //    MFEM_ABORT("NaN val encountered."); 
      // }
      mv_gf[index] = vel[i];
   }
}


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::get_node_velocity(const Vector &S, const int & node, Vector & vel)
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


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::get_vi_geo(const int & node, Vector & vel)
{
   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      vel[i] = v_geo_gf[index];
   }
}


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::get_node_position(const Vector &S, const int & node, Vector & x)
{
   // cout << "=======================================\n"
   //      << "          Get Node Positions           \n"
   //      << "=======================================\n";
   // cout << "node: " << node << ", NDofs_H1: " << NDofs_H1 << endl;
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction x_gf;
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);

   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      // cout << "i: " << i << ", val: " << x_gf[index] << endl;
      x[i] = x_gf[index];
   }
}


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::SaveStateVecsToFile(const Vector &S, const string &output_file_prefix, const string &output_file_suffix)
{
   // Get gfs
   ParGridFunction rho_gf, sv_gf, v_gf, ste_gf;
   Vector* sptr = const_cast<Vector*>(&S);
   sv_gf.MakeRef(&L2, *sptr, block_offsets[2]);
   v_gf.MakeRef(&L2V, *sptr, block_offsets[3]);
   ste_gf.MakeRef(&L2, *sptr, block_offsets[4]);

   Vector center(dim);

   // Fill density

   // Form filenames and ofstream objects
   std::string sv_file = output_file_prefix + "sv_" + output_file_suffix;
   std::ofstream fstream_sv(sv_file.c_str());
   fstream_sv << "x,rho,v\n";

   
   for (int i = 0; i < NDofs_L2; i++)
   {
      pmesh->GetElementCenter(i, center);

      fstream_sv << center[0] << ","
                 << 1./sv_gf[i] << ","
                 << v_gf[i] << "\n";
   }

   // fstream_v.close();

   // const int nqp = ir.GetNPoints();
   // Vector pos(dim);
   // for (int e = 0; e < NE; e++)
   // {
   //    ElementTransformation &Tr = *L2.GetElementTransformation(e);
   //    for (int q = 0; q < nqp; q++)
   //    {
   //       const IntegrationPoint &ip = ir.IntPoint(q);
   //       Tr.SetIntPoint(&ip);
   //       Tr.Transform(ip, pos);

   //       double r = sqrt(pos(0)*pos(0) + pos(1)*pos(1));

   //       double rho = qdata.rho0DetJ0w(e*nqp + q) / Tr.Weight() / ip.weight;
   //       fstream_rho << r << " " << rho << "\n";
   //       fstream_rho.flush();
   //    }
   // }
   // fstream_rho.close();
}

/* Explicit n_vectantiation */
template class LagrangianLOOperator<1>;
template class LagrangianLOOperator<2>;
template class LagrangianLOOperator<3>;

} // end ns hydrodynamics


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


/****************************************************************************************************
* Function: 
* Parameters:
*
* Purpose:
****************************************************************************************************/
void random_mesh_v(const Vector &x, Vector &res)
{
   const double r_max = 1.;
   res[0] = fRand(0, r_max);
   res[1] = fRand(0, r_max);
   return;
}

} // end ns mfem
