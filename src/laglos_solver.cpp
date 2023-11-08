#include "laglos_solver.hpp"
#include <cassert>


namespace mfem
{

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
*                  the cell center is an averaged valued computed in FillCenterVelocitiesWithAvg
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
                                                ParFiniteElementSpace &h1_l,
                                                ParFiniteElementSpace &l2,
                                                ParFiniteElementSpace &l2v,
                                                ParFiniteElementSpace &cr,
                                                ParLinearForm *m,
                                                ProblemBase<dim> *_pb,
                                                bool use_viscosity,
                                                bool mm, 
                                                double CFL) :
   H1(h1),
   // H1_L(h1_l),
   // H1c(H1_L.GetParMesh(), H1_L.FEColl(), 1),
   H1c(H1.GetParMesh(), H1.FEColl(), 1),
   L2(l2),
   L2V(l2v),
   CR(cr),
   CRc(CR.GetParMesh(), CR.FEColl(), 1),
   v_CR_gf(&CR),
   v_CR_gf_corrected(&CR), 
   v_CR_gf_fluxes(&CR),
   v_geo_gf(&H1),
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
   num_elements(L2.GetNE()),
   num_faces(L2.GetNF()),
   num_vertices(pmesh->GetNV()),
   num_edges(pmesh->GetNEdges()),
   vertex_element(pmesh->GetVertexToElementTable()),
   face_element(pmesh->GetFaceToElementTable()),
   edge_vertex(pmesh->GetEdgeVertexTable()),
   // Options
   use_viscosity(use_viscosity),
   mm(mm),
   CFL(CFL)
{
   // Transpose face_element to get element_face
   Transpose(*face_element, element_face);
   Transpose(*edge_vertex, vertex_edge);

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
      lambda_max_vec.SetSize(num_faces);
   }
   
   // Initialize values of intermediate face velocities
   v_geo_gf = 0.;
   v_CR_gf = 0.;
   v_CR_gf_corrected = 0.;
   v_CR_gf_fluxes = 0.;
   // assert(v_CR_gf.Size() == dim * num_faces);
   // cout << "v_CR_gf.Size after: " << v_CR_gf.Size() << endl;

   // Initialize Dij sparse
   InitializeDijMatrix();   

   // Set integration rule for Rannacher-Turek space
   IntegrationRules IntRulesLo(0, Quadrature1D::GaussLobatto);
   RT_ir = IntRules.Get(CR.GetFE(0)->GetGeomType(), RT_ir_order);

   // Reset chronos
   chrono_mm.Clear();
   chrono_state.Clear();

   // Print some dimension information
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
   cout << "each element in the mesh has " << el_num_faces << " faces." << endl;
   cout << "num_elements: " << num_elements << endl;
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
* Function: InitializeDijMatrix
*
* Purpose:
*  The purpose of this function is to initialize a SparseMatrix with the sparsity pattern
*  that corresponds to the adjacency matrix of the mesh.  I.E. this functions constructs a matrix
*  of size num_elements x num_elements where the nonzero elements are for each element row
*  are the cell indices in the stencil. This SparseMatrix will encode the max wave speed of
*  each corresponding local Riemann problem.  An example of cells in the stencil of cell c is given
*  below.
*
*                  -----------
*                 | 0 | 1 | 0 |
*                 |-----------|
*                 | 1 | c | 1 |
*                 |-----------|
*                 | 0 | 1 | 0 |
*                  -----------
* 
* This is accomplished in the following steps:
*  1) Create ParBilinearForm based on L2FESpace
*  2) Create DGTraceIntegrator(DG state coeff, CG coeff (ones), double, double)
*  3) AddInteriorFaceIntegrator using ^
*  4) Assemble, Finalize, ParallelAssemble calls
*  5) HypreParMatrix::MergeDiagAndOffd --> SparseMatrix
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::InitializeDijMatrix()
{
   switch (dim)
   {
      case 1: 
      {
         el_num_faces = pmesh->GetElement(0)->GetNVertices();
         break;
      }
      case 2: 
      {
         el_num_faces = pmesh->GetElement(0)->GetNEdges();
         break;
      }
      case 3: 
      {
         el_num_faces = pmesh->GetElement(0)->GetNFaces();
         break;
      }
      default: MFEM_ABORT("Wrong dimension given.\n");
   }

   // Initialize SparseMatrix
   dij_sparse = new SparseMatrix(num_elements, num_elements, el_num_faces);

   // Create dummy coefficients
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> rho0_static = 
      std::bind(&ProblemBase<dim>::rho0, pb, std::placeholders::_1, std::placeholders::_2);

   FunctionCoefficient rho_coeff(rho0_static);

   Vector ones(dim); ones = 1.;
   VectorConstantCoefficient ones_coeff(ones);

   // Assemble SparseMatrix object
   ParBilinearForm k(&L2);
   k.AddInteriorFaceIntegrator(new DGTraceIntegrator(rho_coeff, ones_coeff, 1.0, 1.0));

   k.Assemble();
   k.Finalize();

   HypreParMatrix * k_hpm = k.ParallelAssemble();
   k_hpm->MergeDiagAndOffd(*dij_sparse);

   // From here, we can modify the sparse matrix according to the sparsity pattern
}


/****************************************************************************************************
* Function: BuildDijMatrix
* Parameters:
*  S - BlockVector that stores mesh information, mesh velocity, and state variables.
*
* Purpose:
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::BuildDijMatrix(const Vector &S)
{
   // cout << "=======================================\n"
   //      << "           Build Dij Matrix            \n"
   //      << "=======================================\n";
   mfem::Mesh::FaceInformation FI;
   int c, cp;
   Vector Uc(dim+2), Ucp(dim+2), n_int(dim), c_vec(dim), n_vec(dim);
   double F, lambda_max, d;

   for (int face = 0; face < num_faces; face++) // face iterator
   {
      // cout << "face: " << face << endl;
      FI = pmesh->GetFaceInformation(face);
      c = FI.element[0].index;
      cp = FI.element[1].index;

      GetCellStateVector(S, c, Uc);
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

         /* Compute max wave speed */ 
         // Some cases require an update to b_covolume at every interface.  This can be done through
         // the function ProblemBase::lm_update(), which is overridden only in the functions it is used.
         double b_covolume = .1 / (max(1./Uc[0], 1./Ucp[0]));
         pb->lm_update(b_covolume);

         // Compute pressure with given EOS
         double pl = pb->pressure(Uc);
         double pr = pb->pressure(Ucp);

         // Finally compute lambda max
         // cout << "pre compute lambda max\n";
         lambda_max = pb->compute_lambda_max(Uc, Ucp, n_vec, pl, pr, b_covolume);
         d = lambda_max * c_norm; 

         double ss = pb->sound_speed(Uc);

         // cout << "c: " << c 
         //      << ", pl: " << pl 
         //      << ", pr: " << pr
         //      << ", cp: " << cp << endl;
         // cout << std::setprecision(12) 
         //      << "sound speed: " << ss 
         //      << ", lambda_max: " << lambda_max 

         //      << ", c_norm: " << c_norm 
         //      << ", d: " << d << endl;

         dij_sparse->Elem(c,cp) = d;
         dij_sparse->Elem(cp,c) = d;

         if (dim == 1)
         {
            lambda_max_vec[face] = lambda_max; // TODO: remove, only temporary
         }
      }
   } // End face iterator

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
   // cout << "CalculateTimestep\n";
   double t_min = 1.;
   double t_temp = 0;
   double mi = 0;
   int cj = 0;

   Array<int> fids, oris;
   Vector U_i(dim+2), U_j(dim+2);
   double d=0., temp_sum = 0.;

   mfem::Mesh::FaceInformation FI;

   for (int ci = 0; ci < NDofs_L2; ci++) // Cell iterator
   {
      // cout << "\tcell: " << ci << endl;
      temp_sum = 0.;
      mi = m_hpv->Elem(ci); 

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

            d = dij_sparse->Elem(ci, cj); 

            // cout << "face: " << fids[j] << endl;
            // cout << "d for ci " << ci << " and cj " << cj << ": " << d << endl;

            temp_sum += d;
         }
      }
      
      t_temp = 0.5 * ((CFL * mi) / temp_sum );

      // cout << "CFL: " << CFL << ", mi: " << mi << ", temp_sum: " << temp_sum << endl;
      // cout << "t_temp: " << t_temp << ", t_min: " << t_min << endl;

      if (t_temp < t_min && t_temp > 1.e-12) { 
         // cout << "timestep reduced\n";
         t_min = t_temp;
      }
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
         pmesh->GetEdgeVertices(fids[j], verts);
         for (int k = 0; k < verts.Size(); k++)
         {
            int index = verts[k];
            if (pb->get_indicator() == "saltzmann")
            {
               // Replace the bdr attribute in the array as long as it is not
               // the dirichlet condition (For Saltzman Problem)
               if (BdrVertexIndexingArray[index] != 1)
               {
                  BdrVertexIndexingArray[index] = bdr_attr;
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
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::MakeTimeStep(Vector &S, const double & t, double & dt)
{
   // cout << "========================================\n"
   //      << "MakeTimeStep\n"
   //      << "========================================\n";
   // chrono_mm.Start();
   if (mm)
   {
      // Compute mesh velocities
      // This function may change the timestep
      switch (mv_option)
      {
      case 1:
         ComputeMeshVelocitiesRaviart(S,t,dt);
         break;

      case 2:
         ComputeMeshVelocitiesNormal(S,t,dt);
         break;
      
      default:
         MFEM_ABORT("Invalid mesh velocity option.\n");
      }
      
      // Enforce boundary conditions
      EnforceMVBoundaryConditions(S,t,dt);
   }
   // chrono_mm.Stop();

   // Update state variables contained in S_new
   // chrono_state.Start();
   ComputeStateUpdate(S, t, dt);
   // chrono_state.Stop();

   // Move the mesh
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction x_gf, mv_gf;
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);
   mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);

   add(x_gf, dt, mv_gf, x_gf);

   pmesh->NodesUpdated();
   // cout << "mm computation took " << chrono_mm.RealTime() << "s.\n";
   // cout << "state update computation took " << chrono_state.RealTime() << "s.\n";
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
   // cout << "========================================\n"
   //      << "ComputeStateUpdate\n"
   //      << "========================================\n";

   // We need a place to store the new state variables
   Vector S_new = S;

   Vector val(dim+2), c(dim), n_int(dim), U_i(dim+2), U_j(dim+2), sums(dim+2);
   Array<int> fids, oris; 

   int cj = 0;
   double d;

   mfem::Mesh::FaceInformation FI;
   H1.ExchangeFaceNbrData();

   bool is_boundary_cell = false;

   Vector sum_validation(dim);
   for (int ci = 0; ci < NDofs_L2; ci++) // Cell iterator
   {
      is_boundary_cell = false;

      sum_validation = 0.;

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

      DenseMatrix F_i = pb->flux(U_i);
      sums = 0.;

      for (int j=0; j < fids.Size(); j++) // Face iterator
      {
         CalcOutwardNormalInt(S, ci, fids[j], n_int);
         c = n_int;
         c /= 2.;
         sum_validation.Add(1., c);

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
               d = dij_sparse->Elem(ci, cj); 

               Vector z = U_j;
               z -= U_i;
               sums.Add(d, z);
            }
            
         }
         else
         {
            assert(FI.IsBoundary());
            is_boundary_cell = true;
            Vector y_temp(dim+2), y_temp_bdry(dim+2), U_i_bdry(dim+2);

            F_i.Mult(c, y_temp);
            U_i_bdry = U_i;
            y_temp_bdry = 0.;

            /* Enforce Boundary Conditions */
            // if (pb->get_indicator() == "noh")
            // {
            //    // Since we will enforce Dirichlet conditions on all boundary cells, 
            //    // this enforcement will be done after the face iterator
            // }
            if (pb->get_indicator() == "saltzmann")
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
                  DenseMatrix F_i_slip = pb->flux(U_i_bdry);
                  F_i_slip.Mult(c, y_temp_bdry);
                  // y_temp *= 2.;
               }
               else
               {
                  cout << "invalid boundary attribute: " << bdr_attribute << endl;
                  cout << "cell : " << ci << endl;
                  cout << "face: " << fids[j] << endl;  
                  y_temp *= 2.; 
               }

               y_temp += y_temp_bdry;
            }
            else
            {
               y_temp *= 2.;
            }
            
            // Add in boundary contribution
            sums -= y_temp;
         } // End boundary face        

      } // End Face iterator

      // Enforce exact condition on boundary for Noh and Isentropic Vortex
      if (pb->has_boundary_conditions() && is_boundary_cell)
      {
         if (pb->get_indicator() == "IsentropicVortex") // || pb->get_indicator() == "Noh")
         {
            EnforceExactBCOnCell(S, ci, t, dt, val);
         }
         
      }
      else
      {
         // Do the normal thing
         assert(sum_validation.Norml2() < 1e-12);

         sums *= dt;
         double k = pmesh->GetElementVolume(ci);
         double _mass = k / U_i[0];
         sums /= _mass;
         val += sums;
      }

      // In either case, update the cell state vector
      SetCellStateVector(S_new, ci, val);
   } // End cell iterator

   // Finally, update S
   S = S_new;
}


/****************************************************************************************************
* Function: EnforceExactBCOnCell
* Parameters:
*  cell      - Index representing the cell
*  state_val - Returned vector containing the exacthydrodynamic state variables
*
* Purpose:
*  Given a cell, return the exact hydrodynamic state variables of that cell.
*  This function is used to enforce the exact solution on the boundary at a
*  particular cell.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::EnforceExactBCOnCell(const Vector &S, const int & cell, const double &t, 
                                                     const double &dt, Vector & state_val)
{
   // cout << "========================================\n"
   //      << "EnforceExactBCOnCell\n"
   //      << "========================================\n";

   // Compute cell center and corresponding velocity at this location
   // Get center node dof
   int cell_vdof;
   Vector cell_x(dim);

   switch (dim)
   {
      case 1:
      {
         cell_vdof = num_faces + cell;
         break;
      }
      case 2:
      {
         cell_vdof = NVDofs_H1 + num_faces + cell;
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

   GetNodePosition(S, cell_vdof, cell_x);

   // Fill in val with exact solution
   state_val[0] = pb->sv0(cell_x, t+dt);
   Vector v_exact(dim);
   pb->v0(cell_x, t+dt, v_exact);
   for (int i = 0; i < dim; i++)
   {
      state_val[1 + i] = v_exact[i];
   }
   state_val[dim + 1] = pb->ste0(cell_x, t+dt);
}


/****************************************************************************************************
* Function: EnforceMVBoundaryConditions
* Parameters:
*     S_new - BlockVector that stores mesh information, mesh velocity, and state variables.
*     t     - Current time
*     dt    - Current timestep
*
* Purpose:
*  Enforces boundary conditions on the mesh velocities according to the problem indicator.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::EnforceMVBoundaryConditions(Vector &S, const double &t, const double &dt)
{
   // cout << "========================================\n"
   //      << "EnforcingMVBoundaryConditions\n"
   //      << "========================================\n";

   if (pb->get_indicator() == "saltzmann")
   {
      // Enforce
      Vector* sptr = const_cast<Vector*>(&S);
      ParGridFunction x_gf, mv_gf;
      mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);

      Vector ex(dim);
      ex = 0.;
      /* Ramping up to ex */
      // if (timestep_first == 0.)
      // {
      //    timestep_first = timestep;
      // }
      // double _xi = t / (2*timestep_first);
      // double _psi = (4 - (_xi + 1) * (_xi - 2) * ((_xi - 2) - (abs(_xi-2) + (_xi-2)) / 2)) / 4.;
      // ex[0] = 1. * _psi;
      ex[0] = 1.;

      VectorConstantCoefficient left_wall_coeff(ex);
      Array<int> left_wall_attr(1);
      left_wall_attr[0] = 1;

      mv_gf.ProjectBdrCoefficient(left_wall_coeff, left_wall_attr);

      /* TODO: Vladimir 
         How to enforce v.n=0 on rest of boundary?
      */
      // Vector vZero(dim);
      // vZero = 0.;
      // VectorConstantCoefficient zero(vZero);
      // Array<int> other_bdr_attr(1);
      // other_bdr_attr[0] = 2;

      // mv_gf.ProjectBdrCoefficientNormal(zero, other_bdr_attr);

      mv_gf.SyncAliasMemory(S);
   }
   else if (pb->get_indicator() == "IsentropicVortex")
   {
      Vector* sptr = const_cast<Vector*>(&S);
      ParGridFunction x_gf, mv_gf;
      mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);

      Vector zero(dim);
      zero = 0.;
      VectorConstantCoefficient zero_coeff(zero);

      mv_gf.ProjectBdrCoefficient(zero_coeff, pmesh->bdr_attributes);
   }
}


/****************************************************************************************************
* Function: GetCellStateVector
* Parameters:
*  S    - BlockVector that stores mesh information, mesh velocity, and state variables.
*  cell - Index representing the cell
*  U    - Returned vector contained the hydrodynamic state variables
*
* Purpose:
*  Given a cell, return the hydrodynamic state variables of that cell, i.e. in vector for get
*  (specific volume, velocity, specific total energy)
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::GetCellStateVector(const Vector &S, const int cell, Vector &U)
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
* Function: SetCellStateVector
* Parameters:
*  S    - BlockVector that stores mesh information, mesh velocity, and state variables.
*  cell - Index representing the cell
*  U    - Vector containing state variables to be put into BlockVector object S.
*
* Purpose:
*  This function is used to update the BlockVector S with the newly computed hydrodynamic state
*  variables on the given cell.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::SetCellStateVector(Vector &S, const int cell, const Vector &U)
{
   Array<int> dofs;
   Array<int> sub_dofs;
   dofs.Append(cell);
   sub_dofs.Append(1);

   // Get current state grid functions
   Vector* sptr = const_cast<Vector*>(&S);
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
   sv_gf.SyncAliasMemory(S);
   v_gf.SyncAliasMemory(S);
   ste_gf.SyncAliasMemory(S);
}


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
   // cout << "=======================================\n"
   //      << "  Calculating outward normal integral  \n"
   //      << "=======================================\n";
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

         vertex_element->GetRow(face, row);

         int cell_gdof = cell + num_faces;
         GetNodePosition(S, cell_gdof, cell_center_x);
         GetNodePosition(S, face, face_x);

         cout << "cell: " << cell << ", location: ";
         cell_center_x.Print(cout);
         cout << "face: " << face << ", location: ";
         face_x.Print(cout);


         subtract(face_x, cell_center_x, res);
         res *= 1. / res.Norml2();

         cout << "outward normal: ";
         res.Print(cout);
         break;
      }
      case 2:
      {
         H1.GetFaceDofs(face, row);

         int face_vdof1 = row[1], face_vdof2 = row[0];
         face_dof = row[2];

         Vector face_x(dim), vdof1_x(dim), vdof2_x(dim), vdof1_v(dim), vdof2_v(dim);
         GetNodePosition(S, face_dof, face_x);
         GetNodePosition(S, face_vdof1, vdof1_x);
         GetNodePosition(S, face_vdof2, vdof2_x);

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
               res *= -1.;
            }
         }

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


/****************************************************************************************************
* Function: ComputeIntermediateFaceVelocities
* Parameters:
*  S        - BlockVector representing FiniteElement information
*  t        - Current time
*  flag     - Flag used to indicate testing, can be "testing" or "NA"
*  test_vel - Velocity used for testing
*
* Purpose:
*  This function computes the intermediate face velocities as outlined in equation (5.7).  If
*  the flag is set to "testing" and a test_vel function is passed in, the velocity at the faces
*  will be the pointwise evaluation of the test_vel function.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::
   ComputeIntermediateFaceVelocities(const Vector &S,
                                        const double t,
                                        const string flag, // Default NA
                                        void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
   /***
    * Compute intermediate face velocities according to (5.7) 
    ***/
   // cout << "=======================================\n"
   //      << "   ComputeIntermediateFaceVelocities   \n"
   //      << "=======================================\n";
   
   mfem::Mesh::FaceInformation FI;
   int c, cp;
   Vector Uc(dim+2), Ucp(dim+2), n_int(dim), c_vec(dim), n_vec(dim), Vf(dim), Vf_flux(dim); 
   double d, c_norm, F;

   Array<int> row;

   for (int face = 0; face < num_faces; face++) // face iterator
   {
      // cout << "face: " << face << endl;
      Vf = 0.;
      FI = pmesh->GetFaceInformation(face);
      c = FI.element[0].index;
      cp = FI.element[1].index;

      GetCellStateVector(S, c, Uc);

      if (flag == "NA")
      {
         // Get normal, d, and |F|
         CalcOutwardNormalInt(S, c, face, n_int);
         n_vec = n_int;
         F = n_vec.Norml2();
         n_vec /= F;
         assert(1. - n_vec.Norml2() < 1e-12);

         if (FI.IsInterior())
         {
            GetCellStateVector(S, cp, Ucp);

            // Get max wave speed
            d = dij_sparse->Elem(c, cp);

            // Compute intermediate face velocity
            Vf = pb->velocity(Uc);
            Vf += pb->velocity(Ucp);
            Vf *= 0.5;

            double coeff = d * (Ucp[0] - Uc[0]) / F;
            Vf.Add(coeff, n_vec);

            // Switch orientation of Vf based on orientation of faces in mfem
            // Vf_flux = Vf;
            // Vf_flux *= F;
            // int ind = -1;
            // int ori = 0;

            // Array<int> element_face_row, element_face_oris;
            // pmesh->GetElementEdges(c, element_face_row, element_face_oris);
            // for (int i = 0; i < 4; i++) // 3DTODO
            // {
            //    if (element_face_row[i] == face)
            //    {
            //       ind = i;
            //       ori = element_face_oris[i];
            //    }
            // }
            // if ((ind == 0 && ori == 1)  || 
            //       (ind == 1 && ori == -1) ||
            //       (ind == 2 && ori == -1) ||
            //       (ind == 3 && ori == 1))
            // {
            //    cout << "flipping orientation of Vf for face: " << face << endl;
            //    Vf_flux *= -1.;
            // }
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

         GetNodePosition(S, face_dof, face_x);
         test_vel(face_x, t, Vf);
      }

      // Put face velocity into object (for state var update)
      for (int i = 0; i < dim; i++)
      {
         int index = face + i * num_faces;
         v_CR_gf[index] = Vf[i];
      }

      // Put face velocity into object [for MFEM implementation, has errors]
      // for (int i = 0; i < dim; i++)
      // {
      //    int index = face + i * num_faces;
      //    v_CR_gf_fluxes[index] = Vf_flux[i];
      // }

   } // End face iterator

}


/****************************************************************************************************
* Function: SetMVOption
* Parameters:
*  option - integer parameter indicating which mesh velocity to use.
*
* Purpose:
*  To set the mesh velocity motion computation type:
*     1) Raviart-Thomas reconstruction method.
*     2) Normal vector type
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::SetMVOption(const int & option)
{
   this->mv_option = option;
}


/****************************************************************************************************
* Function: GetIntermediateFaceVelocity
* Parameters:
*  face - Face index
*  vel  - V_{F_m^n} corresponding to the given face
*
* Purpose:
*  Return intermediate face velocity corresponding to the given face.
*  NOTE: This function requires that ComputeIntermediateFaceVelocities has been run.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::GetIntermediateFaceVelocity(const int & face, Vector & vel)
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
* Function: SetCorrectedFaceVelocity
* Parameters:
*  face - Face index
*  vel  - V_{b}^3 corresponding to the given face
*
* Purpose:
*  Sets corrected face velocity corresponding to the given face.
*  NOTE: This function requires that ComputeIntermediateFaceVelocities has been run,
*        though this value will change after the corrective velocities have been computed.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::SetCorrectedFaceVelocity(const int & face, const Vector & vel)
{
   assert(face < num_faces);
   // Retrieve face velocity from object
   for (int i = 0; i < dim; i++)
   {
      int index = face + i * num_faces;
      v_CR_gf_corrected[index] = vel[i];
   }
}


/****************************************************************************************************
* Function: GetCorrectedFaceVelocity
* Parameters:
*  face - Face index
*  vel  - V_{b}^3 corresponding to the given face
*
* Purpose:
*  Return corrected face velocity corresponding to the given face.
*  NOTE: This function requires that ComputeIntermediateFaceVelocities has been run,
*        though this value will change after the corrective velocities have been computed.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::GetCorrectedFaceVelocity(const int & face, Vector & vel)
{
   assert(face < num_faces);
   // Retrieve face velocity from object
   for (int i = 0; i < dim; i++)
   {
      int index = face + i * num_faces;
      vel[i] = v_CR_gf_corrected[index];
   }
}


/****************************************************************************************************
* Function: SetCorrectedFaceFlux
* Parameters:
*  face - Face index
*  flux - V_{b}^3 corresponding to the given face
*
* Purpose:
*  Sets corrected face flux corresponding to the given face.
*  NOTE: This function requires that ComputeIntermediateFaceVelocities has been run, and that 
*        that v_CR_gf_corrected has been filled (i.e. that ComputeCorrectiveFaceVelocities() 
*        has been run).
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::SetCorrectedFaceFlux(const int & face, const Vector & flux)
{
   assert(face < num_faces);
   // Retrieve face velocity from object
   for (int i = 0; i < dim; i++)
   {
      int index = face + i * num_faces;
      // double new_flux_comp = (v_CR_gf_fluxes[index] + flux[i]) / 2.;
      // v_CR_gf_fluxes[index] = new_flux_comp;
      v_CR_gf_fluxes[index] = flux[i];
   }
}


/****************************************************************************************************
* Function: GetCorrectedFaceFlux
* Parameters:
*  face - Face index
*  flux - V_{b}^3 corresponding to the given face
*
* Purpose:
*  Return corrected face flux corresponding to the corrected face velocity on the given face.
*  NOTE: This function requires that ComputeIntermediateFaceVelocities has been run, and that 
*        that v_CR_gf_corrected has been filled (i.e. that ComputeCorrectiveFaceVelocities() 
*        has been run).
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::GetCorrectedFaceFlux(const int & face, Vector & flux)
{
   assert(face < num_faces);
   // Retrieve face velocity from object
   for (int i = 0; i < dim; i++)
   {
      int index = face + i * num_faces;
      flux[i] = v_CR_gf_fluxes[index];
   }
}


/****************************************************************************************************
* Function: CalcMassLoss
* Parameters:
*  S - BlockVector corresponding to nth timestep that stores mesh information, mesh velocity, 
*      and state variables.
*
* Purpose:
*  This function is used when computing values for the convergence tables.  It does this by computing
*  the sum of the local mass loss at each element.
****************************************************************************************************/
template<int dim>
double LagrangianLOOperator<dim>::CalcMassLoss(const Vector &S)
{
   // cout << "=======================================\n"
   //      << "             CalcMassLoss              \n"
   //      << "=======================================\n";
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

   return num / denom;
}


/****************************************************************************************************
* Function: CheckMassConservation
* Parameters:
*  S - BlockVector corresponding to nth timestep that stores mesh information, mesh velocity, 
*      and state variables.
*
* Purpose:
*  This function calculates the percentage of cells in the mesh where mass has not been conserved.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::CheckMassConservation(const Vector &S, ParGridFunction & mc_gf)
{
   // cout << "=======================================\n"
   //      << "         CheckMassConservation         \n"
   //      << "=======================================\n";
   Vector U_i(dim + 2);
   int counter = 0;

   double current_mass = 0., current_mass_sum = 0.;
   double num = 0., denom = 0.;
   double temp_num = 0., temp_denom = 0., val = 0.;
   
   for (int ci = 0; ci < NDofs_L2; ci++)
   {
      const double m = m_hpv->Elem(ci);
      const double k = pmesh->GetElementVolume(ci);
      GetCellStateVector(S, ci, U_i);

      current_mass = k / U_i[0];
      temp_num = abs(current_mass - m);
      temp_denom = abs(m);

      num += temp_num;
      denom += temp_denom;
      current_mass_sum += current_mass;

      val = temp_num / temp_denom;

      if (val > pow(10, -8))
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
      // Fill corresponding cell to indicate graphically the local change in mass, if any
      mc_gf[ci] = val;
   }

   double cell_ratio = (double)counter / (double)NDofs_L2;

   cout << "Percentage of cells where mass conservation was broken: " << cell_ratio << endl;
   cout << "Initial mass sum: " << denom 
        << ", Current mass sum: " << current_mass_sum << endl;
   cout << "Mass Error: " << num / denom << endl;
}


/****************************************************************************************************
* Function: ComputeDeterminant
* Parameters:
*  C     - Dense matrix representing C_i^geo from section 5.
*  dt    - Timestep (assumes full timestep is given, not half step)
*  alpha - Smallest postivie value satisfying equation (5.12) in [1]
*
* Purpose:
*  This function computes the quantity given by equation (5.12) in [1].  Specifically, this function
*  computes the largest positive value alpha_i satisfying
*        det(alpha_i\mathbb{I} - \frac{dt}{2}C_i) = alpha_i^{d-1}
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeDeterminant(const DenseMatrix &C, const double &dt, double & alpha)
{
   // cout << "=====================\n";
   // cout << "Computing determinant\n";
   // cout << "=====================\n";
   double trace = C.Trace();
   double det = C.Det();

   double a = 1.;
   double b = -1. * (1. + (dt / 2.) * trace);
   double c = (pow(dt, 2) / 4.) * det;

   double pm = sqrt(pow(b,2) - 4. * a * c);

   double alpha1 = -1. * b + pm;
   alpha1 /= (2. * a);

   double alpha2 = -1. * b - pm;
   alpha2 /= (2. * a);

   alpha = std::max(alpha1, alpha2);

   if (alpha <= 0.)
   {
      cout << "Dense Matrix:\n";
      C.Print(cout);
      cout << "trace: " << trace << ", det: " << det << endl;
      cout << "dt: " << dt << endl;
      cout << "a: " << a << ", b: " << b << ", c: " << c << endl;
      cout << "d1: " << alpha1 << ", d2: " << alpha2 << endl;
      MFEM_ABORT("Alpha_i should be positive.\n");
   }
}


/***********************************************************************************************************
* Function: ComputeCorrectiveFaceVelocities
* Parameters:
*   S        - BlockVector representing FiniteElement information
*   t        - Current time
*   dt       - Timestep
*   flag     - Flag used to indicate testing
*   test_vel - Velocity used for testing
* Purpose:
*   This function computes the node velocities on the faces which
*   are designed to bubble in the direction of the normal vector 
*   to conserve mass locally.
***********************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::
   ComputeCorrectiveFaceVelocities(Vector &S, const double & t, const double & dt,
                                    const string flag, // Default NA
                                    void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
   assert(dim > 1); // No need to correct face velocities in dim=1
   // cout << "==============================================\n";
   // cout << "Computing corrective interior face velocities.\n";
   // cout << "==============================================\n";
   /* Parameters needed for face velocity calculations */
   mfem::Mesh::FaceInformation FI;
   Vector Vf(dim), n_int(dim), n_vec(dim), face_velocity(dim);
   Vector Uc(dim+2), Ucp(dim+2);
   Vector cell_c_v(dim), cell_cp_v(dim), cell_center_v(dim); // For computation of bmn
   Array<int> row;
   int face_dof = 0, face_vdof2 = 0, c = 0, cp = 0;

   double V3n, V3nperp;

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
      GetNodePosition(S, face_dof, face_x);
      GetNodePosition(S, face_vdof1, vdof1_x);
      GetNodePosition(S, face_vdof2, vdof2_x);

      // retrieve corner velocities
      GetNodeVelocity(S, face_vdof1, vdof1_v);
      GetNodeVelocity(S, face_vdof2, vdof2_v);

      if (FI.IsInterior())
      {
         // Get adjacent cell information
         c = FI.element[0].index;
         cp = FI.element[1].index;
         GetCellStateVector(S, c, Uc);
         GetCellStateVector(S, cp, Ucp);
         cell_c_v = pb->velocity(Uc);
         cell_cp_v = pb->velocity(Ucp);

         // Calculate outer normal
         CalcOutwardNormalInt(S, c, face, n_int);
         n_vec = n_int;
         double F = n_vec.Norml2();
         n_vec /= F;

         assert(1. - n_vec.Norml2() < 1e-12);

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

         // Compute D (A.4c)
         Vector n_vec_R(dim), temp_vec(dim), temp_vec_2(dim);
         n_vec_R = n_vec;
         Orthogonal(n_vec_R);
         subtract(vdof1_v, vdof2_v, temp_vec); // V1 - V2 = temp_vec
         subtract(vdof2_x_half, vdof1_x_half, temp_vec_2); // A2-A1
         Orthogonal(temp_vec_2);
         double D = dt * (temp_vec * n_vec_R) + 2. * (n_vec * temp_vec_2);

         // Compute c1 (A.4a)
         subtract(vdof2_v, vdof1_v, temp_vec); // only change temp_vec, since temp_vec_2 is same from D calculation (half step representation)
         double c1 = ( dt * (temp_vec * n_vec) + 2. * (temp_vec_2 * n_vec_R) ) / D; // TRYING SOMETHING HERE. WILL NEED CORRECTED

         // Compute c0 (A.4b)
         Vector n_vec_half(dim);
         subtract(vdof2_x_half, vdof1_x_half, n_vec_half);
         Orthogonal(n_vec_half);
         GetIntermediateFaceVelocity(face, Vf);
         
         double bmn = Vf * n_vec;
         bmn *= F;

         temp_vec = vdof1_x_half;
         Orthogonal(temp_vec); // A1R
         temp_vec_2 = vdof2_x_half;
         Orthogonal(temp_vec_2); // A2R
         double const1 = vdof1_v * temp_vec - vdof2_v * temp_vec_2; // V1*A1R - V2*A2R
         double const2 = vdof1_v * temp_vec_2 - vdof2_v * temp_vec; // V1*A2R - V2*A1R
         
         temp_vec = face_x;
         Orthogonal(temp_vec);
         subtract(vdof2_v, vdof1_v, temp_vec_2);
         double const3 = temp_vec_2 * temp_vec; // (V2 - V1) * a3nR
         double c2 = (3. / D) * (const1 / 2. + const2 / 6. + 2. * const3 / 3.);
         
         /**
          * Compute V3n_perp
          * This value we choose:
          *    Option 1) new node is places in the normal direction
          *              of the midpoint of the secant line.
          *    Option 2) Choose perpendicular component of the 
          *              Raviart-Thomas velocity.
          * */ 

         // Necessary quantity for either calculation
         Vector n_vec_perp(dim);
         n_vec_perp = n_vec_R;
         n_vec_perp *= -1.;

         /* Option 1 */ 
         subtract(vdof2_x_new, vdof1_x_new, temp_vec);
         subtract(face_x, vdof12_x_new, temp_vec_2);
         temp_vec_2.Add((c2+ (3. * bmn / D))*dt, n_vec);
         const1 = temp_vec * temp_vec_2; // numerator
         temp_vec_2 = n_vec_R;
         temp_vec_2 *= -1.;
         temp_vec_2.Add(c1, n_vec);
         const2 = temp_vec * temp_vec_2;
         const2 *= dt; // denominator
         V3nperp = -1. * const1 / const2;

         /* Option 2 */
         Vector rt_face_vel(dim);
         GetNodeVelocity(S, face_dof, rt_face_vel);
         V3nperp = rt_face_vel * n_vec_perp;

         // Compute V3n (5.11)
         V3n = c1 * V3nperp + c2 + 3. * bmn / D;

         // Add in normal and tangentional contributions
         face_velocity = 0.;
         face_velocity.Add(V3n, n_vec);
         face_velocity.Add(V3nperp, n_vec_perp); // 10/2

         // Check perpendicular only in the case where new 
         // nodal location is exactly above tangent location
         // Vector face_x_new(dim);
         // face_x_new = face_x;
         // face_x_new.Add(dt, face_velocity);
         // subtract(face_x_new, vdof12_x_new, temp_vec);
         // subtract(vdof2_x_new, vdof1_x_new, temp_vec_2);
         // if (abs(temp_vec * temp_vec_2) > pow(10, -12))
         // {
         //    MFEM_ABORT("vectors are not orthogonal!\n");
         // }
      } // End interior face

      // On the boundary, we just need the average of the adjacent 
      // vertices to account for the corrective corner velocities
      // NOTE that even in the test case, the corner velocities will
      // be modified and so prescribing the exact velocities on the
      // boundary faces will not preserve mass.
      else
      {
         assert(FI.IsBoundary());
         // TODO: Add in boundary conditions similar to corner vertex bcs
         for (int j = 0; j < dim; j++)
         {
            face_velocity[j] = (vdof1_v[j] + vdof2_v[j]) / 2.;
         }
        
      } // End boundary face

      // Lastly, put face velocity into gridfunction object
      if (face_velocity[0] != face_velocity[0] || face_velocity[1] != face_velocity[1])
      {
         MFEM_ABORT("NaN values encountered in corrective face velocity calculation.\n");
      }

      UpdateNodeVelocity(S, face_dof, face_velocity);
      SetCorrectedFaceVelocity(face, face_velocity);
   }
}


/***********************************************************************************************************
* Function: ComputeCorrectiveFaceFluxes
* Parameters:
*   S        - BlockVector representing FiniteElement information
*   t        - Current time
*   dt       - Timestep
* Purpose:
*   This function computes the geometric face node velocities which
*   are used in the iterative procedure to recompute corner node 
*   velocities.
***********************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::
   ComputeCorrectiveFaceFluxes(Vector &S, const double & t, const double & dt)
{
   assert(dim > 1); // There is no need for this function in dim=1
   cout << "========================================================\n";
   cout << "Computing corrective interior geometric face velocities.\n";
   cout << "========================================================\n";
   /* Parameters needed for face velocity calculations */
   mfem::Mesh::FaceInformation FI;
   Vector Vf(dim), n_int(dim), n_vec(dim), face_velocity(dim);
   Vector cell_c_v(dim), cell_cp_v(dim), cell_center_v(dim); // For computation of bmn
   Array<int> row;
   int face_dof = 0, face_vdof2 = 0, c = 0, cp = 0;

   double V3n, V3nperp;

   // Iterate over faces
   for (int face = 0; face < num_faces; face++) // face iterator
   {  
      // Get intermediate face velocity, face information, and face normal
      face_velocity = 0.;
      Vf = 0.;

      FI = pmesh->GetFaceInformation(face);
      GetCorrectedFaceVelocity(face, face_velocity);

      /* adjacent corner indices */
      H1.GetFaceDofs(face, row);
      int face_dof = row[2];

      if (FI.IsInterior())
      {
         // Get adjacent cell information
         c = FI.element[0].index;
         cp = FI.element[1].index;

         // Calculate outer normal
         CalcOutwardNormalInt(S, c, face, n_int);
         n_vec = n_int;
         double F = n_vec.Norml2();
         n_vec /= F;

         assert(1. - n_vec.Norml2() < 1e-12);

         /*** Compute Face Flux ***/
         cout << "computing face flux for face: " << face << ", which corresponds to face_dof: " << face_dof << endl;

         // Compute C_i
         DenseMatrix Ci;
         // mfem::Array<int> elements;
         // face_element->GetRow(face_dof, elements);
         ComputeCiGeoRaviart(face_dof, Ci);
         cout << "ci:\n";
         Ci.Print(cout);

         // Compute alpha_i
         double alpha_i;
         ComputeDeterminant(Ci, dt, alpha_i);
         cout << "alpha_i: " << alpha_i << endl;

         // Compute flux
         Vector face_flux(dim);
         DenseMatrix _mat(Ci);
         _mat *= - dt / 2.;
         for (int i = 0; i < dim; i++)
         {
            _mat(i,i) += alpha_i;
         }

         cout << "face flux correction matrix:\n";
         _mat.Print(cout);

         _mat.Mult(face_velocity, face_flux);
         SetCorrectedFaceFlux(face, face_flux); 
      }
   }
}


/****************************************************************************************************
* Function: FillCenterVelocitiesWithAvg
* Parameters:
*   S        - BlockVector representing FiniteElement information
*
* Purpose:
*  This function is only needed since we are not able to implement the Serendipity finite elements.
*  Since the mesh motion does not depends on a value for the node corresponding to the center of the 
*  element, we need to fill this with a velocity which will ensure the center node remains inside the
*  cell.  This is done by average the values at the four corner nodes of the cell.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::FillCenterVelocitiesWithAvg(Vector &S)
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
            GetNodeVelocity(S, verts[j], node_v);
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
            // GetIntermediateFaceVelocity(face, node_v);
            GetNodeVelocity(S, face + NVDofs_H1, node_v);
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
      UpdateNodeVelocity(S, cell_vdof, Vc);
   }
}


/****************************************************************************************************
* Function: FillFaceVelocitiesWithAvg
* Parameters:
*   S        - BlockVector representing FiniteElement information
*   flag     - Flag used to indicate testing
*   test_vel - Velocity used for testing
*
* Purpose:
*  This function is merely for testing purposes and should not be implemented in the full program
*  as the purpose of moving the face nodes is to ensure local mass conservation.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::FillFaceVelocitiesWithAvg(Vector &S, const string flag, void (*test_vel)(const Vector&, const double&, Vector&))
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
      GetNodePosition(S, face_dof, face_x);

      // retrieve corner velocities
      GetNodeVelocity(S, face_vdof1, vdof1_v);
      GetNodeVelocity(S, face_vdof2, vdof2_v);
      

      // Average nodal velocities
      for (int j = 0; j < dim; j++)
      {
         face_velocity[j] = (vdof1_v[j] + vdof2_v[j]) / 2;
      }

      // test_vel(face_x, 0, face_velocity);

      // Lastly, put face velocity into gridfunction object
      UpdateNodeVelocity(S, face_dof, face_velocity);
   }
}


template<int dim>
void LagrangianLOOperator<dim>::FillFaceVelocitiesWithAvg2(Vector &S, const string flag, void (*test_vel)(const Vector&, const double&, Vector&))
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
      GetNodePosition(S, face_dof, face_x);

      // retrieve corner velocities
      GetNodeVelocity(S, face_vdof1, vdof1_v);
      GetNodeVelocity(S, face_vdof2, vdof2_v);
      

      // Average nodal velocities
      for (int j = 0; j < dim; j++)
      {
         face_velocity[j] = (vdof1_v[j] + vdof2_v[j]) / 2;
      }

      // test_vel(face_x, 0, face_velocity);

      // Lastly, put face velocity into gridfunction object
      UpdateNodeVelocity(S, face_dof, face_velocity);
      if (FI.IsInterior())
      {
         SetCorrectedFaceVelocity(face, face_velocity);
      }
      
   }
}


template<int dim>
void LagrangianLOOperator<dim>::UpdateFaceVelocitiesWithAvg(Vector &S)
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
      GetNodePosition(S, face_dof, face_x);

      // retrieve corner velocities
      GetNodeVelocity(S, face_vdof1, vdof1_v);
      GetNodeVelocity(S, face_vdof2, vdof2_v);
      

      // Average nodal velocities
      for (int j = 0; j < dim; j++)
      {
         face_velocity[j] = (vdof1_v[j] + vdof2_v[j]) / 2;
      }

      // test_vel(face_x, 0, face_velocity);

      // Lastly, put face velocity into gridfunction object
      Vector corrective_face_velocity(dim), _temp_vel(dim);
      GetCorrectedFaceVelocity(face, corrective_face_velocity);
      double _theta = 0.;
      _temp_vel.Add(_theta, corrective_face_velocity);
      _temp_vel.Add(1. - _theta, face_velocity);

      UpdateNodeVelocity(S, face_dof, _temp_vel);
   }
}

/****************************************************************************************************
* Function: UpdateNodeVelocity
* Parameters:
*   S    - BlockVector representing FiniteElement information
*   node - Global index of node in question.
*   vel  - Velocity computed at node to be used to move the mesh
* Purpose:
*   This function is used to update the mv_gf ParGridFunction which is used to move the mesh.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::UpdateNodeVelocity(Vector &S, const int & node, const Vector & vel)
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


/****************************************************************************************************
* Function: GetNodeVelocity
* Parameters:
*   S    - BlockVector representing FiniteElement information
*   node - Global index of node in question.
*   vel  - Velocity at given node
*
* Purpose:
*  This function returns the velocity at the given global node.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::GetNodeVelocity(const Vector &S, const int & node, Vector & vel)
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
* Function: GetViGeo
* Parameters:
*   node - Global index of node in question.
*   vel  - Velocity Vi^geo at given node
*
* Purpose:
*  This function returns Vi_geo.
*  NOTE: Function LagrangianLOOperator::ComputeGeoVRaviart must be called first.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::GetViGeo(const int & node, Vector & vel)
{
   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      vel[i] = v_geo_gf[index];
   }
}


/****************************************************************************************************
* Function: GetNodePosition
* Parameters:
*   S    - BlockVector representing FiniteElement information
*   node - Global index of node in question.
*   x    - Coordinate position at given node.
*
* Purpose:
*  This function returns the cartesian location corresponding to a global node.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::GetNodePosition(const Vector &S, const int & node, Vector & x)
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


/****************************************************************************************************
* Function: SaveStateVecsToFile
* Parameters:
*  S                  - BlockVector representing FiniteElement information
*  output_file_prefix - Parameters used to define output file location
*  output_file_suffix - See above.
*
* Purpose:
*  The results of this program can be used to plot the same problem at multiple refinements.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::SaveStateVecsToFile(const Vector &S, 
                                                    const string &output_file_prefix, 
                                                    const string &output_file_suffix)
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


/****************************************************************************************************
* Function: Normal direction based Velocity Functions
****************************************************************************************************/

/****************************************************************************************************
* Function: tensor
* Parameters:
*  
*        
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
* Function: ComputeMeshVelocitiesNormal
* Parameters:
*  S        - BlockVector corresponding to nth timestep that stores mesh information, 
*             mesh velocity, and state variables.
*  t        - Current time
*  dt       - Current timestep
*  flag     - Flag used to indicate testing, can be "testing" or "NA"
*  test_vel - Velocity used for testing
*
* Purpose:
*     This function computes mesh velocities by:
*        1) Constructs Mi, which is the sum of the tensor products of all
*           adjacent normal vectors.
*        2) Constructs Ri, which is the sum of the tensor products of all
*           adjacent normal vectors, multiplied by the intermediate face 
*           velocity.
*        3) Computes Vi as Mi^{-1} * Ri.
*        4) Optionally computes mass corrective face velocities.
*        5) Fills cell center velocities with average.
*        
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeMeshVelocitiesNormal(
   Vector &S, const double & t, double & dt, const string flag, 
   void (*test_vel)(const Vector&, const double&, Vector&))
{
   // cout << "=======================================\n"
   //      << "     ComputeMeshVelocitiesNormal       \n"
   //      << "=======================================\n";
   ComputeIntermediateFaceVelocities(S, t, flag, test_vel);
   ComputeNodeVelocitiesNormal(S, t, dt, flag, test_vel);

   // Optionally mass correct
   if (dim > 1)
   {
      if (do_mass_correction)
      {
         // cout << "Mass correcting\n";
         ComputeCorrectiveFaceVelocities(S, t, dt, flag, test_vel);
      }
      // else
      // {
      //    FillFaceVelocitiesWithAvg(S);
      // }
   }

   // Fill cell centers with average
   FillCenterVelocitiesWithAvg(S);
}


/****************************************************************************************************
* Function: 
* Parameters:
*  S    - 
*
* Purpose:
*  
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeNodeVelocitiesNormal(
   Vector &S, const double & t, double & dt, const string, 
   void (*test_vel)(const Vector&, const double&, Vector&))
{
   // Compute mesh node velocities
   mfem::Mesh::FaceInformation FI;
   Vector node_v(dim), Ri(dim), n_vec(dim), Vf(dim), y(dim);
   DenseMatrix Mi(dim), dm_tmp(dim);

   Array<int> faces_row, oris;
   int faces_length;
   double c_norm;

   // Corner nodes
   for (int node = 0; node < num_vertices; node++) // Vertex iterator
   {
      // Reset Mi, Ri, Vi
      Mi = 0., Ri = 0., node_v = 0.;

      // Get cell faces
      vertex_edge.GetRow(node, faces_row);
      faces_length = faces_row.Size();

      // iterate over adjacent faces
      for (int face_it = 0; face_it < faces_length; face_it++) // Adjacent face iterator
      {
         int face = faces_row[face_it];
         GetIntermediateFaceVelocity(face, Vf);
         FI = pmesh->GetFaceInformation(face);

         // Retrieve corresponding normal (orientation doesn't matter)
         CalcOutwardNormalInt(S, FI.element[0].index, face, n_vec);
         c_norm = n_vec.Norml2();
         n_vec /= c_norm;

         tensor(n_vec, n_vec, dm_tmp);
         Mi += dm_tmp;
         dm_tmp.Mult(Vf, y);
         Ri += y;
      }

      Mi.Invert();
      Mi.Mult(Ri, node_v);

      // Finally, update the node velocity
      UpdateNodeVelocity(S, node, node_v);
   }

   // face nodes just set to ivf
   for (int face = 0; face < num_faces; face++)
   {
      int node = face + num_vertices;
      GetIntermediateFaceVelocity(face, node_v);
      UpdateNodeVelocity(S, node, node_v);
   }
}


/****************************************************************************************************
* Function: Raviart-Thomas Velocity Functions
****************************************************************************************************/

/****************************************************************************************************
* Function: ComputeMeshVelocitiesRaviart
* Parameters:
*  S        - BlockVector corresponding to nth timestep that stores mesh information, 
*             mesh velocity, and state variables.
*  t        - Current time
*  dt       - Current timestep
*  flag     - Flag used to indicate testing, can be "testing" or "NA"
*  test_vel - Velocity used for testing
*
* Purpose:
*     This function accomplishes the following steps:
*        1) Constructs intermediate face velocity objects for all faces at the nth timestep.
*        2) Computes V_i^geo on all nodes using Raviart-Thomas
*        3) Compute the corrected velocity, which uses the above defined velocity field to 
*           then compute C^i on all nodes
*        4) Computes the corrective face velocities.
*        5) Fills unnecessary node at cell center with the average velocity of the corner nodes.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeMeshVelocitiesRaviart(
   Vector &S, const double & t, double & dt, const string flag, 
   void (*test_vel)(const Vector&, const double&, Vector&))
{
   // cout << "=======================================\n"
   //      << "     ComputeMeshVelocitiesRaviart      \n"
   //      << "=======================================\n";

   ComputeIntermediateFaceVelocities(S, t, flag, test_vel);
   v_CR_gf_fluxes = v_CR_gf; // V_F^0 = V_F
   
   for (int face_corr_it = 0; 
        face_corr_it < num_face_correction_iterations + 1;  // Option set in Laglos with -fci
        face_corr_it++)
   {
      StopWatch _chrono;
      if (dim > 1)
      {
         // _chrono.Clear();
         // _chrono.Start();
         ComputeGeoVRaviart(S);
         // ComputeGeoVRaviart2(S);
         // _chrono.Stop();
         // cout << "computegeov took: " << _chrono.RealTime() << " seconds.\n";
      }

      // _chrono.Clear();
      // _chrono.Start();
      ComputeNodeVelocitiesRaviart(S, t, dt);
      // _chrono.Stop();
      // cout << "linearization took: " << _chrono.RealTime() << " seconds.\n";

      if (dim > 1)
      {
         if (do_mass_correction)
         {
            // cout << "Mass correcting\n";
            ComputeCorrectiveFaceVelocities(S, t, dt, flag, test_vel);
            
            // Turn off mass correction if at the end of the face
            // corner node correction iteration
            if (face_corr_it >= num_face_correction_iterations)
            {
               do_mass_correction = false;
            }
            else
            {
               // Since we will continue to iterate, we must
               // compute our iterative face flux and store 
               // the value
               ComputeCorrectiveFaceFluxes(S, t, dt);
            }
            // Otherwise keep mass correcting
            // do_mass_correction = false;
         }
         // else
         // {
         //    cout << "Face averaging\n";
         //    // ComputeCorrectiveFaceVelocities(S,t,dt,flag,test_vel);
         //    // UpdateFaceVelocitiesWithAvg(S); // theta parameter to incremently add newly computed
         //    FillFaceVelocitiesWithAvg(S);

         //    // End the iteration loop
         //    face_corr_it = num_face_correction_iterations + 1;
         // }
      }
   } // End face corner node correction iteration

   FillCenterVelocitiesWithAvg(S);
}


/****************************************************************************************************
* Function: 
* Parameters:
*  S    - 
*
* Purpose:
*  Iterate over nodes, then look at adjacent elements
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeGeoVRaviart(const Vector &S)
{
   // cout << "=======================================\n"
   //      << "          ComputeGeoVRaviart           \n"
   //      << "=======================================\n";

   Vector node_v(dim), n_vec(dim);
   Vector ivf(dim), v_cell(dim);
   Array<int> element_row, cell_faces_row;
   int row_length, cell_faces_length;

   DofEntity entity;
   int EDof;

   for (int node = 0; node < NDofs_H1; node++) // Vertex iterator
   {
      node_v = 0.; // Reset node velocity

      // Check if node corresponds to a face or to a 
      // vertex and get adjacent cell indices and edge indices
      GetEntityDof(node, entity, EDof);

      switch (entity)
      {
         case 0: // corner
         {
            // cout << "ComputeGeoVRaviart::corner node\n";
            vertex_element->GetRow(node, element_row);
            row_length = element_row.Size();

            break;
         }
         case 1: // face
         {
            // cout << "ComputeGeoVRaviart::face node for face: " << EDof << endl;
            face_element->GetRow(EDof, element_row);
            row_length = element_row.Size();
            // cout << "row length: " << row_length << endl;

            break;
         }
         case 2: // Cell Center
         {
            Vector Uc(dim+2);
            GetCellStateVector(S, EDof, Uc);
            node_v = pb->velocity(Uc);

            /* Put this velocity into v_geo_gf*/
            for (int i = 0; i < dim; i++)
            {
               int index = node + i * NDofs_H1;
               v_geo_gf[index] = node_v[i];
            }
            continue;
         }
         default:
         {
            MFEM_ABORT("Invalid entity\n");
         }
      }

      // Once you have the adjacent elements to a node, the computation is the same
      // Iterate over cells that share this node and average their contribution to vertex velocity
      for (int element_it = 0; element_it < row_length; element_it++)
      {
         v_cell = 0.;
         int el_index = element_row[element_it];
         // cout << "\t node: " << node << endl;

         Array<int> oris;
         pmesh->GetElementEdges(el_index, cell_faces_row, oris);
         cell_faces_length = cell_faces_row.Size();

         // Get Transformation and IntegrationPoint
         ElementTransformation * trans = pmesh->GetElementTransformation(el_index);
         Vector node_x(dim);
         GetNodePosition(S, node, node_x);
         IntegrationPoint ip;
         trans->TransformBack(node_x, ip);
         trans->SetIntPoint(&ip);

         // Get vector shape function evaluations
         DenseMatrix shape_funcs(4,2);
         const FiniteElement * CR_fe = CR.GetFE(el_index);
         CR_fe->CalcPhysVShape(*trans, shape_funcs);

         // Iterate over faces of this cell that touch vertex
         for (int face_it = 0; face_it < cell_faces_length; face_it++)
         {
            ivf = 0.;
            int face_index = cell_faces_row[face_it];

            GetCorrectedFaceFlux(face_index, ivf); 
            CalcOutwardNormalInt(S, el_index, face_index, n_vec);

            double coeff = n_vec * ivf;

            Vector _shape_fun(dim);
            shape_funcs.GetRow(face_it, _shape_fun);

            v_cell.Add(coeff, _shape_fun);
         }

         node_v.Add(1., v_cell);
      }
      // Average contributions
      node_v /= row_length;

      /* Put this velocity into v_geo_gf*/
      for (int i = 0; i < dim; i++)
      {
         int index = node + i * NDofs_H1;
         v_geo_gf[index] = node_v[i];
      }
   } // End vertex iterator
}


/****************************************************************************************************
* Function: 
* Parameters:
*  S    - 
*
* Purpose:
*  Iterate over cells and average contribution to nodes
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeGeoVRaviart2(const Vector &S)
{
   // cout << "=======================================\n"
   //      << "          ComputeGeoVRaviart2          \n"
   //      << "=======================================\n";
   Vector n_vec(dim), ivf(dim), v_cell(dim), node_x(dim), _shape_fun(dim);
   Array<int> cell_faces_row, oris, zones_per_dof, vdofs;
   int row_length, cell_faces_length;

   // reset Vgeo
   v_geo_gf = 0.;

   // Element transformation specific
   ElementTransformation * trans;
   IntegrationPoint ip;
   int ldof;
   DenseMatrix shape_funcs(4,2);

   zones_per_dof.SetSize(NDofs_H1);
   zones_per_dof = 0;

   // Iterate over all cells and distribute its contribution to the nodes
   for (int el = 0; el < L2.GetNE(); el++)
   {
      pmesh->GetElementEdges(el, cell_faces_row, oris);
      cell_faces_length = cell_faces_row.Size();
      
      trans = pmesh->GetElementTransformation(el);
      H1.GetElementDofs(el, vdofs);

      for (int i = 0; i < vdofs.Size(); i++)
      {
         v_cell = 0.;
         ldof = vdofs[i];

         // Get integration point corresponding to dof
         GetNodePosition(S, ldof, node_x);
         trans->TransformBack(node_x, ip);
         trans->SetIntPoint(&ip);

         // Get vector shape function evaluations
         const FiniteElement * CR_fe = CR.GetFE(el);
         CR_fe->CalcPhysVShape(*trans, shape_funcs);

         // Iterate over faces of this cell
         for (int face_it = 0; face_it < cell_faces_length; face_it++)
         {
            ivf = 0.;
            int face_index = cell_faces_row[face_it];

            GetCorrectedFaceFlux(face_index, ivf); 
            CalcOutwardNormalInt(S, el, face_index, n_vec);

            double coeff = n_vec * ivf;
            shape_funcs.GetRow(face_it, _shape_fun);
            v_cell.Add(coeff, _shape_fun);
         }

         /* Put this velocity into v_geo_gf*/
         for (int i = 0; i < dim; i++)
         {
            int index = ldof + i * NDofs_H1;
            v_geo_gf[index] += v_cell[i];
         }
         
         // Keep a tally per dof of contributions for later averaging
         zones_per_dof[ldof]++;
      }
   }

   // Average contributions from all cells
   for (int i = 0; i < zones_per_dof.Size(); i++)
   {
      const int nz = zones_per_dof[i];
      if (nz) { v_geo_gf(i) /= nz; }
   }
}

/****************************************************************************************************
* Function: 
* Parameters:
*  S    - 
*
* Purpose:
*  
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeNodeVelocitiesRaviart(
   Vector &S, const double & t, double & dt, const string, 
   void (*test_vel)(const Vector&, const double&, Vector&))
{
   // cout << "=======================================\n"
   //      << "     ComputeNodeVelocitiesRaviart      \n"
   //      << "=======================================\n";
   
   Vector node_v(dim);

   bool is_dt_changed = false;
   // Iterate over vertices
   for (int node = 0; node < NDofs_H1 - NDofs_L2; node++) // Vertex iterator
   {
      ComputeNodeVelocityRaviart(node, dt, node_v, is_dt_changed);
      // GetViGeo(node, node_v);

      // if (node_v[0] != node_v[0] || node_v[1] != node_v[1])
      // {
      //    cout << "NaN velocity encountered in ComputeNodeVelocities at node: " << node << endl;
      //    cout << "node_v in violation: ";
      //    node_v.Print(cout);
      //    MFEM_ABORT("Aborting due to NaNs.\n");
      // }

      UpdateNodeVelocity(S, node, node_v);

      // If we restricted the timestep, we must recompute the vertex velocities that were computed previously
      // if (is_dt_changed)
      // {
      //    vertex = -1;
      //    is_dt_changed = false;
      //    cout << "Restarting vertex iterator\n";
      // }
   } // End Vertex iterator
}


/****************************************************************************************************
* Function: 
* Parameters:
*  S    - 
*
* Purpose:
*  
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeNodeVelocityRaviart(
   const int & node, double & dt, Vector &node_v, bool &is_dt_changed)
{
   // cout << "=======================================\n"
   //      << "      ComputeNodeVelocityRaviart       \n"
   //      << "=======================================\n";
   chrono_temp.Clear();
   switch (dim)
   {
      case 1:
      {
         // Remark 5.3 in the paper states that since each geometric node belongs to one
         // face only, the RT velocity is equal to the IFV previously computed.
         assert(node < num_faces);
         GetIntermediateFaceVelocity(node, node_v);
         break;
      }
      default:
      {
         DenseMatrix Ci(dim);
         Vector Vgeo(dim);
         double d = 0.;
         
         Ci = 0., Vgeo = 0., node_v = 0.;
         // mfem::Array<int> elements;
         // vertex_element->GetRow(node, elements);
         // chrono_temp.Clear();
         // chrono_temp.Start();
         ComputeCiGeoRaviart(node, Ci);
         // chrono_temp.Stop();
         // cout << "ci computation for node " << node << " took: " << chrono_temp.RealTime() << "s\n";


         // Enforce time restriction imposed by calculation of alpha_i
         // if (dim == 2)
         // {
         //    double trace = Ci.Trace();
         //    double det = Ci.Det();
         //    double val = 2. * sqrt(det);
         //    cout << "trace: " << trace << ", det: " << det << endl;
         //
         //    // alpha_i must be positive restriction
         //    if (det <= 1.e-12 && trace < 0.)
         //    {
         //       if (dt > 2. / abs(trace))
         //       {
         //          dt = 2. / (abs(trace) + 1.);
         //          is_dt_changed = true;
         //          cout << "timestep restriction from velocity computation, det = 0 case.\n";
         //       }
         //    }
         //  
         //    // alpha_i must be real
         //    else if (det > 0.)
         //    {
         //       if (trace < val)
         //       {
         //          // Enforce timestep restriction
         //          double zero2 = 1. / (-trace + val);
         //          if (dt > 2 * zero2)
         //          {
         //             cout << "timestep restriction from velocity computation\n";
         //             cout << "dt: " << dt << ", z2: " << zero2 << endl;
         //             dt = 2. * zero2;
         //             is_dt_changed = true;
         //          }
         //       }
         //       else 
         //       {
         //          cout << "positive determinant but no timestep restriction needed.\n";
         //       }
         //    }
         //    else
         //    {
         //       cout << "negative determinant, no timestep restriction needed.\n";
         //    }
         // } // End time restriction from velocity computation

         GetViGeo(node, Vgeo);
         
         // chrono_temp.Clear();
         // chrono_temp.Start();
         ComputeDeterminant(Ci, dt, d);
         // chrono_temp.Stop();
         // cout << "determinant for node " << node << " took: " << chrono_temp.RealTime() << "s\n";

         // Compute V_i^n
         // chrono_temp.Clear();
         // chrono_temp.Start();
         DenseMatrix _mat(Ci);
         _mat *= - dt / 2.;
         for (int i = 0; i < dim; i++)
         {
            _mat(i,i) += d;
         }

         _mat.Invert();
         _mat.Mult(Vgeo, node_v);
         // chrono_temp.Stop();
         // cout << "mult for node " << node << " took: " << chrono_temp.RealTime() << "s\n";

         break;
      }
   }
}


/****************************************************************************************************
* Function: 
* Parameters:
*  S    - 
*
* Purpose:
*  
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeCiGeoRaviart(const int & node, DenseMatrix & res)
{
   // cout << "=======================================\n"
   //      << "         ComputeCiGeoRaviart           \n"
   //      << "=======================================\n";
   assert(node < NDofs_H1); // "Invalid nodal index"
   res.SetSize(dim);
   res = 0.;

   // Check if node corresponds to a face or to a 
   // vertex and get adjacent cell indices
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
         // cout << "ComputeCiGeoRaviart::corner node\n";
         vertex_element->GetRow(node, row);
         break;
      }
      case 1: // face
      {
         // cout << "ComputeCiGeoRaviart::face node\n";
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

   // cout << "computing cigeo for entity: " << entity << ", EDof: " << EDof << endl;
   // cout << "corresponding node dof: " << node << endl;
   // cout << "row: ";
   // row.Print(cout);
   int row_length = row.Size();
   for (int row_it = 0; row_it < row_length; row_it++)
   {
      int row_el = row[row_it];
      double cell_vol = pmesh->GetElementVolume(row_el);
      denom += cell_vol;

      dm_temp = 0.;
      IntGradRaviart(row_el, dm_temp);
      res.Add(cell_vol, dm_temp);
   }
   res *= 1./denom;
}


/****************************************************************************************************
* Function: IntGradRaviart
* Parameters:
*    cell - index corrseponding to the cell (K_c)
*    res  - (dim x dim) DenseMatrix representing the outer product of Vfm with 
*           the gradient of its corresponding scalar RT shape function.
* Purpose:
*    This function calculates the integral of the gradient of the v_geo velocity function
*    on a given cell.

* NOTE: This function assumes that dim > 1.
* NOTE: This function assumes that the functions:
*          - LagrangianLOOperator::ComputeIntermediateFaceVelocities()
*          - LagrangianLOOperator::ComputeGeoVRaviart() 
*       have already been called.  If the functions have not been called, then the 
*       returned velocity will be 0.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::IntGradRaviart(const int cell, DenseMatrix & res)
{
   // cout << "=======================================\n"
   //      << "             IntGradRaviart            \n"
   //      << "=======================================\n";

   ParGridFunction H1c_gf(&H1c);
   const int _size = H1c.GetVSize();

   ElementTransformation * trans = pmesh->GetElementTransformation(cell);
   Vector grad(dim), row(dim);
   res.SetSize(dim);
   res = 0., row = 0.;

   // Iterate over quadrature
   for (int i = 0; i < RT_ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = RT_ir.IntPoint(i);
      trans->SetIntPoint(&ip);

      for (int j = 0; j < dim; j++)
      {
         H1c_gf.MakeRef(&H1c, v_geo_gf, j*_size);
         H1c_gf.GetGradient(*trans, grad);
         // cout << "grad for dim " << j << ": \n";
         // grad.Print(cout);
         // cout << "H1c_gf: \n";
         // H1c_gf.Print(cout);

         // Put information into Dense Matrix
         res.GetRow(j, row);
         row.Add(ip.weight, grad);
         res.SetRow(j, row);
      }
   }
}


/* Explicit n_vectantiation */
template class LagrangianLOOperator<1>;
template class LagrangianLOOperator<2>;
template class LagrangianLOOperator<3>;

} // end ns hydrodynamics

} // end ns mfem
