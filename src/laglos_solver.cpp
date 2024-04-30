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
   H1_L(h1_l),
   H1Lc(H1_L.GetParMesh(), H1_L.FEColl(), 1),
   L2(l2),
   L2V(l2v),
   CR(cr),
   CRc(CR.GetParMesh(), CR.FEColl(), 1),
   v_CR_gf(&CR),
   v_CR_gf_corrected(&CR), 
   v_CR_gf_fluxes(&CR),
   v_geo_gf(&H1_L),
   cell_bdr_flag_gf(&L2),
   mv_gf_prev_it(&h1),
   pmesh(H1.GetParMesh()),
   m_lf(m),
   pb(_pb),
   Vsize_H1(H1.GetVSize()),
   TVSize_H1(H1.TrueVSize()),
   GTVSize_H1(H1.GlobalTrueVSize()),
   NDofs_H1(H1.GetNDofs()),
   NDofs_H1L(H1_L.GetNDofs()),
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

   // Fill cell boundaries
   cell_bdr_flag_gf = -1.;
   FillCellBdrFlag();

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
   cout << "NDofs_H1L: " << NDofs_H1L << endl;
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
   delete dij_sparse;
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
         double pl = pb->pressure(Uc,pmesh->GetAttribute(c));
         double pr = pb->pressure(Ucp,pmesh->GetAttribute(cp));

         // Finally compute lambda max
         // cout << "pre compute lambda max\n";
         lambda_max = pb->compute_lambda_max(Uc, Ucp, n_vec, pl, pr, pb->get_b());
         d = lambda_max * c_norm; 

         double ss = pb->sound_speed(Uc, pmesh->GetAttribute(c));

         dij_sparse->Elem(c,cp) = d;
         dij_sparse->Elem(cp,c) = d;

         if (dim == 1)
         {
            lambda_max_vec[face] = lambda_max; // TODO: remove, only temporary
         }
      }
   } // End face iterator
   // cout << "Printing dij_sparse:\n";
   // dij_sparse->Print(cout);
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

      GetCellStateVector(S, ci, U_i);

      // Compute mass at time tn
      const double k = pmesh->GetElementVolume(ci);
      mi = k / U_i[0];

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
*  BdrElementIndexingArray[face] = bdr_attr, and if the face is an interior face, then we will have
*  BdrElementIndexingArray[face] = -1.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::CreateBdrElementIndexingArray()
{
   cout << "Constructing BdrElementIndexingArray:\n";
   for (int i = 0; i < pmesh->GetNBE(); i++)
   {
      int bdr_attr = pmesh->GetBdrAttribute(i);
      int index = pmesh->GetBdrElementEdgeIndex(i);
      BdrElementIndexingArray[index] = bdr_attr;
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
*
*  Note that this function will label a vertex with 5 if that vertex lies on the corner, indicating 
*  that the vertex velocity should be set to 0 to preserve the slip BCs on both of its faces.
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
      int face = pmesh->GetBdrFace(i);

      pmesh->GetEdgeVertices(face, verts);
      for (int k = 0; k < verts.Size(); k++)
      {
         int index = verts[k];
         // TODO: Get rid of all get_indicator() funcalls, replace with use of class int problem.
         if (pb->get_indicator() == "saltzmann")
         {
            // Replace the bdr attribute in the array as long as it is not
            // the dirichlet condition (For Saltzmann Problem)
            // This ensures the left wall vertices have the proper indicator
            if (BdrVertexIndexingArray[index] != 1)
            {
               BdrVertexIndexingArray[index] = bdr_attr;
            }
         }
         else if (pb->get_indicator() == "TriplePoint" || 
                  pb->get_indicator() == "Sedov" || 
                  pb->get_indicator() == "SodRadial")
         {
            // Mark corner vertices as 5
            // These nodes should not move at all during the simulation
            // Identify these corner vertices as those that already have
            // a value non negative
            if (BdrVertexIndexingArray[index] != -1 && BdrVertexIndexingArray[index] != bdr_attr)
            {
               cout << "vertex " << index << " is a corner.\n";
               BdrVertexIndexingArray[index] = 5;
            }
            else 
            {
               BdrVertexIndexingArray[index] = bdr_attr;
            }
         }
         else
         {
            BdrVertexIndexingArray[index] = bdr_attr;
         }
      } // end vertex iterator

   } // end boundary elements
}


/****************************************************************************************************
* Function: FillCellBdrFlag
*
* Purpose: 
*  To fill an L2 gridfunction of size NDofs_L2 (num cells) with a value indicating if the cell
*  is on the boundary or not.  This ParGridFunction can be retrieved with GetCellBdrFlagGF().
*
*  Note that this function will label a cell with 5 if that cell lies on the corner, indicating 
*  that the cell velocity should be set to 0 to preserve the slip BCs on both of its faces.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::FillCellBdrFlag()
{
   for (int i = 0; i < pmesh->GetNBE(); i++)
   {
      int bdr_attr = pmesh->GetBdrAttribute(i);
      int face = pmesh->GetBdrFace(i);

      Array<int> row;
      // Get cell
      face_element->GetRow(face, row);
      
      // Each face should only have 1 adjacent cell
      assert(row.Size() == 1);

      if (cell_bdr_flag_gf[row[0]] != -1 && 
          (pb->get_indicator() == "TriplePoint" || 
           pb->get_indicator() == "Sedov" || 
           pb->get_indicator() == "SodRadial"))
      {
         // Corner cell
         cout << "we have a corner cell: " << row[0] << endl;
         cell_bdr_flag_gf[row[0]] = 5;
      }
      else
      {
         cell_bdr_flag_gf[row[0]] = bdr_attr;
      }
      
      // cout << "cell: " << row[0] << ", bdr_attr: " << bdr_attr << endl;
   }
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
   //      << "             MakeTimeStep               \n"
   //      << "========================================\n";
   ComputeMeshVelocities(S, t, dt);

   // Update state variables contained in S_new
   // chrono_state.Start();
   ComputeStateUpdate(S, t, dt);
   // chrono_state.Stop();

   // Move the mesh
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction x_gf, mv_gf;
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);
   mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);

   // cout << "Printing select nodal velocities.\n";
   // Vector tmp_vel(dim);
   // Array<int> des_nodes({49,50,51,52,53,54,55});
   // for (int i = 0; i < des_nodes.Size(); i++)
   // {
   //    cout << "nodal velocity at " << des_nodes[i] << " is: ";
   //    GetNodeVelocity(S, des_nodes[i], tmp_vel);
   //    tmp_vel.Print(cout);
   // }

   add(x_gf, dt, mv_gf, x_gf);
   pmesh->NewNodes(x_gf, false);
   // cout << "mm computation took " << chrono_mm.RealTime() << "s.\n";
   // cout << "state update computation took " << chrono_state.RealTime() << "s.\n";
} 


/****************************************************************************************************
* Function: ComputeMeshVelocities
* Parameters:
*     S_new - BlockVector that stores mesh information, mesh velocity, and state variables.
*     t     - Current time
*     dt    - Current timestep
*
* Purpose:
*     Determine the velocity with which to move the mesh.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeMeshVelocities(Vector &S, const double &t, double &dt)
{
   // cout << "========================================\n"
   //      << "         ComputeMeshVelocities          \n"
   //      << "========================================\n";
   // chrono_mm.Start();
   /* Set mesh velocity for the previous iteration before we make 
   any modifications to the mesh velocity object contained in S */
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction mv_gf;
   mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);
   mv_gf_prev_it = mv_gf;

   if (mm)
   {
      ComputeIntermediateFaceVelocities(S, t);
      
      if (dim == 1)
      {
         /* Use intermediate face velocities to move the mesh */
         for (int face = 0; face < num_faces; face++)
         {
            Vector _vel(dim);
            GetIntermediateFaceVelocity(face, _vel);
            UpdateNodeVelocity(S, face, _vel);
         }
      }
      else 
      {
         switch (mv_option)
         {
         case 1:
            ComputeGeoVRaviart(S);
            break;

         case 2:
            ComputeGeoVNormal(S);
            break;
         
         case 3:
            ComputeGeoVCellFaceNormal(S);
            break;

         case 4:
            ComputeGeoVCAVEAT(S);
            break;
         
         case 5:
            ComputeGeoVCAVEATCellFace(S);
            break;
         
         case 6:
            ComputeGeoVCAVEATCellFaceWeighted(S);
            break;
         
         default:
            MFEM_ABORT("Invalid mesh velocity option.\n");
         }

         ComputeNodeVelocitiesFromVgeo(S, t, dt);

         if (this->use_corner_velocity_MC_iteration)
         {
            double val = ComputeIterationNorm(S,dt);
            cout << "Starting val: " << val << endl;
            for (int i = 1; i < corner_velocity_MC_num_iterations+1; i++)
            {
               // cout << "iterating on the corner node velocities\n";
               IterativeCornerVelocityMC(S, dt);
               double val = ComputeIterationNorm(S,dt);
               mv_gf_prev_it = mv_gf;
               cout << "val at iteration " << i << ": " << val
                    << ", mv_prev_norm: " << mv_gf_prev_it.Norml2() <<  endl;
            }
         }

         EnforceMVBoundaryConditions(S,t,dt);

         switch (fv_option)
         {
         case 1:
            ComputeCorrectiveFaceVelocities(S, t, dt);
            break;
         
         case 2:
            FillFaceVelocitiesWithAvg(S);
            break;
         
         default:
            /* Do nothing */
            break;
         } // End face velocity switch case
      }
      
      FillCenterVelocitiesWithAvg(S);
   }
   // chrono_mm.Stop();
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
   //      << "          ComputeStateUpdate            \n"
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
      Array<int> cell_bdr_arr;
      sum_validation = 0.;
      GetCellStateVector(S, ci, U_i);
      val = U_i;
      // cout << "cell state vector at tn:\n";
      // val.Print(cout);

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
            break;
         }
         default:
         {
            MFEM_ABORT("Incorrect dimension provided.\n");
         }
      }

      const DenseMatrix F_i = pb->flux(U_i, pmesh->GetAttribute(ci));
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
            DenseMatrix dm = pb->flux(U_j, pmesh->GetAttribute(cj));
            // cout << "flux at cj = " << cj << ":\n";
            // dm.Print(cout);
            dm += F_i; 
            // cout << "flux at ci = " << ci << ":\n";
            // F_i.Print(cout);
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

               /* REMOVE, checking visocity */
               // Array<int> tmp_dofs(2);
               // tmp_dofs[0] = 1, tmp_dofs[1]=2;
               // Vector tmp_vel(dim);
               // z.GetSubVector(tmp_dofs, tmp_vel);
               // if (tmp_vel[0] < 0.)
               // {
               //    cout << "\tcells ci " << ci << " and cj " << cj << " have a neg vel contribution:\n";
               //    tmp_vel.Print(cout);
               //    cout << "corresponding d: " << d << endl;
               // }
            }
         }
         else
         {
            assert(FI.IsBoundary());
            is_boundary_cell = true;
            int bdr_attr = BdrElementIndexingArray[fids[j]];
            // cout << "boundary attribute for face " << fids[j] << ": " << bdr_attr << endl;
            cell_bdr_arr.Append(bdr_attr);

            Vector y_temp(dim+2), y_temp_bdry(dim+2), U_i_bdry(dim+2);
            F_i.Mult(c, y_temp);
            U_i_bdry = U_i;

            /* Enforce Boundary Conditions */
            // if (pb->get_indicator() == "Sod" && (bdr_attr == 1 || bdr_attr == 3))
            // {
            //    /** 
            //     * numerical flux on boundary faces should be 
            //     *            |  0^T |
            //     *     f(U) = | pI_D |
            //     *            |  0^T |
            //    */
            //    y_temp[0] = 0.;
            //    y_temp[dim+1] = 0.;
            //    y_temp *= 2.;
            // }
            if (pb->get_indicator() == "saltzmann")
            {
               // Check bdry flag
               int bdr_attribute = BdrElementIndexingArray[fids[j]]; 

               switch (bdr_attribute)
               {
               case 1:
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
                  break;
               }
               // case 2:
               // case 3:
               // case 4:
               // {
               //    // Negate velocity
               //    for (int _it = 0; _it < dim; _it++)
               //    {
               //       U_i_bdry[_it + 1] = U_i_bdry[_it + 1] * -1;
               //    }
               //    DenseMatrix F_i_slip = pb->flux(U_i_bdry, pmesh->GetAttribute(i));
               //    F_i_slip.Mult(c, y_temp_bdry);
               //    // y_temp *= 2.;
               //    break;
               // }
               
               default:
               {
                  // cout << "invalid boundary attribute: " << bdr_attribute << endl;
                  // cout << "cell : " << ci << endl;
                  // cout << "face: " << fids[j] << endl;  
                  y_temp *= 2.; 
                  break;
               }
               } // switch (bdr_attr)
               y_temp += y_temp_bdry;
            } // if saltzmann
            else
            {
               y_temp *= 2.;
            }
            
            // Add in boundary contribution
            sums -= y_temp;
         } // End boundary face        

      } // End Face iterator

      // Enforce exact condition on boundary for Isentropic Vortex
      if (pb->get_indicator() == "IsentropicVortex" && is_boundary_cell)
      {
         EnforceExactBCOnCell(S, ci, t, dt, val);
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
      
      if (pb->has_boundary_conditions())
      {
         /* Iterate over cell boundaries, i.e. corner cells could have multiple bc flags */
         // Post processing modify computed values to enforce BCs
         if (is_boundary_cell)
         {
            for (int cell_bdr_it=0; cell_bdr_it < cell_bdr_arr.Size(); cell_bdr_it++)
            {
               int cell_bdr = cell_bdr_arr[cell_bdr_it];
               Array<int> tmp_dofs(2);
               tmp_dofs[0] = 1, tmp_dofs[1]=2;
               Vector tmp_vel(dim), normal(dim);
               normal = 0.;
               val.GetSubVector(tmp_dofs, tmp_vel);

               if (pb->get_indicator() == "Sod" && dim > 1)
               {
                  if (cell_bdr == 1)
                  {
                     // bottom
                     normal[1] = -1.;
                  }
                  else if (cell_bdr == 3)
                  {
                     // top
                     normal[1] = 1.;
                  }

                  double coeff = tmp_vel * normal;
                  normal *= coeff;
                  subtract(tmp_vel, normal, tmp_vel);
                  val.SetSubVector(tmp_dofs, tmp_vel);
               }
               else if (pb->get_indicator() == "saltzmann")
               {
                  if (cell_bdr == 2)
                  {
                     // bottom
                     normal[1] = -1.;
                  }
                  else if (cell_bdr == 3)
                  {
                     // right
                     normal[0] = 1.;
                  }
                  else if (cell_bdr == 4)
                  {
                     // top
                     normal[1] = 1.;
                  }

                  double coeff = tmp_vel * normal;
                  normal *= coeff;
                  subtract(tmp_vel, normal, tmp_vel);

                  if (cell_bdr == 1)
                  {
                     tmp_vel = 0.;
                     /* Ramping up to ex */
                     if (timestep_first == 0.)
                     {
                        timestep_first = timestep;
                     }
                     double _xi = t / (2*timestep_first);
                     double _psi = (4 - (_xi + 1) * (_xi - 2) * ((_xi - 2) - (abs(_xi-2) + (_xi-2)) / 2)) / 4.;
                     tmp_vel[0] = 1. * _psi;
                  }

                  val.SetSubVector(tmp_dofs, tmp_vel);


               }
               else if (pb->get_indicator() == "TriplePoint" || 
                        pb->get_indicator() == "SodRadial" ||
                        pb->get_indicator() == "Sedov")
               {
                  switch (cell_bdr)
                  {
                  case 0:
                     // not a boundary cell
                     continue;
                  case 1:
                     // bottom
                     // cout << "indic: " << 1 << endl;
                     normal[1] = -1.;
                     break;
                  
                  case 2:
                     // right
                     // cout << "indic: " << 2 << endl;
                     normal[0] = 1.;
                     break;
                  
                  case 3:
                     // top
                     // cout << "indic: " << 3 << endl;
                     normal[1] = 1.;
                     break;
                  
                  case 4:
                     // left
                     // cout << "indic: " << 4 << endl;
                     normal[0] = -1;
                     break;
                  
                  default:
                     MFEM_ABORT("Not a valid cell_bdr index.\n");
                  }

                  // cout << "normal: ";
                  // normal.Print(cout);
                  // cout << "old vel: ";
                  // tmp_vel.Print(cout);
                  double coeff = tmp_vel * normal;
                  normal *= coeff;
                  subtract(tmp_vel, normal, tmp_vel);
                  // cout << "new vel: ";
                  // tmp_vel.Print(cout);
                  
                  coeff = tmp_vel * normal;
                  if (coeff > 1E-12)
                  {
                     cout << "coeff: " << coeff << endl;
                     cout << "Normal component not removed: " << coeff << endl;
                     tmp_vel.Print(cout);
                     assert(false);
                  }      
                  val.SetSubVector(tmp_dofs, tmp_vel);
               }
               else if (pb->get_indicator() == "TestBCs" ||
                        pb->get_indicator() == "Vdw1" ||
                        pb->get_indicator() == "Vdw2" ||
                        pb->get_indicator() == "Vdw3" ||
                        pb->get_indicator() == "Vdw4")
               {
                  switch (cell_bdr)
                  {
                  case 0:
                     // not a boundary cell
                     continue;
                  case 1:
                     // bottom
                     // cout << "indic: " << 1 << endl;
                     normal[1] = -1.;
                     break;
                  
                  case 2:
                     // right
                     // cout << "indic: " << 2 << endl;
                     normal[0] = 1.;
                     break;
                  
                  case 3:
                     // top
                     // cout << "indic: " << 3 << endl;
                     normal[1] = 1.;
                     break;
                  
                  case 4:
                     // left
                     // cout << "indic: " << 4 << endl;
                     normal[0] = -1;
                     break;
                  
                  default:
                     MFEM_ABORT("Not a valid cell_bdr index.\n");
                  }

                  double coeff = tmp_vel * normal;
                  normal *= coeff;
                  subtract(tmp_vel, normal, tmp_vel);
                  val.SetSubVector(tmp_dofs, tmp_vel);
               }
            } // End post processing BC application
         
         } // End is_bdry_cell
      } // End pb->has_boundary_conditions()

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
   //      << "         EnforceExactBCOnCell           \n"
   //      << "========================================\n";

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
*
*  Note that this function assumes dim > 1.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::EnforceMVBoundaryConditions(Vector &S, const double &t, const double &dt)
{
   // cout << "========================================\n"
   //      << "     EnforcingMVBoundaryConditions      \n"
   //      << "========================================\n";
   assert(dim > 1);

   if (pb->get_indicator() == "Sod")
   {
      int bdr_ind = 0;
      Vector normal(dim), node_v(dim);
      for (int i = 0; i < NVDofs_H1; i++)
      {
         bdr_ind = BdrVertexIndexingArray[i];
         /* Get corresponding normal vector */
         switch (bdr_ind)
         {
         case 1:
            // bottom
            normal[0] = 0., normal[1] = -1.;
            break;
         
         case 2:
            // right
            continue;
         
         case 3:
            // top
            normal[0] = 0., normal[1] = 1.;
            break;
         
         case 4:
            // left
            continue;
         
         case -1:
            // Not a boundary vertex
            continue;

         default:
            MFEM_ABORT("Incorrect bdr attribute encountered while enforcing mesh velocity BCs.\n");
            break;
         }

         if (bdr_ind == -1) 
         {
            MFEM_ABORT("Dont correct interior vertices.\n");
         }
         /* Correct node velocity accordingly */
         GetNodeVelocity(S, i, node_v);

         double coeff = node_v * normal;
         node_v.Add(-coeff, normal);

         UpdateNodeVelocity(S, i, node_v);
      }
   }
   else if (pb->get_indicator() == "TriplePoint" ||
            pb->get_indicator() == "Sedov" ||
            pb->get_indicator() == "SodRadial")
   {
      int bdr_ind = 0;
      Vector normal(dim), node_v(dim);
      for (int i = 0; i < NVDofs_H1; i++)
      {
         // cout << "enforcing MV BCs on vertex: " << i << endl;
         bdr_ind = BdrVertexIndexingArray[i];
         /* Get corresponding normal vector */
         switch (bdr_ind)
         {
         case 1:
            // bottom
            normal[0] = 0., normal[1] = -1.;
            /* code */
            break;
         
         case 2:
            // right
            normal[0] = 1., normal[1] = 0.;
            break;
            
         case 3:
            // top
            normal[0] = 0., normal[1] = 1.;
            break;
         
         case 4:
            // left
            normal[0] = -1., normal[1] = 0.;
            break;

         case 5:
            node_v = 0.;
            UpdateNodeVelocity(S, i, node_v);
            continue;
         
         case -1:
            // Not a boundary vertex
            continue;

         default:
            MFEM_ABORT("Incorrect bdr attribute encountered while enforcing mesh velocity BCs.\n");
            break;
         }

         if (bdr_ind == -1) 
         {
            MFEM_ABORT("Dont correct interior vertices.\n");
         }
         /* Correct node velocity accordingly */
         GetNodeVelocity(S, i, node_v);
         double coeff = node_v * normal;
         node_v.Add(-coeff, normal);
         UpdateNodeVelocity(S, i, node_v);
      }
   }
   else if (pb->get_indicator() == "saltzmann")
   {
      // Enforce
      // Vector* sptr = const_cast<Vector*>(&S);
      // ParGridFunction x_gf, mv_gf;
      // mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);

      Vector ex(dim);
      ex = 0.;
      /* Ramping up to ex */
      if (timestep_first == 0.)
      {
         timestep_first = timestep;
      }
      double _xi = t / (2*timestep_first);
      double _psi = (4 - (_xi + 1) * (_xi - 2) * ((_xi - 2) - (abs(_xi-2) + (_xi-2)) / 2)) / 4.;
      ex[0] = 1. * _psi;
      // ex[0] = 1.;

      // VectorConstantCoefficient left_wall_coeff(ex);
      // Array<int> left_wall_attr(1);
      // left_wall_attr[0] = 1;

      // mv_gf.ProjectBdrCoefficient(left_wall_coeff, left_wall_attr);

      /* Enforce v.n=0 on rest of boundary */
      int bdr_ind = 0;
      Vector normal(dim), node_v(dim);
      for (int i = 0; i < NVDofs_H1; i++)
      {
         bdr_ind = BdrVertexIndexingArray[i];

         /* Get corresponding normal vector */
         switch (bdr_ind)
         {
         case 1:
            // left
            UpdateNodeVelocity(S,i,ex);
            continue;
         
         case 2:
            // bottom
            normal[0] = 0., normal[1] = -1.;
            break;
         
         case 3:
            // right
            normal[0] = 1., normal[1] = 0.;
            break;
         
         case 4:
            // top
            normal[0] = 0., normal[1] = 1.;
            break;
         
         case -1:
            // Not a boundary vertex
            continue;

         default:
            MFEM_ABORT("Incorrect bdr attribute encountered while enforcing mesh velocity BCs.\n");
            break;
         }

         if (bdr_ind == -1) 
         {
            MFEM_ABORT("Dont correct interior vertices.\n");
         }
         else if (bdr_ind == 1)
         {
            MFEM_ABORT("Left wall BCs already handled.\n");
         }
         /* Correct node velocity accordingly */
         GetNodeVelocity(S, i, node_v);

         double coeff = node_v * normal;

         node_v.Add(-coeff, normal);
         UpdateNodeVelocity(S, i, node_v);
      }
   }
   else if (pb->get_indicator() == "IsentropicVortex")
   {
      // Boundary vertices should not move at all
      Vector* sptr = const_cast<Vector*>(&S);
      ParGridFunction x_gf, mv_gf;
      mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);

      Vector zero(dim);
      zero = 0.;
      VectorConstantCoefficient zero_coeff(zero);

      mv_gf.ProjectBdrCoefficient(zero_coeff, pmesh->bdr_attributes);
   }
   // Remove normal velocity at ALL boundaries
   else if (pb->get_indicator() == "TestBCs" ||
            pb->get_indicator() == "Vdw1" ||
            pb->get_indicator() == "Vdw2" ||
            pb->get_indicator() == "Vdw3" ||
            pb->get_indicator() == "Vdw4")
   {
      int bdr_ind = 0;
      Vector normal(dim), node_v(dim);
      for (int i = 0; i < NVDofs_H1; i++)
      {
         bdr_ind = BdrVertexIndexingArray[i];
         /* Get corresponding normal vector */
         switch (bdr_ind)
         {
         case 1:
            // bottom
            normal[0] = 0., normal[1] = -1.;
            break;
         
         case 2:
            // right
            normal[0] = 1., normal[1] = 0.;
            break;
         
         case 3:
            // top
            normal[0] = 0., normal[1] = 1.;
            break;
         
         case 4:
            // left
            normal[0] = -1., normal[1] = 0.;
            break;
         
         case -1:
            // Not a boundary vertex
            continue;

         default:
            MFEM_ABORT("Incorrect bdr attribute encountered while enforcing mesh velocity BCs.\n");
            break;
         }

         if (bdr_ind == -1) 
         {
            MFEM_ABORT("Dont correct interior vertices.\n");
         }
         /* Correct node velocity accordingly */
         GetNodeVelocity(S, i, node_v);

         double coeff = node_v * normal;
         node_v.Add(-coeff, normal);

         UpdateNodeVelocity(S, i, node_v);
      }
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
*
*  The definition of the intermediate face velocity is given by eq (5.7).
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::
   ComputeIntermediateFaceVelocities(const Vector &S,
                                        const double t,
                                        const string flag, // Default NA
                                        void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
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
            pb->velocity(Uc, Vf);
            Vector vcp; pb->velocity(Ucp, vcp);
            Vf.Add(1., vcp);
            Vf *= 0.5;
            if (use_viscosity)
            {
               double pcp = pb->pressure(Ucp,pmesh->GetAttribute(cp));
               double pc = pcp - pb->pressure(Uc,pmesh->GetAttribute(c));
               // double coeff = d * pc * (Ucp[0] - Uc[0]) / F; // This fixes the upward movement in tp
               double coeff = d * (Ucp[0] - Uc[0]) / F; // This is how 5.7b is defined.
               Vf.Add(coeff, n_vec);

               // if (pmesh->GetAttribute(cp) != pmesh->GetAttribute(c))
               // {
               //    cout << "cells " << c << " and " << cp << " have different attributes.\n";
               //    cout << "pressure differential: " << pc << endl;
               //    cout << "coeff: " << coeff << endl;
               // }
            }
            
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
            pb->velocity(Uc, Vf);             
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
*     3) Adjacent cell face based mesh movement (LS based)
*     4) CAVEAT Weighted LS
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::SetMVOption(const int & option)
{
   this->mv_option = option;
}


/****************************************************************************************************
* Function: SetFVOption
* Parameters:
*  option - integer parameter indicating which face velocity to use.
*
* Purpose:
*  To set the mesh velocity motion computation type:
*     1) Mass corrective.
*     2) Average of corners, Q1 type.
*     Default) Use computed velocity, RT or Vf
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::SetFVOption(const int & option)
{
   this->fv_option = option;
}


/****************************************************************************************************
* Function: SetMVIteration
* Parameters:
*  num_iterations - integer parameter indicating the number of times to iterate
*                   on the mesh velocities
*
* Purpose:
*  To enable the use of and set the number of iterations for the 
*  IterativeCornerVelocityMC() function.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::SetMVIteration(const int num_iterations) { 
   this->use_corner_velocity_MC_iteration = true;
   this->corner_velocity_MC_num_iterations = num_iterations; 
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
      if (pb->get_indicator() == "Noh" && cell_bdr_flag_gf[ci] != -1)
      {
         // Skip boundary cells
         continue;
      }
      const double m = m_hpv->Elem(ci);
      const double k = pmesh->GetElementVolume(ci);
      GetCellStateVector(S, ci, U_i);
      num += abs((k / U_i[0]) - m);
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

   cout << "Percentage of cells where mass conservation was broken: " << cell_ratio << endl
        << "Initial mass sum: " << denom 
        << ", Current mass sum: " << current_mass_sum << endl
        << "Mass Error: " << num / denom << endl
        << "--------------------------------------\n";
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
   // cout << "==============================================\n";
   // cout << "Computing corrective interior face velocities.\n";
   // cout << "==============================================\n";

   assert(dim > 1); // No need to correct face velocities in dim=1

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
         pb->velocity(Uc, cell_c_v);
         pb->velocity(Ucp, cell_cp_v);

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
         // Vector rt_face_vel(dim);
         // GetNodeVelocity(S, face_dof, rt_face_vel);
         // V3nperp = rt_face_vel * n_vec_perp;

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
   // cout << "========================================================\n";
   // cout << "Computing corrective interior geometric face velocities.\n";
   // cout << "========================================================\n";
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

         // Compute C_i
         DenseMatrix Ci;
         // mfem::Array<int> elements;
         // face_element->GetRow(face_dof, elements);
         ComputeCiGeo(face_dof, Ci);

         // Compute alpha_i
         double alpha_i;
         ComputeDeterminant(Ci, dt, alpha_i);

         // Compute flux
         Vector face_flux(dim);
         DenseMatrix _mat(Ci);
         _mat *= - dt / 2.;
         for (int i = 0; i < dim; i++)
         {
            _mat(i,i) += alpha_i;
         }

         _mat.Mult(face_velocity, face_flux);
         SetCorrectedFaceFlux(face, face_flux); 
      }
   }
}

/****************************************************************************************************
* Function: SetCellCenterAsCenter
* Parameters:
*   S        - BlockVector representing FiniteElement information
*
* Purpose:
*  This function is only needed since we are not able to implement the Serendipity finite elements.
*  Since the mesh motion does not depends on a value for the node corresponding to the center of the 
*  element, we fill the actual grid function with the average of the moved adjacent corners
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::SetCellCenterAsCenter(Vector &S)
{
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction x_gf;
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);

   // Since we cannot use Serendipity elements, we must update cell center velocities
   Vector cell_x(dim), node_x(dim);

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

      cell_x = 0.;
      pmesh->GetElementVertices(ci, verts);

      for (int j = 0; j < verts.Size(); j++)
      {
         GetNodePosition(S, verts[j], node_x);
         cell_x += node_x;
      }

      cell_x /= verts.Size();

      UpdateNodePosition(S, cell_vdof, cell_x);
   }
}


/****************************************************************************************************
* Function: FillCenterVelocitiesWithL2
* Parameters:
*   S        - BlockVector representing FiniteElement information
*
* Purpose:
*  This function is only needed since we are not able to implement 
*  the Serendipity finite elements. Since the mesh motion does not 
*  depends on a value for the node corresponding to the center of 
*  the element, we need to fill this with a velocity which will 
*  ensure the center node remains inside the cell.  This is done 
*  by taking the hydrodynamic velocity at the cell.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::FillCenterVelocitiesWithL2(Vector &S)
{
   // Since we cannot use Serendipity elements, we must update cell center velocities
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction v_gf;
   v_gf.MakeRef(&H1, *sptr, block_offsets[3]);

   Vector Uc(dim+2), node_v(dim);

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

      GetCellStateVector(S, ci, Uc);
      pb->velocity(Uc, node_v);
      
      // Get corresponding velocity from ParGridFunction
      UpdateNodeVelocity(S, cell_vdof, node_v);
   }
}


/****************************************************************************************************
* Function: FillCenterVelocitiesWithAvg
* Parameters:
*   S        - BlockVector representing FiniteElement information
*
* Purpose:
*  This function is only needed since we are not able to implement 
*  the Serendipity finite elements. Since the mesh motion does not 
*  depends on a value for the node corresponding to the center of 
*  the element, we need to fill this with a velocity which will 
*  ensure the center node remains inside the cell.  This is done 
*  by average the values at the four corner nodes of the cell.
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
         pmesh->GetElementVertices(ci, verts);
         for (int j = 0; j < verts.Size(); j++)
         {
            GetNodeVelocity(S, verts[j], node_v);
            Vc += node_v;
         }

         Vc /= verts.Size();
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
   Vector Uc(dim+2), Ucp(dim+2), face_velocity(dim);
   Vector cell_c_v(dim), cell_cp_v(dim); // For computation of bmn
   Array<int> row;
   int face_dof = 0, face_vdof1 = 0, face_vdof2 = 0;

   // Iterate over faces
   for (int face = 0; face < num_faces; face++) // face iterator
   {  
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

      // Lastly, put face velocity into gridfunction object
      UpdateNodeVelocity(S, face_dof, face_velocity);
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
* Function: SetViGeo
* Parameters:
*   node - Global index of node in question.
*   vel  - Velocity Vi^geo at given node
*
* Purpose:
*  This function sets ViGeo.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::SetViGeo(const int &node, const Vector &vel)
{
   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NVDofs_H1;
      v_geo_gf[index] = vel[i];
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
      int index = node + i * NVDofs_H1;
      vel[i] = v_geo_gf[index];
   }
}


/****************************************************************************************************
* Function: UpdateNodePosition
* Parameters:
*   S    - BlockVector representing FiniteElement information
*   node - Global index of node in question.
*   x    - Coordinate position at given node.
*
* Purpose:
*  This function updates the cartesian location corresponding to a global node.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::UpdateNodePosition(Vector &S, const int & node, const Vector &x)
{
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction x_gf;
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);

   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      x_gf[index] = x[i];
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
   Vector center(dim), U(dim+2);
   double pressure=0., ss=0., x_val=0.;

   // Form filenames and ofstream objects
   std::string sv_file = output_file_prefix + output_file_suffix;
   std::ofstream fstream_sv(sv_file.c_str());
   fstream_sv << "x,rho,v,ste,p,ss,cell_type\n";

   for (int i = 0; i < NDofs_L2; i++)
   {
      // compute pressure and sound speed on the fly
      GetCellStateVector(S, i, U);
      pressure = pb->pressure(U, pmesh->GetAttribute(i));
      ss = pb->sound_speed(U, pmesh->GetAttribute(i));
      
      pmesh->GetElementCenter(i, center);

      // We fill the x-coordinate depending on the problem.
      // This really depends on what kind of visualization is needed
      switch (problem)
      {
      // Any radial plot, we take the L2 norm of the cell center
      case 4:  // Noh
      case 6:  // Sedov
      case 13: // Sod Radial
         x_val = center.Norml2();
         break;
      
      // For stacked Sod problem, we only need the x-coord
      case 0:  // Smooth
      case 1:  // Sod
      case 2:  // Lax
      case 3:  // Leblanc
      case 8:  // Vdw1
      case 9:  // Vdw2
      case 10: // Vdw3
      case 11: // Vdw4
      case 20:
      default:
         x_val = center[0];
         break;
      }

      // if (abs(x_val - 0.6) < 0.02)
      // {
      //    cout << "=============\n"
      //         << "cell: " << i << endl 
      //         << "attr: " << BdrElementIndexingArray[i] << endl
      //         << "x_val: " << x_val << endl
      //         << "coords: ";
      //    center.Print(cout);
      //    cout << "rho: " << 1./U[0] << endl << endl;
      // }
      
      fstream_sv << x_val << ","    // x
               << 1./U[0] << ","  // rho
               << U[1] << ","     // v_x
               << U[dim+1] << "," // ste
               << pressure << "," // pressure
               << ss << ",";     // sound speed

      // Print flag if interior or bdr
      if (cell_bdr_flag_gf[i] == -1.)
      {
         fstream_sv << "int\n";
      }
      else 
      {
         fstream_sv << "bdr\n";
      }
   }
}


/****************************************************************************************************
* Function: Normal direction based Velocity Functions
****************************************************************************************************/

/****************************************************************************************************
* Function: tensor
* Parameters:
*  v1     - Vector object
*  v2     - Vector object
*  dm     - Resultant DenseMatrix
*
* Purpose:
*  Compute the tensor product of two vectors.
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
* Function: ComputeGeoVNormal
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function computes the geometric velocity from the adjacent face normals
*  and face velocities as follows:
*
*  Define Mi    = Sum (nf x nf)
*         Ri    = Sum [(nf x nf) Vf]
*         Vigeo = Mi^-1 Ri
*
*  Note: Vigeo is only the geometric velocity, and is not the velocity used 
*        to move the mesh.  For the mesh velocity, one must solve for the 
*        alpha_i, which process is found in ComputeNodeVelocitiesFromVgeo().
*
*  Note: This function fills the Q1 velocity field v_geo_gf and also the mv_gf
*        in the BlockVector S.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeGeoVNormal(Vector &S)
{
   mfem::Mesh::FaceInformation FI;
   Vector node_v(dim), Ri(dim), n_vec(dim), Vf(dim), y(dim);
   DenseMatrix Mi(dim), dm_tmp(dim);

   Array<int> faces_row, oris;
   int faces_length;
   double c_norm;

   DofEntity entity;
   int EDof;

   // Iterate over all mesh velocity dofs
   for (int node = 0; node < NDofs_H1; node++) // Vertex iterator
   {
      node_v = 0.; // Reset node velocity
      GetEntityDof(node, entity, EDof);

      switch (entity)
      {
         case 0: // corner
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

            // Only update the v_geo_gf in the case of the corners
            // Using Q1 for this velocity field
            SetViGeo(node, node_v);

            break;
         }
         case 1: // face
         {
            // face nodes just set to ivf
            assert(face >= 0);
            GetIntermediateFaceVelocity(EDof, node_v);

            break;
         }
         case 2: // Cell Center
         {
            Vector Uc(dim+2);
            GetCellStateVector(S, EDof, Uc);
            pb->velocity(Uc, node_v);

            break;
         }
         default:
         {
            MFEM_ABORT("Invalid entity\n");
         }
      }
      // cout << "node velocity from mv2 for node (" << node << "): ";
      // node_v.Print(cout);
      // In every case, update the corresponding nodal velocity in S
      UpdateNodeVelocity(S, node, node_v);

   } // End node iterator
}


/****************************************************************************************************
* Function: ComputeGeoVCellFaceNormal
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function computes the geometric velocity from the face normals and face 
*  velocities from all adjacent cells as follows:
*
*  CellFaceNormalTODO:
*        | num cells | num faces |
*        |-----------|-----------|
*        |     1     |     4     |
*        |     2     |     7     |
*        |     4     |    12     |
*
*  Note: Vigeo is only the geometric velocity, and is not the velocity used 
*        to move the mesh.  For the mesh velocity, one must solve for the 
*        alpha_i, which process is found in ComputeNodeVelocitiesFromVgeo().
*
*  Note: This function fills the Q1 velocity field v_geo_gf and also the mv_gf
*        in the BlockVector S.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeGeoVCellFaceNormal(Vector &S)
{
   cout << "ComputeGeoVCellFaceNormal call\n";
   mfem::Mesh::FaceInformation FI;
   
   Vector node_v(dim), node_x(dim), n_vec(dim), Vf(dim);
   Vector vdof1_x(dim), vdof2_x(dim), face_x(dim), vdof_avg(dim);
   double F;

   Array<int> element_row, cell_faces_row, oris, face_dofs;
   int row_index, length_element_row, length_cell_faces, num_rows_A, num_cols_A = 6;

   DenseMatrix A, AT, ATA(num_cols_A);
   Vector RHS, RHS_mod(num_cols_A), res(num_cols_A), temp_vec;

   DofEntity entity;
   int EDof;

   std::unordered_set<int> uniqueFaceIndices;

   // Iterate over all mesh velocity dofs
   // cout << "ndofs_h1: " << NDofs_H1 << endl;
   for (int node = 0; node < NDofs_H1; node++) // Vertex iterator
   {
      // cout << "NODE: " << node << endl;
      node_v = 0.; // Reset node velocity
      GetEntityDof(node, entity, EDof);
      // cout << "EDof: " << EDof << endl;

      switch (entity)
      {
         case 0: // corner
         {
            // Reset vals
            uniqueFaceIndices.clear();
            F = 0.;
            GetNodePosition(S, node, node_x);

            // Get adjacent cells
            vertex_element->GetRow(node, element_row);
            // cout << "element_row: ";
            // element_row.Print(cout);
            length_element_row = element_row.Size();

            // Set size of overdetermined matrix A depending on number of adjacent cells
            switch (length_element_row)
            {
               case 1:
               {
                  // TODO: What is to be done at the corners of our domain?
                  node_v = 1.;
                  continue;
               }
               case 2:
               {
                  // TODO: This results in A^TA being a singular matrix
                  num_rows_A = 7;
                  cout << "===== boundary node: " << node << "\n";
                  continue;
               }
               case 4:
               {
                  num_rows_A = 12;
                  break;
               }
               default:
               {
                  MFEM_ABORT("Invalid number of adjacent ceells.\n");
                  break;
               }
            }
            A.SetSize(num_rows_A, num_cols_A);
            RHS.SetSize(num_rows_A);
            temp_vec.SetSize(num_rows_A);

            // cout << "===== new corner node ===== [at: " << node << "]\n";
            /*** Fill matrix representing overdetermined system and RHS vector ***/
            // Iterate over adjacent cells
            for (int element_it = 0; element_it < length_element_row; element_it++)
            {
               int el_index = element_row[element_it];

               pmesh->GetElementEdges(el_index, cell_faces_row, oris);
               length_cell_faces = cell_faces_row.Size();

               // Iterate over cell faces
               for (int face_it = 0; face_it < length_cell_faces; face_it++)
               {
                  int face_index = cell_faces_row[face_it];

                  // Check if cell face is unique
                  if (uniqueFaceIndices.find(face_index) == uniqueFaceIndices.end())
                  {
                     // New face, fill corresponding row in A and RHS
                     uniqueFaceIndices.insert(face_index);
                     row_index = uniqueFaceIndices.size() - 1;
                     // cout << "row_index: " << row_index << endl;
                     
                     // Get face velocity, outward normal, and face position
                     // n_vec is the outward normal times |F|
                     GetIntermediateFaceVelocity(face_index, Vf);
                     CalcOutwardNormalInt(S, el_index, face_index, n_vec);

                     double indicator = Vf * n_vec;
                     // if (indicator < 0.)
                     // {
                     //    n_vec *= -1.;
                     // }

                     // Get average of face edge nodes
                     H1.GetFaceDofs(face_index, face_dofs);
                     GetNodePosition(S, face_dofs[1], vdof1_x);
                     GetNodePosition(S, face_dofs[0], vdof2_x);
                     face_x = vdof1_x;
                     face_x.Add(1., vdof2_x);
                     face_x *= 0.5;
                     
                     // Cout info
                     // cout << "Unique face corresponds to el: " << el_index << ", face: " << face_index << endl;
                     // cout << "Vf: ";
                     // Vf.Print(cout);
                     // cout << "nf: ";
                     // n_vec.Print(cout);
                     // cout << "face_x: ";
                     // face_x.Print(cout);
                     // GetNodePosition(S, face_dofs[2], face_x);

                     // Fill rows of A
                     A(row_index, 0) = n_vec(0);
                     A(row_index, 1) = n_vec(1);
                     A(row_index, 2) = n_vec(0)*face_x(0);
                     A(row_index, 3) = n_vec(1)*face_x(0);
                     A(row_index, 4) = n_vec(0)*face_x(1);
                     A(row_index, 5) = n_vec(1)*face_x(1);

                     // Fill rows of RHS
                     RHS(row_index) = Vf * n_vec;
                     // cout << endl;
                  } // End unique new face
               } // End face iterator
            } // End cell iterator

            // cout << "RHS: ";
            // RHS.Print(cout);

            // Solve for vgeo on this node
            AT = A;
            AT.Transpose();
            Mult(AT, A, ATA);
            AT.Mult(RHS, RHS_mod);
            // cout << "A: " << endl;
            // A.Print(cout);
            // cout << "ATA: " << endl;
            // ATA.Print(cout);
            ATA.Invert();
            ATA.Mult(RHS_mod,res);

            // Compute residual
            A.Mult(res, temp_vec);
            // cout << "temp_vec: ";
            // temp_vec.Print(cout);

            temp_vec -= RHS;
            // cout << "residual: " << temp_vec.Norml2() << endl;

            // cout << "coefficient vector: ";
            // res.Print(cout);

            node_v(0) = res[0] + res[2]*node_x[0] + res[4]*node_x[1];
            node_v(1) = res[1] + res[3]*node_x[0] + res[5]*node_x[1];

            // cout << "node_x: ";
            // node_x.Print(cout);
            // cout << "node_v: ";
            // node_v.Print(cout);

            // Only update the v_geo_gf in the case of the corners
            // Using Q1 for this velocity field
            SetViGeo(node, node_v);

            break;
         } // End corner
         case 1: // face
         {
            // face nodes just set to ivf
            assert(face >= 0);
            GetIntermediateFaceVelocity(EDof, node_v);

            break;
         }
         case 2: // Cell Center
         {
            Vector Uc(dim+2);
            GetCellStateVector(S, EDof, Uc);
            pb->velocity(Uc, node_v);

            break;
         }
         default:
         {
            MFEM_ABORT("Invalid entity\n");
         }
      }
      if (node_v.Norml2() > 0.0001)
      {
         cout << "node: " << node << endl       
              << "nodev: ";
         node_v.Print(cout);
      }

      // In every case, update the corresponding nodal velocity in S
      UpdateNodeVelocity(S, node, node_v);

   } // End node iterator
}


/****************************************************************************************************
* Function: ComputeGeoVCAVEAT
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function computes the geometric velocity from the face normals and face 
*  velocities from all adjacent faces following the Caveat weighted least squares
*  algorithm described in https://www.osti.gov/biblio/5366795.
*
*  The weights used for the weighted least squares method is simply the sum of 
*  the densities across any given face.
*
*
*  Note: Vigeo is only the geometric velocity, and is not the velocity used 
*        to move the mesh.  For the mesh velocity, one must solve for the 
*        alpha_i, which process is found in ComputeNodeVelocitiesFromVgeo().
*
*  Note: This function fills the Q1 velocity field v_geo_gf and also the mv_gf
*        in the BlockVector S.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeGeoVCAVEAT(Vector &S)
{
   // cout << "==========ComputeGeoVCaveat==========\n";
   mfem::Mesh::FaceInformation FI;
   Vector node_v(dim), Ri(dim), n_vec(dim), Vf(dim), Uc(dim+2), Ucp(dim+2);
   DenseMatrix A, AT, WA, ATWA(dim,dim), W;
   Vector RHS, RHS_mod;

   Array<int> faces_row, oris;
   int faces_length, c, cp;
   double c_norm, weight;

   DofEntity entity;
   int EDof;

   // Iterate over all mesh velocity dofs
   for (int node = 0; node < NDofs_H1; node++) // Vertex iterator
   {
      // cout << "\n\t $$$node: " << node << endl;
      node_v = 0.; // Reset node velocity
      GetEntityDof(node, entity, EDof);

      switch (entity)
      {
         case 0: // corner
         {
            node_v = 0.;
            // Get cell faces
            vertex_edge.GetRow(node, faces_row);
            faces_length = faces_row.Size();

            // Set sizes of nodal matrices
            A.SetSize(faces_length, dim);
            AT.SetSize(dim, faces_length);
            W.SetSize(faces_length);
            WA.SetSize(faces_length, dim);
            RHS.SetSize(faces_length);
            RHS_mod.SetSize(dim);
            A = 0., AT = 0., ATWA = 0., W = 0., WA = 0., RHS = 0., RHS_mod = 0.;

            // iterate over adjacent faces
            for (int face_it = 0; face_it < faces_length; face_it++) // Adjacent face iterator
            {
               int face = faces_row[face_it];
               GetIntermediateFaceVelocity(face, Vf);
               FI = pmesh->GetFaceInformation(face);
               c = FI.element[0].index;
               cp = FI.element[1].index;

               // Retrieve corresponding normal (orientation doesn't matter)
               CalcOutwardNormalInt(S, FI.element[0].index, face, n_vec);
               c_norm = n_vec.Norml2();
               n_vec /= c_norm;

               // Compute weight
               double rhoc = 0., rhocp = 0.;
               GetCellStateVector(S, c, Uc);
               
               rhoc = 1. / Uc[0];
               if (FI.IsInterior())
               {
                  GetCellStateVector(S, cp, Ucp);
                  rhocp = 1 / Ucp[0];
               }
               weight = rhoc + rhocp;

               // Add in contribution to system
               A.SetRow(face_it, n_vec);
               W(face_it, face_it) = weight;
               double val = n_vec * Vf;
               RHS[face_it] = weight * val;

               // Print stuff
               // cout << "face: " << face << ", vf: ";
               // Vf.Print(cout);
               // cout << "normal: ";
               // n_vec.Print(cout);
               // cout << "weight: " << weight << endl;
            }
            // Matrix setup
            AT = A;
            AT.Transpose();
            Mult(W, A, WA);
            // cout << "WA: ";
            // WA.Print(cout);
            Mult(AT,WA,ATWA);
            // cout << "ATWA:\n";
            // ATWA.Print(cout);
            AT.Mult(RHS, RHS_mod);

            // Solve for best fit
            ATWA.Invert();
            ATWA.Mult(RHS_mod, node_v);

             // Print stuff
            // cout << "A:\n";
            // A.Print(cout);
            // cout << "AT:\n";
            // AT.Print(cout);
            // cout << "W:\n";
            // W.Print(cout);
            // cout << "RHS: ";
            // RHS.Print(cout);
            // cout << "RHS_mod: ";
            // RHS_mod.Print(cout);
            // cout << "node_v: ";
            // node_v.Print(cout); 

            // Only update the v_geo_gf in the case of the corners
            // Using Q1 for this velocity field
            SetViGeo(node, node_v);

            break;
         }
         case 1: // face
         {
            // face nodes just set to ivf
            assert(face >= 0);
            GetIntermediateFaceVelocity(EDof, node_v);

            break;
         }
         case 2: // Cell Center
         {
            Vector Uc(dim+2);
            GetCellStateVector(S, EDof, Uc);
            pb->velocity(Uc, node_v);

            break;
         }
         default:
         {
            MFEM_ABORT("Invalid entity\n");
         }
      }

      // In every case, update the corresponding nodal velocity in S
      UpdateNodeVelocity(S, node, node_v);

   } // End node iterator
}


/****************************************************************************************************
* Function: ComputeGeoVCAVEATCellFace
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function computes the geometric velocity from the face normals and face 
*  velocities from all adjacent faces.
*
*  On the boundary vertices, the above CAVEAT algorithm is used, which is 
*  described in https://www.osti.gov/biblio/5366795.
*
*  On the interior vertices, the above adjacent Cell face based method is used.
*
*  Note: Vigeo is only the geometric velocity, and is not the velocity used 
*        to move the mesh.  For the mesh velocity, one must solve for the 
*        alpha_i, which process is found in ComputeNodeVelocitiesFromVgeo().
*
*  Note: This function fills the Q1 velocity field v_geo_gf and also the mv_gf
*        in the BlockVector S.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeGeoVCAVEATCellFace(Vector &S)
{
   // cout << "ComputeGeoVCAVEATCellFace call\n";
   mfem::Mesh::FaceInformation FI;
   
   Vector node_v(dim), node_x(dim), n_vec(dim), Vf(dim);
   Vector vdof1_x(dim), vdof2_x(dim), face_x(dim), vdof_avg(dim);
   double F, c_norm;

   Array<int> element_row, cell_faces_row, oris, face_dofs;
   int row_index, length_element_row, length_cell_faces, num_rows_A, num_cols_A = 6;
   int c, cp;
   Vector Uc(dim + 2), Ucp(dim+2);

   DenseMatrix A, AT, ATA(num_cols_A), W, WA, ATWA(dim);
   Vector RHS, RHS_mod(num_cols_A), res(num_cols_A), temp_vec;

   DofEntity entity;
   int EDof;

   std::unordered_set<int> uniqueFaceIndices;

   // Iterate over all mesh velocity dofs
   for (int node = 0; node < NDofs_H1; node++) // Vertex iterator
   {
      node_v = 0.; // Reset node velocity
      GetEntityDof(node, entity, EDof);

      switch (entity)
      {
         case 0: // corner
         {
            // Reset vals
            uniqueFaceIndices.clear();
            F = 0.;
            node_v = 0.;
            GetNodePosition(S, node, node_x);

            // Get adjacent cells
            vertex_element->GetRow(node, element_row);
            length_element_row = element_row.Size();

            // Set size of overdetermined matrix A depending on number of adjacent cells
            switch (length_element_row)
            {
               case 1: // Domain corner vertex
               case 2: // Boundary vertex
               {
                  // Since in the Cell Face based method we get a singular matrix for
                  // the domain corner and boundary vertices, we use the CAVEAT method
                  // of projecting the vertex velocity onto the adjacent cell faces
                  // and solve the system using weighted least squares.
                  if (length_element_row == 1) {
                     num_rows_A = 4;
                  } else {
                     num_rows_A = 7;
                  }

                  // Set sizes of nodal matrices
                  A.SetSize(num_rows_A, dim);
                  AT.SetSize(dim, num_rows_A);
                  W.SetSize(num_rows_A);
                  WA.SetSize(num_rows_A, dim);
                  RHS.SetSize(num_rows_A);
                  RHS_mod.SetSize(dim);
                  A = 0., AT = 0., ATWA = 0., W = 0., WA = 0., RHS = 0., RHS_mod = 0.;

                  /*** Fill matrix representing overdetermined system and RHS vector ***/
                  // Iterate over adjacent cells
                  for (int element_it = 0; element_it < length_element_row; element_it++)
                  {
                     int el_index = element_row[element_it];

                     pmesh->GetElementEdges(el_index, cell_faces_row, oris);
                     length_cell_faces = cell_faces_row.Size();

                     // Iterate over cell faces
                     for (int face_it = 0; face_it < length_cell_faces; face_it++)
                     {
                        int face_index = cell_faces_row[face_it];

                        // Check if cell face is unique
                        if (uniqueFaceIndices.find(face_index) == uniqueFaceIndices.end())
                        {
                           // New face, fill corresponding row in A and RHS
                           uniqueFaceIndices.insert(face_index);
                           row_index = uniqueFaceIndices.size() - 1;

                           // Get adjacent cell information for the weighting
                           FI = pmesh->GetFaceInformation(face_index);
                           c = FI.element[0].index;
                           cp = FI.element[1].index;
                           
                           // Get face velocity, outward normal, and face position
                           // n_vec is the outward normal times |F|
                           GetIntermediateFaceVelocity(face_index, Vf);
                           CalcOutwardNormalInt(S, el_index, face_index, n_vec);

                           // Compute weight
                           double rhoc = 0., rhocp = 0.;
                           GetCellStateVector(S, c, Uc);
                           
                           // rhoc = 1. / Uc[0];
                           // if (FI.IsInterior())
                           // {
                           //    GetCellStateVector(S, cp, Ucp);
                           //    rhocp = 1 / Ucp[0];
                           // }
                           // double weight = rhoc + rhocp;
                           double weight = 1.;

                           // Fill rows of A
                           A.SetRow(row_index, n_vec);
                           // Fill weight matrix
                           W(row_index, row_index) = weight;
                           // Fill rows of RHS
                           double val = n_vec * Vf;
                           RHS(row_index) = weight * val;
                           // cout << endl;
                        } // End unique new face
                     } // End face iterator
                  } // End cell iterator

                  // Matrix setup
                  AT = A;
                  AT.Transpose();
                  Mult(W, A, WA);
                  Mult(AT,WA,ATWA);
                  AT.Mult(RHS, RHS_mod);

                  // Solve for best fit
                  ATWA.Invert();
                  ATWA.Mult(RHS_mod, node_v);

                  // Compute residual
                  // temp_vec.SetSize(num_rows_A);
                  // A.Mult(node_v, temp_vec);
                  // temp_vec -= RHS;
                  // cout << "Boundary/Corner Node: " << node << ", residual: " << temp_vec.Norml2() << endl;

                  break;
               } // End corner and boundary vertex
               case 4:
               {
                  num_rows_A = 12;

                  A.SetSize(num_rows_A, num_cols_A);
                  RHS.SetSize(num_rows_A);
                  temp_vec.SetSize(num_rows_A);

                  // cout << "===== new corner node ===== [at: " << node << "]\n";
                  /*** Fill matrix representing overdetermined system and RHS vector ***/
                  // Iterate over adjacent cells
                  for (int element_it = 0; element_it < length_element_row; element_it++)
                  {
                     int el_index = element_row[element_it];

                     pmesh->GetElementEdges(el_index, cell_faces_row, oris);
                     length_cell_faces = cell_faces_row.Size();

                     // Iterate over cell faces
                     for (int face_it = 0; face_it < length_cell_faces; face_it++)
                     {
                        int face_index = cell_faces_row[face_it];

                        // Check if cell face is unique
                        if (uniqueFaceIndices.find(face_index) == uniqueFaceIndices.end())
                        {
                           // New face, fill corresponding row in A and RHS
                           uniqueFaceIndices.insert(face_index);
                           row_index = uniqueFaceIndices.size() - 1;
                           
                           // Get face velocity, outward normal, and face position
                           // n_vec is the outward normal times |F|
                           GetIntermediateFaceVelocity(face_index, Vf);
                           CalcOutwardNormalInt(S, el_index, face_index, n_vec);

                           // Get average of face edge nodes
                           H1.GetFaceDofs(face_index, face_dofs);
                           GetNodePosition(S, face_dofs[1], vdof1_x);
                           GetNodePosition(S, face_dofs[0], vdof2_x);
                           face_x = vdof1_x;
                           face_x.Add(1., vdof2_x);
                           face_x *= 0.5;

                           // Fill rows of A
                           A(row_index, 0) = n_vec(0);
                           A(row_index, 1) = n_vec(1);
                           A(row_index, 2) = n_vec(0)*face_x(0);
                           A(row_index, 3) = n_vec(1)*face_x(0);
                           A(row_index, 4) = n_vec(0)*face_x(1);
                           A(row_index, 5) = n_vec(1)*face_x(1);

                           // Fill rows of RHS
                           RHS(row_index) = Vf * n_vec;
                        } // End unique new face
                     } // End face iterator
                  } // End cell iterator

                  // Solve for vgeo on this node
                  AT = A;
                  AT.Transpose();
                  Mult(AT, A, ATA);
                  AT.Mult(RHS, RHS_mod);
                  ATA.Invert();
                  ATA.Mult(RHS_mod,res);

                  // Compute residual
                  // A.Mult(res, temp_vec);
                  // temp_vec -= RHS;
                  // cout << "Interior Node: " << node << ", residual: " << temp_vec.Norml2() << endl;

                  // Set node velocity from coefficient vector
                  node_v(0) = res[0] + res[2]*node_x[0] + res[4]*node_x[1];
                  node_v(1) = res[1] + res[3]*node_x[0] + res[5]*node_x[1];
                  break;
               }
               default:
               {
                  MFEM_ABORT("Invalid number of adjacent ceells.\n");
                  break;
               }
            } // end switch between corner, boundary, and interior vertex

            // Only update the v_geo_gf in the case of the corners
            // Using Q1 for this velocity field
            SetViGeo(node, node_v);

            break;
         } // End corner
         case 1: // face
         {
            // face nodes just set to ivf
            assert(face >= 0);
            GetIntermediateFaceVelocity(EDof, node_v);

            break;
         } // end face
         case 2: // Cell Center
         {
            Vector Uc(dim+2);
            GetCellStateVector(S, EDof, Uc);
            pb->velocity(Uc, node_v);

            break;
         } // End cell center
         default:
         {
            MFEM_ABORT("Invalid entity\n");
         }
      } // end entity switch

      // In every case, update the corresponding nodal velocity in S
      UpdateNodeVelocity(S, node, node_v);

   } // End node iterator
}


/****************************************************************************************************
* Function: ComputeGeoVCAVEATCellFaceWeighted
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function computes the geometric velocity from the face normals and face 
*  velocities from all adjacent faces.
*
*  On the boundary vertices, the above CAVEAT algorithm is used, which is 
*  described in https://www.osti.gov/biblio/5366795.
*
*  On the interior vertices, the above adjacent Cell face based method is used.
*
*  Note: Vigeo is only the geometric velocity, and is not the velocity used 
*        to move the mesh.  For the mesh velocity, one must solve for the 
*        alpha_i, which process is found in ComputeNodeVelocitiesFromVgeo().
*
*  Note: This function fills the Q1 velocity field v_geo_gf and also the mv_gf
*        in the BlockVector S.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeGeoVCAVEATCellFaceWeighted(Vector &S)
{
   // cout << "ComputeGeoVCAVEATCellFaceWeighted call\n";
   mfem::Mesh::FaceInformation FI;
   
   Vector node_v(dim), node_x(dim), n_vec(dim), Vf(dim);
   Vector vdof1_x(dim), vdof2_x(dim), face_x(dim), vdof_avg(dim);
   double F, c_norm;

   Array<int> element_row, cell_faces_row, oris, face_dofs;
   int row_index, length_element_row, length_cell_faces, num_rows_A, num_cols_A = 6;
   int c, cp;
   Vector Uc(dim + 2), Ucp(dim+2);

   DenseMatrix A, AT, ATA(num_cols_A), W, WA, ATWA(dim);
   Vector RHS, RHS_mod(num_cols_A), res(num_cols_A), temp_vec;

   DofEntity entity;
   int EDof;

   std::unordered_set<int> uniqueFaceIndices;

   // Iterate over all mesh velocity dofs
   for (int node = 0; node < NDofs_H1; node++) // Vertex iterator
   {
      node_v = 0.; // Reset node velocity
      GetEntityDof(node, entity, EDof);

      switch (entity)
      {
         case 0: // corner
         {
            // Reset vals
            uniqueFaceIndices.clear();
            F = 0.;
            node_v = 0.;
            GetNodePosition(S, node, node_x);

            // Get adjacent cells
            vertex_element->GetRow(node, element_row);
            length_element_row = element_row.Size();

            // Set size of overdetermined matrix A depending on number of adjacent cells
            switch (length_element_row)
            {
               case 1: // Domain corner vertex
               case 2: // Boundary vertex
               {
                  // Since in the Cell Face based method we get a singular matrix for
                  // the domain corner and boundary vertices, we use the CAVEAT method
                  // of projecting the vertex velocity onto the adjacent cell faces
                  // and solve the system using weighted least squares.
                  if (length_element_row == 1) {
                     num_rows_A = 4;
                  } else {
                     num_rows_A = 7;
                  }

                  // Set sizes of nodal matrices
                  A.SetSize(num_rows_A, dim);
                  AT.SetSize(dim, num_rows_A);
                  W.SetSize(num_rows_A);
                  WA.SetSize(num_rows_A, dim);
                  RHS.SetSize(num_rows_A);
                  RHS_mod.SetSize(dim);
                  A = 0., AT = 0., ATWA = 0., W = 0., WA = 0., RHS = 0., RHS_mod = 0.;

                  /*** Fill matrix representing overdetermined system and RHS vector ***/
                  // Iterate over adjacent cells
                  for (int element_it = 0; element_it < length_element_row; element_it++)
                  {
                     int el_index = element_row[element_it];

                     pmesh->GetElementEdges(el_index, cell_faces_row, oris);
                     length_cell_faces = cell_faces_row.Size();

                     // Iterate over cell faces
                     for (int face_it = 0; face_it < length_cell_faces; face_it++)
                     {
                        int face_index = cell_faces_row[face_it];

                        // Check if cell face is unique
                        if (uniqueFaceIndices.find(face_index) == uniqueFaceIndices.end())
                        {
                           // New face, fill corresponding row in A and RHS
                           uniqueFaceIndices.insert(face_index);
                           row_index = uniqueFaceIndices.size() - 1;

                           // Get adjacent cell information for the weighting
                           FI = pmesh->GetFaceInformation(face_index);
                           c = FI.element[0].index;
                           cp = FI.element[1].index;
                           
                           // Get face velocity, outward normal, and face position
                           // n_vec is the outward normal times |F|
                           GetIntermediateFaceVelocity(face_index, Vf);
                           CalcOutwardNormalInt(S, el_index, face_index, n_vec);

                           // Compute weight
                           double rhoc = 0., rhocp = 0.;
                           GetCellStateVector(S, c, Uc);
                           
                           rhoc = 1. / Uc[0];
                           if (FI.IsInterior())
                           {
                              GetCellStateVector(S, cp, Ucp);
                              rhocp = 1. / Ucp[0];
                           }
                           double weight = rhoc + rhocp;
                           // double weight = 1.;

                           // Fill rows of A
                           A.SetRow(row_index, n_vec);
                           // Fill weight matrix
                           W(row_index, row_index) = weight;
                           // Fill rows of RHS
                           double val = n_vec * Vf;
                           RHS(row_index) = weight * val;
                           // cout << endl;

                           // Print stuff
                           // cout << "face: " << face << ", vf: ";
                           // Vf.Print(cout);
                           // cout << "normal: ";
                           // n_vec.Print(cout);
                           // cout << "weight: " << weight << endl;
                        } // End unique new face
                     } // End face iterator
                  } // End cell iterator

                  // Matrix setup
                  AT = A;
                  AT.Transpose();
                  Mult(W, A, WA);
                  Mult(AT,WA,ATWA);
                  AT.Mult(RHS, RHS_mod);

                  // Solve for best fit
                  ATWA.Invert();
                  ATWA.Mult(RHS_mod, node_v);

                  // Compute residual
                  // temp_vec.SetSize(num_rows_A);
                  // A.Mult(node_v, temp_vec);
                  // temp_vec -= RHS;
                  // cout << "Boundary/Corner Node: " << node << ", residual: " << temp_vec.Norml2() << endl;

                  break;
               } // End corner and boundary vertex
               case 4:
               {
                  num_rows_A = 12;

                  A.SetSize(num_rows_A, num_cols_A);
                  RHS.SetSize(num_rows_A);
                  temp_vec.SetSize(num_rows_A);

                  // cout << "===== new corner node ===== [at: " << node << "]\n";
                  /*** Fill matrix representing overdetermined system and RHS vector ***/
                  // Iterate over adjacent cells
                  for (int element_it = 0; element_it < length_element_row; element_it++)
                  {
                     int el_index = element_row[element_it];

                     pmesh->GetElementEdges(el_index, cell_faces_row, oris);
                     length_cell_faces = cell_faces_row.Size();

                     // Iterate over cell faces
                     for (int face_it = 0; face_it < length_cell_faces; face_it++)
                     {
                        int face_index = cell_faces_row[face_it];

                        // Check if cell face is unique
                        if (uniqueFaceIndices.find(face_index) == uniqueFaceIndices.end())
                        {
                           // New face, fill corresponding row in A and RHS
                           uniqueFaceIndices.insert(face_index);
                           row_index = uniqueFaceIndices.size() - 1;
                           
                           // Get face velocity, outward normal, and face position
                           // n_vec is the outward normal times |F|
                           GetIntermediateFaceVelocity(face_index, Vf);
                           CalcOutwardNormalInt(S, el_index, face_index, n_vec);

                           // Compute weight
                           double rhoc = 0., rhocp = 0., weight = 0.;
                           GetCellStateVector(S, c, Uc);
                           
                           rhoc = 1. / Uc[0];
                           if (FI.IsInterior())
                           {
                              GetCellStateVector(S, cp, Ucp);
                              rhocp = 1. / Ucp[0];
                              weight = rhoc + rhocp;
                           }
                           else 
                           {
                              weight = rhoc;
                           }
                           
                           // double weight = 1.;

                           // Get average of face edge nodes
                           H1.GetFaceDofs(face_index, face_dofs);
                           GetNodePosition(S, face_dofs[1], vdof1_x);
                           GetNodePosition(S, face_dofs[0], vdof2_x);
                           face_x = vdof1_x;
                           face_x.Add(1., vdof2_x);
                           face_x *= 0.5;

                           // Fill rows of A
                           A(row_index, 0) = n_vec(0);
                           A(row_index, 1) = n_vec(1);
                           A(row_index, 2) = n_vec(0)*face_x(0);
                           A(row_index, 3) = n_vec(1)*face_x(0);
                           A(row_index, 4) = n_vec(0)*face_x(1);
                           A(row_index, 5) = n_vec(1)*face_x(1);

                           // Fill weight matrix
                           W(row_index, row_index) = weight;

                           // Fill rows of RHS
                           double rhs_val = Vf * n_vec;
                           rhs_val *= weight;
                           RHS(row_index) = rhs_val;

                           // Print stuff
                           // cout << "face: " << face << ", vf: ";
                           // Vf.Print(cout);
                           // cout << "normal: ";
                           // n_vec.Print(cout);
                           // cout << "weight: " << weight << endl;
                        } // End unique new face
                     } // End face iterator
                  } // End cell iterator

                  // Matrix setup
                  AT = A;
                  AT.Transpose();
                  Mult(W, A, WA);
                  Mult(AT,WA,ATWA);
                  AT.Mult(RHS, RHS_mod);

                  // Solve for best fit
                  ATWA.Invert();
                  ATWA.Mult(RHS_mod, res);

                  // Compute residual
                  // A.Mult(res, temp_vec);
                  // temp_vec -= RHS;
                  // cout << "Interior Node: " << node << ", residual: " << temp_vec.Norml2() << endl;

                  // Set node velocity from coefficient vector
                  node_v(0) = res[0] + res[2]*node_x[0] + res[4]*node_x[1];
                  node_v(1) = res[1] + res[3]*node_x[0] + res[5]*node_x[1];
                  break;
               }
               default:
               {
                  MFEM_ABORT("Invalid number of adjacent ceells.\n");
                  break;
               }
            } // end switch between corner, boundary, and interior vertex

            // Only update the v_geo_gf in the case of the corners
            // Using Q1 for this velocity field
            SetViGeo(node, node_v);

            break;
         } // End corner
         case 1: // face
         {
            // face nodes just set to ivf
            assert(face >= 0);
            GetIntermediateFaceVelocity(EDof, node_v);

            break;
         } // end face
         case 2: // Cell Center
         {
            Vector Uc(dim+2);
            GetCellStateVector(S, EDof, Uc);
            pb->velocity(Uc, node_v);

            break;
         } // End cell center
         default:
         {
            MFEM_ABORT("Invalid entity\n");
         }
      } // end entity switch

      // In every case, update the corresponding nodal velocity in S
      UpdateNodeVelocity(S, node, node_v);

   } // End node iterator
}


/****************************************************************************************************
* Function: IterativeCornerVelocityMC
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function constitutes one iteration on the previously computed corner node 
*  velocity to reduce total face motion while still preserving mass conservation.
*  The idea is that we correct the velocity at the vertices to ensure mass 
*  conservation and dampen the parachuting effect seen at the faces. This method 
*  was proposed on 04/03/2024.
*
*  Note: This function relies on the mesh velocities having already been computed
*  and linearized.  One can use whichever mesh velocity computation before 
*  iteration.
*
*  Note: This function will modify mv_gf in the BlockVector S
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::IterativeCornerVelocityMC(Vector &S, const double & dt)
{
   // cout << "=====IterativeCornerVelocityMC=====\n";
   bool do_theta_averaging = false;
   double theta = 1.;

   // We need a palce to store the new mesh velocities that we compute
   Vector S_new = S;

   // Values needed during iteration
   Array<int> faces_row, face_dofs_row;
   int faces_length;
   mfem::Mesh::FaceInformation FI;
   H1.ExchangeFaceNbrData();
   int Vadj_index, c;
   bool is_v2 = false, int_face_bdr_node = false;
   double const1 = 0., const2 = 0., F = 0.;

   Vector predicted_node_v(dim), Vf(dim), Vadj(dim), Vnode(dim), Vnode_prev_it(dim);
   double Vnode_n_comp = 0.; // The quantity we solve for in each face iteration
   double Vnode_prev_it_nR_comp = 0., Vadj_n_comp = 0., Vadj_nR_comp = 0.;
   Vector n_int(dim), n_vec(dim), n_vec_R(dim);
   Vector Vnode_x(dim), Vadj_x(dim), face_x(dim);
   Vector Badj(dim), Bnode(dim);

   // Averaging
   DenseMatrix Mi(dim), dm_tmp(dim);
   Vector Ri(dim), v_tmp(dim);
   Vector Vghosted(dim);

   int bdr_ind = 0;

   /* Iterate over corner nodes */
   for (int node = 0; node < NDofs_H1L; node++) // TODO: Is NDofs_H1L == NVDofs_H1?
   {
      // Reset new nodal velocity, and averaging objects
      Mi = 0., Ri = 0.;
      predicted_node_v = 0.;

      // Set corresponding boundary indicator
      // If node is boundary node, will have to add corresponding correction
      // node is interior node if bdr_ind == -1
      bdr_ind = BdrVertexIndexingArray[node];
      // cout << "\t node: " << node << ", bdr_ind: " << bdr_ind << endl;
      // if (bdr_ind != -1) { continue; }

      // Get index of face to be reflected

      // Get current nodal velocity from S
      // GetNodeVelocity(S, node, node_v);
      // Get nodal velocity at previous iteration from mv_gf_prev_it
      for (int i = 0; i < dim; i++)
      {
         int index = node + i * NDofs_H1;
         Vnode_prev_it[i] = mv_gf_prev_it[index];
      }

      /* Get node position */
      GetNodePosition(S, node, Vnode_x);

      // Get cell faces
      vertex_edge.GetRow(node, faces_row);
      faces_length = faces_row.Size();

      /* Iterate over cell faces */
      for (int face_it = 0; face_it < faces_length; face_it++) // Adjacent face iterator
      {
         // Reset flag
         is_v2 = false;

         // Get face information
         int face = faces_row[face_it];
         GetIntermediateFaceVelocity(face, Vf);
         FI = pmesh->GetFaceInformation(face);

         // cout << "face: " << face << endl;

         // Calculate outer normal
         c = FI.element[0].index;
         CalcOutwardNormalInt(S, c, face, n_int);
         n_vec = n_int;
         F = n_vec.Norml2();
         n_vec /= F;
         n_vec_R = n_vec;
         Orthogonal(n_vec_R);

         // Calculate bmn
         double bmn = Vf * n_vec;
         bmn *= F;

         // Get Vnode_prev_it component in tangent direction from previous iteration
         Vnode_prev_it_nR_comp = Vnode_prev_it * n_vec_R;
         // if (abs(Vnode_prev_it_nR_comp) > 0.)
         // {
         //    cout << "Vnode_prev_it: ";
         //    Vnode_prev_it.Print(cout);
         //    cout << "n_vec_R: ";
         //    n_vec_R.Print(cout);
         //    cout << "Vnode_prev_it_nR_comp: " << Vnode_prev_it_nR_comp << endl;
         // }

         /* adjacent corner indices */
         H1.GetFaceDofs(face, face_dofs_row);
         int face_vdof1 = face_dofs_row[1], 
             face_vdof2 = face_dofs_row[0], 
             face_dof = face_dofs_row[2]; // preserve node orientation where cell to right of face is with lower cell index
         // cout << "face: " << face << endl
         //      << "interior cell: " << c << endl
         //      << "normal: ";
         // n_vec.Print(cout);
         // cout << "nR: ";
         // n_vec_R.Print(cout); 
         // cout << ", face_vdof1: " << face_vdof1 
         //      << ", face_vdof2: " << face_vdof2 
         //      << ", face_dof: " << face_dof << endl;

         // Grab corresponding vertex velocity from S
         if (node == face_vdof1) {
            Vadj_index = face_vdof2;
         } else {
            Vadj_index = face_vdof1;
            is_v2 = true;
         }
         GetNodeVelocity(S, Vadj_index, Vadj);

         /* Get adj node and face locations */
         GetNodePosition(S, Vadj_index, Vadj_x);
         GetNodePosition(S, face_dof, face_x);

         // Get normal and rotated components of Vadj
         Vadj_n_comp = Vadj * n_vec;
         Vadj_nR_comp = Vadj * n_vec_R;

         // Compute geometrical vectors
         Bnode = face_x;
         Bnode *= -2.;
         Bnode.Add(-.5, Vadj_x);
         Bnode.Add(2.5, Vnode_x);
         Orthogonal(Bnode);

         Badj = face_x;
         Badj *= -2.;
         Badj.Add(-.5, Vnode_x);
         Badj.Add(2.5, Vadj_x);
         Orthogonal(Badj);

         // Compute dot products
         double badjn = Badj * n_vec;
         double badjnr = Badj * n_vec_R;
         double bnoden = Bnode * n_vec;
         double bnodenr = Bnode * n_vec_R;

         // Evaluate numerator, which depends on geometric orientation
         // The center of the cell should be on the right when 
         // traversing from node 1 to node 2
         double numer = 0.;

         // Solve for Vnode_n_comp
         if (is_v2) {
            numer += 3*bmn; 
         } else {
            numer -= 3*bmn;  
         }
         numer += (badjn - dt * Vadj_nR_comp) * Vadj_n_comp;
         numer += badjnr * Vadj_nR_comp;
         numer += ((dt/2) * Vadj_n_comp - bnodenr)*Vnode_prev_it_nR_comp;

         double denom = -dt*Vnode_prev_it_nR_comp + (dt/2)*Vadj_nR_comp + bnoden;
   
         Vnode_n_comp = numer / denom;

         Vnode = 0.;
         Vnode.Add(Vnode_n_comp, n_vec);
         Vnode.Add(Vnode_prev_it_nR_comp, n_vec_R);



         // Add contribution to the averaging objects
         tensor(n_vec, n_vec, dm_tmp);
         dm_tmp.Mult(Vnode, v_tmp);
         Mi.Add(1., dm_tmp);
         Ri.Add(1., v_tmp);

         // In the case of boundary nodes, check if this is an interior face
         if (bdr_ind != -1 && FI.IsInterior())
         {
            Vghosted = Vnode;
            // We must add the ghost node contribution
            // cout << "we have an interior face " << face << " on a boundary node " << node << endl;
            // What is done here depends on which boundary face we have 
            // We will be employing reflective boundary conditions

            // Compute vector from node to face 
            Vector face_dx(dim);
            subtract(Vnode_x, face_x, face_dx);

            // Compute vector from node to adjacent node
            Vector Vadj_dx(dim);
            subtract(Vnode_x, Vadj_x, Vadj_dx);

            // Sod
            switch (bdr_ind)
            {
            case 1: // bottom
            case 3: // top
               n_vec[1] *= -1.;
               Vadj[1] *= -1.;
               face_dx[0] *= -1.;
               Vadj_dx[0] *= -1.;
               break;
            
            case 2: // right
            case 4: // left
               n_vec[0] *= -1.;
               Vadj[0] *= -1.;
               face_dx[1] *= -1.;
               Vadj_dx[1] *= -1.;
               break;
            
            case -1:
               // Not a boundary vertex
               continue;

            default:
               MFEM_ABORT("Incorrect bdr attribute encountered while enforcing mesh velocity BCs.\n");
               break;
            }

            n_vec_R = n_vec;
            Orthogonal(n_vec_R);
            Vnode_prev_it_nR_comp = Vnode_prev_it * n_vec_R;
            // Negate is_v2
            is_v2 = !is_v2;

            // Get normal and rotated components of Vadj
            Vadj_n_comp = Vadj * n_vec;
            Vadj_nR_comp = Vadj * n_vec_R;

            /* Get flipped face_x and vadj_x */
            add(Vnode_x, face_dx, face_x);
            add(Vnode_x, Vadj_dx, Vadj_x);

            Bnode = face_x;
            Bnode *= -2.;
            Bnode.Add(-.5, Vadj_x);
            Bnode.Add(2.5, Vnode_x);
            Orthogonal(Bnode);

            Badj = face_x;
            Badj *= -2.;
            Badj.Add(-.5, Vnode_x);
            Badj.Add(2.5, Vadj_x);
            Orthogonal(Badj);

            double badjn = Badj * n_vec;
            double badjnr = Badj * n_vec_R;
            double bnoden = Bnode * n_vec;
            double bnodenr = Bnode * n_vec_R;

            double numer = 0.;

            // Solve for Vnode_n_comp
            if (is_v2) {
               numer += 3*bmn; 
            } else {
               numer -= 3*bmn;  
            }
            numer += (badjn - dt * Vadj_nR_comp) * Vadj_n_comp;
            numer += badjnr * Vadj_nR_comp;
            numer += ((dt/2) * Vadj_n_comp - bnodenr)*Vnode_prev_it_nR_comp;

            double denom = -dt*Vnode_prev_it_nR_comp + (dt/2)*Vadj_nR_comp + bnoden;
      
            Vnode_n_comp = numer / denom;

            Vnode = 0.;
            Vnode.Add(Vnode_n_comp, n_vec);
            Vnode.Add(Vnode_prev_it_nR_comp, n_vec_R);

            // Add contribution to the averaging objects
            tensor(n_vec, n_vec, dm_tmp);
            dm_tmp.Mult(Vnode, v_tmp);
            Mi.Add(1., dm_tmp);
            Ri.Add(1., v_tmp);
            
         } // End ghost node
      }

      // Average node_v
      Mi.Invert();
      Mi.Mult(Ri, predicted_node_v);
      
      // cout << "predicted velocity at node " << node << ": ";
      // predicted_node_v.Print(cout);

      // Theta average with previous velocity
      if (do_theta_averaging)
      {
         predicted_node_v *= theta; 
         predicted_node_v.Add(1. - theta, Vnode_prev_it);
      }
      

      // Put velocity in S_new
      UpdateNodeVelocity(S_new, node, predicted_node_v);
   }

   // Once all the predicted node velocities have been computed, set them in S
   S = S_new;
}


/****************************************************************************************************
* Function: ComputeIterationNorm
* Parameters:
*  S        - BlockVector representing FiniteElement information
*  dt       - Current time step
*
* Purpose:
* 
*  NOTE: Interior faces. 
*  
****************************************************************************************************/
template<int dim>
double LagrangianLOOperator<dim>::ComputeIterationNorm(Vector &S, const double & dt)
{
   double val = 0., denom_val = 0.;
   int num_broken = 0, total_num = 0;
   
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

      // if (FI.IsInterior())
      // {
      // Get adjacent cell information
      c = FI.element[0].index;
      cp = FI.element[1].index;
      GetCellStateVector(S, c, Uc);
      GetCellStateVector(S, cp, Ucp);
      pb->velocity(Uc, cell_c_v);
      pb->velocity(Ucp, cell_cp_v);

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
      // subtract(vdof2_v, vdof1_v, temp_vec); // only change temp_vec, since temp_vec_2 is same from D calculation (half step representation)
      // double c1 = ( dt * (temp_vec * n_vec) + 2. * (temp_vec_2 * n_vec_R) ) / D; // TRYING SOMETHING HERE. WILL NEED CORRECTED

      // Compute c0 (A.4b)
      // Vector n_vec_half(dim);
      // subtract(vdof2_x_half, vdof1_x_half, n_vec_half);
      // Orthogonal(n_vec_half);
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
      double c0 = (3. / D) * (bmn + const1 / 2. + const2 / 6. + 2. * const3 / 3.);
      
      // Add to val for interior faces only
      add(vdof1_v, vdof2_v, temp_vec);
      double fval = (temp_vec * n_vec) / 2.;
      denom_val += abs(fval);
      // if (fval != c0)
      // {
      //    cout << "face: " << face 
      //         << ", fval: " << fval 
      //         << ", c0: " << c0 << endl;
      // }
      fval -= c0;

      if (abs(fval) > 1.e-8)
      {
         num_broken++;
         // cout << "face " << face << " is broken. fval: " << fval << ".\n";

      }

      val += abs(fval);
      total_num++;
      // }
   }
   double perc_broken = double(num_broken)/double(total_num);
   // cout << "percentage of broken faces: " << perc_broken << endl;
   return val / denom_val;
}


/****************************************************************************************************
* Function: ComputeNodeVelocitiesFromVgeo
* Parameters:
*  S        - BlockVector representing FiniteElement information
*  t        - Current time
*  dt       - Current time step
*  flag     - Flag used to indicate testing, can be "testing" or "NA"
*  test_vel - Velocity used for testing
*
* Purpose:
*  This function iterates over all Serendipity nodes and computes the 
*  velocity with which to move the mesh from the geometric velocity.
*  
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeNodeVelocitiesFromVgeo(
   Vector &S, const double & t, double & dt, const string, 
   void (*test_vel)(const Vector&, const double&, Vector&))
{
   Vector node_v(dim);
   bool is_dt_changed = false;

   // Iterate over vertices and faces
   for (int node = 0; node < NVDofs_H1; node++) // Vertex and face iterator
   {
      ComputeNodeVelocityFromVgeo(S, node, dt, node_v, is_dt_changed);

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
* Function: ComputeNodeVelocityFromVgeo
* Parameters:
*  S             - BlockVector representing FiniteElement information
*  node          - Global DOF to compute velocity at
*  dt            - Current time step
*  node_v        - The computed velocity to be returned
*  is_dt_changed - Boolean representing if the timestep was modified
*                  due to a timestep restriction from the computation
*                  of alpha_i
*
* Purpose:
*  This function is responsible for computing the nodal velocity
*  from the geometric representation of the velocity through the
*  algorithm outlined in the paper:
*     [alpha_i I - (dt/2) C_i]^-1 Vigeo.
*  
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeNodeVelocityFromVgeo(
   Vector &S, const int & node, double & dt, Vector &node_v, bool &is_dt_changed)
{
   // cout << "=======================================\n"
   //      << "      ComputeNodeVelocityFromVgeo      \n"
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

         // chrono_temp.Clear();
         // chrono_temp.Start();
         ComputeCiGeo(node, Ci);
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

         GetNodeVelocity(S, node, Vgeo);
         
         // chrono_temp.Clear();
         // chrono_temp.Start();
         // cout << "node: " << node << endl;
         ComputeDeterminant(Ci, dt, d);
         // chrono_temp.Stop();

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
* Function: Raviart-Thomas Velocity Functions
****************************************************************************************************/

/****************************************************************************************************
* Function: ComputeGeoVRaviart
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function computes the geometric velocity using the Raviart-Thomas elements,
*  which are face-based vector finite elements.  This function must take into account
*  the Piola transformations, which lucky for us MFEM will handle by use of the 
*  ElementTransformation class.
*
*  Note: Vigeo is only the geometric velocity, and is not the velocity used 
*        to move the mesh.  For the mesh velocity, one must solve for the 
*        alpha_i, which process is found in ComputeNodeVelocitiesFromVgeo().
*
*  Note: This function fills the Q1 velocity field v_geo_gf and also the mv_gf
*        in the BlockVector S.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeGeoVRaviart(Vector &S)
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

   for (int node = 0; node < NDofs_H1 - NDofs_L2; node++) // Vertex iterator
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
            pb->velocity(Uc, node_v);

            /* Put this velocity into mv_gf */
            UpdateNodeVelocity(S,node,node_v);
            
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

            // GetCorrectedFaceFlux(face_index, ivf); 
            GetIntermediateFaceVelocity(face_index, ivf);
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
      if (entity == 0) // Only for corners
      {
         SetViGeo(node, node_v);
      }

      UpdateNodeVelocity(S, node, node_v);
      
   } // End vertex iterator
}


/****************************************************************************************************
* Function: ComputeGeoVRaviart2
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  Iterate over cells and average contribution to nodes. Alternate to above
*  function.
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
            int index = ldof + i * NVDofs_H1;
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
* Function: ComputeCiGeo
* Parameters:
*  node  - Global DOF corresponding to the node the CiGeo should be evaluated at.
*          DOF should correspond to a corner or face node, as there is no reason 
*          evaluate CiGeo at a cell center.
*  res   - DenseMatrix representing CiGeo evaluated at the provided node.
*
* Purpose:
*  This function serves to evaluate CiGeo at a given global dof.  Recall that
*  CiGeo is defined as
*     CiGeo =  [1 / Sum |K|] Sum int grad v_geo_gf
*
*  This function is essential to go from ViGeo to Vi.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeCiGeo(const int & node, DenseMatrix & res)
{
   // cout << "=======================================\n"
   //      << "             ComputeCiGeo              \n"
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

   switch (entity)
   {
      case 0: // corner
      {
         vertex_element->GetRow(node, row);
         break;
      }
      case 1: // face
      {
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

   int row_length = row.Size();
   for (int row_it = 0; row_it < row_length; row_it++)
   {
      int row_el = row[row_it];
      double cell_vol = pmesh->GetElementVolume(row_el);
      denom += cell_vol;

      dm_temp = 0.;
      IntGrad(row_el, dm_temp);
      res.Add(cell_vol, dm_temp);
   }
   res *= 1./denom;
}


/****************************************************************************************************
* Function: IntGrad
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
void LagrangianLOOperator<dim>::IntGrad(const int cell, DenseMatrix & res)
{
   // cout << "=======================================\n"
   //      << "             IntGrad            \n"
   //      << "=======================================\n";

   ParGridFunction H1Lc_gf(&H1Lc);
   const int _size = H1Lc.GetVSize();

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
         H1Lc_gf.MakeRef(&H1Lc, v_geo_gf, j*_size);
         H1Lc_gf.GetGradient(*trans, grad);

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
