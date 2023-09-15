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
   cout << "v_CR_gf.Size before: " << v_CR_gf.Size() << endl;

   // resize v_CR_gf to correspond to the number of faces
   if (dim == 1)
   {
      v_CR_gf.SetSize(num_faces);
      lambda_max_vec.SetSize(num_faces);
   }
   
   // Initialize values of intermediate face velocities
   v_CR_gf = 0.;
   // assert(v_CR_gf.Size() == dim * num_faces);
   // cout << "v_CR_gf.Size after: " << v_CR_gf.Size() << endl;

   // Initialize Dij sparse
   InitializeDijMatrix();   

   // Set integration rule for Rannacher-Turek space
   IntegrationRules IntRulesLo(0, Quadrature1D::GaussLobatto);
   RT_ir = IntRules.Get(CR.GetFE(0)->GetGeomType(), RT_ir_order);

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
   mfem::Mesh::FaceInformation FI;
   int c, cp;
   Vector Uc(dim+2), Ucp(dim+2), n_int(dim), c_vec(dim), n_vec(dim);
   double F, lambda_max, d;

   for (int face = 0; face < num_faces; face++) // face iterator
   {
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

         // Compute max wave speed
         // double b_covolume = .1 / (max(1./Uc[0], 1./Ucp[0]));
         double b_covolume = 1.;
         double pl = pb->pressure(Uc);
         double pr = pb->pressure(Ucp);

         lambda_max = pb->compute_lambda_max(Uc, Ucp, n_vec, pl, pr, b_covolume);

         d = lambda_max * c_norm; 

         double ss = pb->sound_speed(Uc);

         cout << "c: " << c 
              << ", pl: " << pl 
              << ", pr: " << pr
              << ", cp: " << cp << endl;
         cout << std::setprecision(12) 
              << "sound speed: " << ss 
              << ", lambda_max: " << lambda_max 

              << ", c_norm: " << c_norm 
              << ", d: " << d << endl;

         lambda_max_vec[face] = lambda_max; // TODO: remove, only temporary
         dij_sparse->Elem(c,cp) = d;
         dij_sparse->Elem(cp,c) = d;
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
   cout << "CalculateTimestep\n";
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
      cout << "\tcell: " << ci << endl;
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

            cout << "face: " << fids[j] << endl;
            cout << "d for ci " << ci << " and cj " << cj << ": " << d << endl;

            temp_sum += d;
         }
      }
      
      t_temp = 0.5 * ((CFL * mi) / temp_sum );

      cout << "CFL: " << CFL << ", mi: " << mi << ", temp_sum: " << temp_sum << endl;
      cout << "t_temp: " << t_temp << ", t_min: " << t_min << endl;

      if (t_temp < t_min && t_temp > 1.e-12) { 
         cout << "timestep reduced\n";
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
   cout << "========================================\n"
        << "MakeTimeStep\n"
        << "========================================\n";
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
            if (pb->get_indicator() == "noh")
            {
               // Since we will enforce Dirichlet conditions on all boundary cells, 
               // this enforcement will be done after the face iterator
            }
            else if (pb->get_indicator() == "saltzmann")
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
         if (pb->get_indicator() == "IsentropicVortex" || pb->get_indicator() == "Noh")
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
   cout << "========================================\n"
        << "EnforceExactBCOnCell\n"
        << "========================================\n";

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
   cout << "========================================\n"
        << "EnforcingMVBoundaryConditions\n"
        << "========================================\n";

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
   cout << "=======================================\n"
        << "   ComputeIntermediateFaceVelocities   \n"
        << "=======================================\n";
   
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
            GetCellStateVector(S, cp, Ucp);

            // Get normal, d, and |F|
            CalcOutwardNormalInt(S, c, face, n_int);
            n_vec = n_int;
            double F = n_vec.Norml2();
            n_vec /= F;

            assert(1. - n_vec.Norml2() < 1e-12);

            // Get max wave speed
            d = dij_sparse->Elem(c, cp);

            // Compute intermediate face velocity
            Vf = pb->velocity(Uc);
            Vf += pb->velocity(Ucp);
            Vf *= 0.5;

            cout << "d: " << d
                 << ", tau cp: " << Ucp[0]
                 << ", tau c: " << Uc[0] << endl;

            double coeff = d * (Ucp[0] - Uc[0]) / F;
            Vf.Add(coeff, n_vec);
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

      // Put face velocity into object
      for (int i = 0; i < dim; i++)
      {
         int index = face + i * num_faces;
         v_CR_gf[index] = Vf[i];
      }

   } // End face iterator

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
* Function: ComputeMeshVelocities
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
*        2) Computes corrected node velocities on the corner vertices.
*        3) Computes the corrective face velocities.
*        4) Fills unnecessary node at cell center with the average velocity of the corner nodes.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeMeshVelocities(Vector &S, 
                                                      const double & t, 
                                                      double & dt,
                                                      const string flag, // Default NA
                                                      void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL)
{
   cout << "=======================================\n"
        << "         ComputeMeshVelocities         \n"
        << "=======================================\n";

   ComputeIntermediateFaceVelocities(S, t, flag, test_vel);  

   // Compute v_geo_gf
   ComputeGeoV();

   ComputeNodeVelocities(S, t, dt);

   // EnforceMVBoundaryConditions(S, t, dt);

   if (dim > 1)
   {
      // Don't need to fill face velocities with average in dim=1
      // FillFaceVelocitiesWithAvg(S);
      ComputeCorrectiveFaceVelocities(S, t, dt, flag, test_vel);

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
         val = 1.;
         k = pmesh->GetElementVolume(ci);
         temp_sum = 0.;
         pmesh->GetElementVertices(ci, fids);
         for (int j=0; j < fids.Size(); j++) // Face iterator
         {
            CalcOutwardNormalInt(S, ci, fids[j], n);
            GetIntermediateFaceVelocity(fids[j], Vf);
            temp_sum += n[0] * Vf[0];
         }
         val += (dt / k) * temp_sum;
         assert(val > 0.);
      }
   }

   FillCenterVelocitiesWithAvg(S);
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
   cout << "=======================================\n"
        << "             CalcMassLoss              \n"
        << "=======================================\n";
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
         
         // Fill corresponding cell to indicate graphically mass was broken
         mc_gf[ci] = 1.;
      }
      else
      {
         // Fill corresponding cell to indicate graphically mass was not broken
         mc_gf[ci] = 0.;
      }
   }
   double cell_ratio = (double)counter / (double)NDofs_L2;

   cout << "Percentage of cells where mass conservation was broken: " << cell_ratio << endl;
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

   double alpha1 = -1. * b + pm;
   alpha1 /= (2. * a);

   double alpha2 = -1. * b - pm;
   alpha2 /= (2. * a);

   alpha = std::max(alpha1, alpha2);

   if (alpha <= 0.)
   {
      cout << "d1: " << alpha1 << ", d2: " << alpha2 << endl;
      MFEM_ABORT("Alpha_i should be positive.\n");
   }

   cout << "alpha1: " << alpha1 << endl;
   cout << "alpha2: " << alpha2 << endl;
   cout << "alpha: " << alpha << endl;
}


/****************************************************************************************************
* Function: ComputeNodeVelocityRT
* Parameters:
*  node          - global index of node to compute velocity at
*  dt            - timestep
*  node_v        - computed velocity at node
*  is_dt_changed - Boolean representing if the timestep must be changed to accomodate for the
*                  values computed at this node.  If is_dt_changed is set to true, the loop 
*                  that computes the nodal velocities at the corner nodes must restart to 
*                  ensure the node velocities are in sync.
*
* Purpose:
*  This function computes the nodal velocity as explained in section 5 of [1].  This function 
*  also enforces the corresponding time step restriction imposed by the given node.
****************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::
   ComputeNodeVelocityRT(const int & node, double & dt, Vector &node_v, bool &is_dt_changed)
{
   cout << "=======================================\n"
        << "         ComputeNodeVelocityRT         \n"
        << "=======================================\n";
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
         ComputeCiGeo(node, Ci);

         // Enforce time restriction imposed by calculation of alpha_i
         if (dim == 2)
         {
            double trace = Ci.Trace();
            double det = Ci.Det();
            double val = 2. * sqrt(det);
            cout << "trace: " << trace << ", det: " << det << endl;

            // alpha_i must be positive restriction
            if (det <= 1.e-12 && trace < 0.)
            {
               if (dt > 2. / abs(trace))
               {
                  dt = 2. / (abs(trace) + 1.);
                  is_dt_changed = true;
                  cout << "timestep restriction from velocity computation, det = 0 case.\n";
               }
            }
            
            // alpha_i must be real
            else if (det > 0.)
            {
               if (trace < val)
               {
                  // Enforce timestep restriction
                  double zero2 = 1. / (-trace + val);
                  if (dt > 2 * zero2)
                  {
                     cout << "timestep restriction from velocity computation\n";
                     cout << "dt: " << dt << ", z2: " << zero2 << endl;
                     dt = 2. * zero2;
                     is_dt_changed = true;
                  }
               }
               else 
               {
                  cout << "positive determinant but no timestep restriction needed.\n";
               }
            }
            else
            {
               cout << "negative determinant, no timestep restriction needed.\n";
            }
         } // End time restriction from velocity computation

         GetViGeo(node, Vgeo);
         cout << "cnv_RT Vgeo: ";
         Vgeo.Print(cout);

         ComputeDeterminant(Ci, dt, d);

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

         // TODO: Raviart-Thomas implementation
         // mfem::Mesh::FaceInformation FI;
         // Array<int> row, cell_faces, vertex_edge_row;
         // int row_length, cell_faces_length, vertex_edge_row_length;

         // // Get adjacent edges
         // vertex_edge.GetRow(node, vertex_edge_row);
         // vertex_edge_row_length = vertex_edge_row.Size();

         // Vector n_vec(dim);
         // double mi = 0., Ti = 0., Ki=0.;
         // Vector Uc(dim+2), Vi(dim), ivf(dim), v_cell(dim);
         // int vertex_bdr_attribute = 0;

         // node_v = 0.; // Reset vertex velocity 

         // // Get adjacent elements
         // vertex_element->GetRow(node, row);
         // row_length = row.Size();

         // cout << "node: " << node << endl;
         // // Iterate over cells that share this vertex and average their contribution to vertex velocity
         // for (int element_it = 0; element_it < row_length; element_it++)
         // {
         //    v_cell = 0.;
         //    int el_index = row[element_it];
         //    cout << "el: " << el_index << endl;

         //    // Get element faces
         //    element_face.GetRow(el_index, cell_faces);
         //    cell_faces_length = cell_faces.Size();

         //    // Iterate over faces of this cell that touch vertex
         //    for (int face_it = 0; face_it < cell_faces_length; face_it++)
         //    {
         //       ivf = 0.;
         //       int face_index = cell_faces[face_it];
         //       cout << "face: " << face_index << endl;

         //       if (vertex_edge_row.Find(face_index) != -1)
         //       {
         //          cout << "element " << el_index << " has adjacent face: " << face_index << endl;
         //          GetIntermediateFaceVelocity(face_index, ivf);
         //          cout << "ifv: ";
         //          ivf.Print(cout);
         //          v_cell.Add(1., ivf);
         //       }
         //       else
         //       {
         //          cout << "face " << face_index << " is not adjacent to node " << node << endl;
         //       }
         //    }

         //    node_v.Add(1., v_cell);
         // }
         // // Average contributions
         // cout << "node_v pre division: ";
         // node_v.Print(cout);
         // cout << "row_length: " << row_length << endl;
         // node_v /= row_length;
         // cout << "node_v: ";
         // node_v.Print(cout);
      }
   }
}

/***********************************************************************************************************
Function: IntGradRT
Parameters:
   cell - index corrseponding to the cell (K_c)
   res  - (dim x dim) DenseMatrix representing the outer product of Vfm with 
          the gradient of its corresponding scalar RT shape function.
Purpose:
   This function calculates the integral of the gradient of the RT velocity function
   on a given cell.

   NOTE: This function assumes that dim > 1.
   NOTE: This function assumes that the function LagrangianLOOperator:ComputeIntermediateFaceVelocities()
   has already been called.  If this function has not been called, then the returned velocity will be 0.
***********************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::IntGradRT(const int cell, DenseMatrix & res)
{
   // cout << "=======================================\n"
   //      << "               IntGradRT               \n"
   //      << "=======================================\n";
   ParGridFunction CRc_gf(&CRc);
   const int _size = CRc.GetVSize();

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
         CRc_gf.MakeRef(&CRc, v_CR_gf, j*_size);
         CRc_gf.GetGradient(*trans, grad);

         // Put information into Dense Matrix
         res.GetRow(j, row);
         row.Add(ip.weight, grad);
         res.SetRow(j, row);
      }
   }
}

/***********************************************************************************************************
Function: ComputeCiGeo
Parameters:
   node - Index corresponding to the global node index, to include face nodes and cell center nodes.  
          Note that we should not be computing any sort of averaged velocity at the cell centers.
   res  - Averaged DenseMatrix resulting from the average gradient of the cell wise velocities of all
          adjacent cells.
Purpose:
   Part of the process to reconstructing the continuous geometric velocity field from the discontinuous
   RT representation.  The exact equation to be solved is given by equation (5.11)

   NOTE: This function assumes dim > 1.
   NOTE: This function assumes that the function LagrangianLOOperator:ComputeIntermediateFaceVelocities()
   has already been called.  If this function has not been called, then the returned velocity will be 0.
***********************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeCiGeo(const int &node, DenseMatrix & res)
{
   // cout << "=======================================\n"
   //      << "             ComputeCiGeo              \n"
   //      << "=======================================\n";
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
         cout << "ComputeCiGeo::corner node\n";
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
      IntGradRT(row_el, dm_temp);
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
Function: ComputeGeoV
Purpose:
   Reconstruct and average RT velocity field.  This projects the discontinuous GF from the 
   Crouzeix-Raviart FE space onto the H1 space, filling the valies into v_geo_gf. (5.11)

   NOTE: This function assumes that the function LagrangianLOOperator:ComputeIntermediateFaceVelocities()
   has already been called.  If this function has not been called, then the returned velocity will be 0.
***********************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::ComputeGeoV()
{
   cout << "=======================================\n"
        << "             ComputeViGeo              \n"
        << "=======================================\n";
   // cout << "v_CR_gf:\n";
   // v_CR_gf.Print(cout);

   VectorGridFunctionCoefficient vcr_coeff(&v_CR_gf);
   v_geo_gf.ProjectDiscCoefficient(vcr_coeff, mfem::ParGridFunction::AvgType::ARITHMETIC);  
   
   // cout << "v_geo gf:\n";
   // v_geo_gf.Print(cout);
}


/***********************************************************************************************************
* Function: ComputeNodeVelocities
* Parameters:
*   S        - BlockVector representing FiniteElement information
*   t        - Current time
*   dt       - Timestep
*   flag     - Flag used to indicate testing
*   test_vel - Velocity used for testing
*
* Purpose: 
*  This function computes the velocities on the corner nodes of the mesh via the process outlined in
*  section 5 of the Lagrangian paper.
***********************************************************************************************************/
template<int dim>
void LagrangianLOOperator<dim>::
   ComputeNodeVelocities(Vector &S, 
                           const double & t, 
                           double & dt, 
                           const string flag, // Default NA
                           void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
   cout << "=======================================\n"
        << "         ComputeNodeVelocities         \n"
        << "=======================================\n";
   Vector vertex_v(dim);

   bool is_dt_changed = false;
   // Iterate over vertices
   for (int vertex = 0; vertex < num_vertices; vertex++) // Vertex iterator
   {
      ComputeNodeVelocityRT(vertex, dt, vertex_v, is_dt_changed);

      if (vertex_v[0] != vertex_v[0] || vertex_v[1] != vertex_v[1])
      {
         cout << "NaN velocity encountered in ComputeNodeVelocities at vertex: " << vertex << endl;
         MFEM_ABORT("Aborting due to NaNs.\n");
      }

      UpdateNodeVelocity(S, vertex, vertex_v);

      // If we restricted the timestep, we must recompute the vertex velocities that were computed previously
      if (is_dt_changed)
      {
         vertex = -1;
         is_dt_changed = false;
         cout << "Restarting vertex iterator\n";
      }
   } // End Vertex iterator

   // Enforce bounndary conditions on mesh velocity
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
   cout << "==============================================\n";
   cout << "Computing corrective interior face velocities.\n";
   cout << "==============================================\n";
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
      GetNodePosition(S, face_dof, face_x);
      GetNodePosition(S, face_vdof1, vdof1_x);
      GetNodePosition(S, face_vdof2, vdof2_x);

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
      GetNodeVelocity(S, face_vdof1, vdof1_v);
      GetNodeVelocity(S, face_vdof2, vdof2_v);
      
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
         
         GetIntermediateFaceVelocity(face, Vf);
         
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
         cout << "Houston we have a problem.\n";
         MFEM_ABORT("Exit");
      }
      UpdateNodeVelocity(S, face_dof, face_velocity);
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
      // if (vel[i] != vel[i]) { 
      //    cout << "Node: " << node << endl;
      //    MFEM_ABORT("NaN val encountered."); 
      // }
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
* Function: get_vi_deo
* Parameters:
*   node - Global index of node in question.
*   vel  - Velocity Vi^geo at given node
*
* Purpose:
*  This function returns Vi_geo.
*  NOTE: Function LagrangianLOOperator::ComputeGeoV must be called first.
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
      // cout << "i: " << i << ", val: " << x_gf[index] << endl;
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

/* Explicit n_vectantiation */
template class LagrangianLOOperator<1>;
template class LagrangianLOOperator<2>;
template class LagrangianLOOperator<3>;

} // end ns hydrodynamics

} // end ns mfem
