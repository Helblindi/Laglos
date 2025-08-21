#include "mfem.hpp"
#include "laglos_solver.hpp"
#include "test_problems_include.h"
#include "var-config.h"
#include <cassert>

using namespace std;
using namespace mfem;
using namespace hydroLO;

/* Parameters */
int dim = 2;
int order_mv = 1;          // Order of mesh movement approximation space
int order_t = 0;
string result = string(LAGLOS_DIR) + "data/ref-square.mesh";
const char* mesh_file_location = result.c_str();

bool use_viscosity = true; // Doesn't matter
const bool _mm = true;          // Doesn't matter
double CFL = 0.5;          // Doesn't matter
int rs_levels = 2;
double t_init = 0.0;
const int elastic_eos = 0;

static int num_procs = 0;
static int myid = -1;
/* ---------------- End Parameters ---------------- */

int CompareCijComputation();
int CheckCijSkewSymmetry();
int TestDGNormalIntegrator();
int ValidateOT2();

int main(int argc, char *argv[])
{
   // Initialize MPI.
   Mpi::Init();
   num_procs = Mpi::WorldSize();
   myid = Mpi::WorldRank();

   int d = 0;
   d += CompareCijComputation();
   d += CheckCijSkewSymmetry();
   d += TestDGNormalIntegrator();
   d += ValidateOT2();

   return d;
}


/**
 * @brief Test function to validate the correctness of the Cij matrix construction in the LagrangianLOOperator.
 *
 * This function sets up a parallel mesh and finite element spaces, initializes problem-specific coefficients,
 * and constructs the Lagrangian low-order operator. It then calls BuildCijMatrices() to compute the Cij matrices,
 * which are essential for the Lagrangian hydrodynamics solver. The function iterates over all interior faces
 * of the mesh, computes the outward normal and the corresponding Cij vector for each face, and checks that
 * the computed Cij matches the expected outward normal (up to a tolerance). If any mismatch is detected,
 * the function prints diagnostic information and returns 1 to indicate failure. Otherwise, it returns 0
 * to indicate success.
 *
 * This test is intended to be called from main() to automatically verify the correctness of Cij matrix
 * construction after any code changes.
 *
 * @return int Returns 0 if all Cij computations match the expected normals, 1 otherwise.
 */
int CompareCijComputation()
{
   cout << "========================================================================\n"
        << "\tCompare Cij computation to outward normal for all interior faces (2D).\n";
   // Create mesh
   Mesh *mesh = new Mesh(mesh_file_location, true, true);
   dim = mesh->Dimension();

   // Refine the mesh in serial to increase the resolution. In this example
   // we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
   // a command-line parameter. If the mesh is of NURBS type, we convert it
   // to a (piecewise-polynomial) high-order mesh.
   for (int lev = 0; lev < rs_levels; lev++)
   {
      mesh->UniformRefinement();
   }

   // Define the parallel mesh by a partitioning of the serial mesh. Refine
   // this mesh further in parallel to increase the resolution. Once the
   // parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;

   int NE = pmesh->GetNE(), ne_min, ne_max;
   MPI_Reduce(&NE, &ne_min, 1, MPI_INT, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&NE, &ne_max, 1, MPI_INT, MPI_MAX, 0, pmesh->GetComm());
   double hmin, hmax, kmin, kmax;
   pmesh->GetCharacteristics(hmin, hmax, kmin, kmax);

   // Set up problem
   ProblemBase * problem_class = new SodProblem(dim);

   // Define the parallel finite element spaces. We use:
   // - H1 (Q2, continuous) for mesh movement.
   // - L2 (Q0, discontinuous) for state variables
   // - CR/RT for mesh reconstruction at nodes
   H1_FECollection H1FEC(order_mv, dim);
   H1_FECollection H1FEC_L(1, dim);
   L2_FECollection L2FEC(order_t, dim, BasisType::GaussLobatto);
   FiniteElementCollection * CRFEC;
   if (dim == 1)
   {
      CRFEC = new CrouzeixRaviartFECollection();
   }
   else
   {
      CRFEC = new RT_FECollection(0, dim);
   }

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace H1FESpace_L(pmesh, &H1FEC_L, dim);
   /* Finite element space solely constructed for continuous representation of density field */
   ParFiniteElementSpace H1cFESpace(pmesh, &H1FEC, 1);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);
   ParFiniteElementSpace CRFESpace(pmesh, CRFEC, dim);

   /* The monolithic BlockVector stores unknown fields as:
   *   - 0 -> position
   *   - 1 -> specific volume
   *   - 2 -> velocity (L2V)
   *   - 3 -> speific total energy
   */
   const int Vsize_l2 = L2FESpace.GetVSize();
   const int Vsize_l2v = L2VFESpace.GetVSize();
   const int Vsize_h1 = H1FESpace.GetVSize();
   Array<int> offset(5);
   offset[0] = 0;
   offset[1] = offset[0] + Vsize_h1;
   offset[2] = offset[1] + Vsize_l2;
   offset[3] = offset[2] + Vsize_l2v;
   offset[4] = offset[3] + Vsize_l2;
   BlockVector S(offset, Device::GetMemoryType());

   // Define GridFunction objects for the position, mesh velocity and specific
   // volume, velocity, and specific internal energy. At each step, each of 
   // these values will be updated.
   ParGridFunction x_gf, sv_gf, v_gf, ste_gf;
   ParGridFunction rho0_gf(&L2FESpace);
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   sv_gf.MakeRef(&L2FESpace, S, offset[1]);
   v_gf.MakeRef(&L2VFESpace, S, offset[2]);
   ste_gf.MakeRef(&L2FESpace, S, offset[3]);

   // Initialize x_gf using starting mesh positions
   pmesh->SetNodalGridFunction(&x_gf);
   // Sync the data location of x_gf with its base, S
   x_gf.SyncAliasMemory(S);

   // Change class variables into static std::functions since virtual static member functions are not an option
   // and Coefficient class requires std::function arguments
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> sv0_static = 
      std::bind(&ProblemBase::sv0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<void(const Vector &, const double, Vector &)> v0_static = 
      std::bind(&ProblemBase::v0, problem_class, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
   std::function<double(const Vector &,const double)> ste0_static = 
      std::bind(&ProblemBase::ste0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> rho0_static = 
      std::bind(&ProblemBase::rho0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> p0_static = 
      std::bind(&ProblemBase::p0, problem_class, std::placeholders::_1, std::placeholders::_2);

   // Initialize specific volume, velocity, and specific total energy
   FunctionCoefficient sv_coeff(sv0_static);
   sv_coeff.SetTime(t_init);
   sv_gf.ProjectCoefficient(sv_coeff);
   sv_gf.SyncAliasMemory(S);

   VectorFunctionCoefficient v_coeff(dim, v0_static);
   v_coeff.SetTime(t_init);
   v_gf.ProjectCoefficient(v_coeff);
   v_gf.SyncAliasMemory(S);

   // While the ste_coeff is not used for initialization in the Sedov case,
   // it is necessary for plotting the exact solution
   FunctionCoefficient ste_coeff(ste0_static);
   ste_coeff.SetTime(t_init);
   ste_gf.ProjectCoefficient(ste_coeff);
   ste_gf.SyncAliasMemory(S);

   // PLF to build mass vector
   FunctionCoefficient rho_coeff(rho0_static); 
   rho_coeff.SetTime(t_init);
   IntegrationRule ir = IntRules.Get(pmesh->GetElementBaseGeometry(0), 3 * H1FESpace.GetOrder(0) + L2FESpace.GetOrder(0) - 1);
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff, &ir));
   m->Assemble();

   /* Initialize rho0_gf */
   rho0_gf.ProjectCoefficient(rho_coeff);

   FunctionCoefficient p_coeff(p0_static);
   p_coeff.SetTime(t_init);

   /* Create Lagrangian Low Order Solver Object */
   LagrangianLOOperator hydro(dim, S.Size(), H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, rho_coeff, rho0_gf, m, ir, problem_class, offset, use_viscosity, elastic_eos, _mm, CFL);
   hydro.BuildCijMatrices();
   // assert(false);
   Array<int> fids, oris;
   mfem::Mesh::FaceInformation FI;
   Vector n_int(dim), c(dim), cij(dim), _t(dim);
   int cj = 0;
   for (int ci = 0; ci < NE; ci++)
   {
      pmesh->GetElementEdges(ci, fids, oris);
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
               assert(ci == FI.element[1].index);
               cj = FI.element[0].index; 
            }

            hydro.CalcOutwardNormalInt(S, ci, fids[j], n_int);
            c = n_int;
            c /= 2.;

            hydro.GetLocalCij(ci, cj, cij);
            
            subtract(c, cij, _t);
            if (_t.Norml2() > 1.e-12)
            {
               cout << "cij computation does not match outward normal!\n";
               cout << "i: " << ci << ", j: " << cj << endl;
               cout << "Outward normal: ";
               c.Print(cout);
               cout << "cij: ";
               cij.Print(cout);
               return 1;
            }

         }
      }
   }

   /* free memory */
   delete problem_class;
   delete CRFEC;
   delete m;
   delete pmesh;

   cout << "Passed Cij comparison to outward normal.\n"
        << "========================================================================\n\n";

   return 0;
}

/*
The values of Cij should match the definition given in equation (4.17) in the paper:
Guermond, Popov, Tomas: 
"Invariant domain preserving discretization-independent schemes and convex limiting for hyperbolic systems."
*/
int CheckCijSkewSymmetry()
{
   // This function is a placeholder for validating the Cij computation order 1.
   // It can be implemented with specific tests or assertions as needed.
   cout << "========================================================================\n"
        << "\tValidating Cij skew symmetry..." << endl;

   Mesh *mesh = new Mesh(Mesh::MakeCartesian2D(1, 1, Element::QUADRILATERAL, true));
   dim = mesh->Dimension();
   /* Optional serial refinement */
   int num_refinements = 1;
   for (int lev = 0; lev < num_refinements; lev++)
   {
      mesh->UniformRefinement();
   }
   
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   ParMesh *smesh = new ParMesh(ParMesh::MakeRefined(*pmesh, order_t+1, BasisType::ClosedUniform));;
   delete mesh;

   // Set up problem
   ProblemBase * problem_class = new SodProblem(dim);

   // Define the parallel finite element spaces. We use:
   // - H1 (Q2, continuous) for mesh movement.
   // - L2 (Q0, discontinuous) for state variables
   // - CR/RT for mesh reconstruction at nodes
   H1_FECollection H1FEC(order_mv, dim);
   H1_FECollection H1FEC_L(1, dim);
   L2_FECollection L2FEC(order_t, dim, BasisType::GaussLobatto);
   FiniteElementCollection * CRFEC;
   if (dim == 1)
   {
      CRFEC = new CrouzeixRaviartFECollection();
   }
   else
   {
      CRFEC = new RT_FECollection(0, dim);
   }

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace H1FESpace_L(pmesh, &H1FEC_L, dim);
   /* Finite element space solely constructed for continuous representation of density field */
   ParFiniteElementSpace H1cFESpace(pmesh, &H1FEC, 1);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);
   ParFiniteElementSpace CRFESpace(smesh, CRFEC, dim);
   const int NDofs_L2 = L2FESpace.GetNDofs();

   /* The monolithic BlockVector stores unknown fields as:
   *   - 0 -> position
   *   - 1 -> specific volume
   *   - 2 -> velocity (L2V)
   *   - 3 -> speific total energy
   */
   const int Vsize_l2 = L2FESpace.GetVSize();
   const int Vsize_l2v = L2VFESpace.GetVSize();
   const int Vsize_h1 = H1FESpace.GetVSize();
   Array<int> offset(5);
   offset[0] = 0;
   offset[1] = offset[0] + Vsize_h1;
   offset[2] = offset[1] + Vsize_l2;
   offset[3] = offset[2] + Vsize_l2v;
   offset[4] = offset[3] + Vsize_l2;
   BlockVector S(offset, Device::GetMemoryType());

   // Define GridFunction objects for the position, mesh velocity and specific
   // volume, velocity, and specific internal energy. At each step, each of 
   // these values will be updated.
   ParGridFunction x_gf, sv_gf, v_gf, ste_gf;
   ParGridFunction rho0_gf(&L2FESpace);
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   sv_gf.MakeRef(&L2FESpace, S, offset[1]);
   v_gf.MakeRef(&L2VFESpace, S, offset[2]);
   ste_gf.MakeRef(&L2FESpace, S, offset[3]);

   // Initialize x_gf using starting mesh positions
   pmesh->SetNodalGridFunction(&x_gf);
   // Sync the data location of x_gf with its base, S
   x_gf.SyncAliasMemory(S);

   // Change class variables into static std::functions since virtual static member functions are not an option
   // and Coefficient class requires std::function arguments
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> sv0_static = 
      std::bind(&ProblemBase::sv0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<void(const Vector &, const double, Vector &)> v0_static = 
      std::bind(&ProblemBase::v0, problem_class, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
   std::function<double(const Vector &,const double)> ste0_static = 
      std::bind(&ProblemBase::ste0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> rho0_static = 
      std::bind(&ProblemBase::rho0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> p0_static = 
      std::bind(&ProblemBase::p0, problem_class, std::placeholders::_1, std::placeholders::_2);

   // Initialize specific volume, velocity, and specific total energy
   FunctionCoefficient sv_coeff(sv0_static);
   sv_coeff.SetTime(t_init);
   sv_gf.ProjectCoefficient(sv_coeff);

   VectorFunctionCoefficient v_coeff(dim, v0_static);
   v_coeff.SetTime(t_init);
   v_gf.ProjectCoefficient(v_coeff);

   // While the ste_coeff is not used for initialization in the Sedov case,
   // it is necessary for plotting the exact solution
   FunctionCoefficient ste_coeff(ste0_static);
   ste_coeff.SetTime(t_init);
   ste_gf.ProjectCoefficient(ste_coeff);


   // PLF to build mass vector
   FunctionCoefficient rho_coeff(rho0_static); 
   rho_coeff.SetTime(t_init);
   int ir_order = 3 * H1FESpace.GetOrder(0) + L2FESpace.GetOrder(0) - 1;
   // ir_order = 0;
   IntegrationRules _IntRules(0, Quadrature1D::GaussLobatto);
   IntegrationRule ir = _IntRules.Get(pmesh->GetElementBaseGeometry(0), ir_order);
   // for (int p = 0; p < ir.GetNPoints(); p++)
   // {
   //    const IntegrationPoint &ip = ir.IntPoint(p);
   //    cout << "ip.x: " << ip.x << ", ip.y: " << ip.y << "\n";
   // }
   // assert(false);
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff, &ir));
   m->Assemble();

   /* Initialize rho0_gf */
   rho0_gf.ProjectCoefficient(rho_coeff);

   FunctionCoefficient p_coeff(p0_static);
   p_coeff.SetTime(t_init);

   /* Create Lagrangian Low Order Solver Object */
   LagrangianLOOperator hydro(dim, S.Size(), H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, rho_coeff, rho0_gf, m, ir, problem_class, offset, use_viscosity, elastic_eos, _mm, CFL);
   hydro.BuildCijMatrices();
   Table L2Connectivity;
   hydro.GetL2ConnectivityTable(L2Connectivity);
   delete pmesh;
   bool _is_ss = hydro.IsSkewSymmetric();
   cout << "Done\n"
        << "========================================================================\n\n";
   
   return _is_ss ? 0 : 1; // Return 0 if skew-symmetric, 1 otherwise
}

int TestDGNormalIntegrator()
{
   // This function is a placeholder for testing the DGNormalIntegrator.
   // It can be implemented with specific tests or assertions as needed.
   cout << "Testing DGNormalIntegrator..." << endl;

   Mesh *mesh = new Mesh(Mesh::MakeCartesian2D(1, 1, Element::QUADRILATERAL, true));
   dim = mesh->Dimension();
   /* Optional serial refinement */
   int num_refinements = 0;
   for (int lev = 0; lev < num_refinements; lev++)
   {
      mesh->UniformRefinement();
   }

   L2_FECollection L2FEC(1, dim, BasisType::Positive);
   FiniteElementSpace L2(mesh, &L2FEC);
   const int NDofs_L2 = L2.GetNDofs();
   BilinearForm xbf(&L2), ybf(&L2);
   ConstantCoefficient one(1.0);

   // xbf.AddDomainIntegrator(new DerivativeIntegrator(one, 0));
   xbf.AddInteriorFaceIntegrator(new DGNormalIntegrator(-1., 0));
   xbf.AddBdrFaceIntegrator(new DGNormalIntegrator(-1., 0));
   xbf.Assemble();
   xbf.Finalize();
   SparseMatrix xmat_sparse = xbf.SpMat();

   // ybf.AddDomainIntegrator(new DerivativeIntegrator(one, 1));
   ybf.AddInteriorFaceIntegrator(new DGNormalIntegrator(-1., 1));
   ybf.AddBdrFaceIntegrator(new DGNormalIntegrator(-1., 1));
   ybf.Assemble();
   ybf.Finalize();
   SparseMatrix ymat_sparse = ybf.SpMat();
   
   // cout << "xmat_sparse:\n";
   // xmat_sparse.Print(cout);
   // cout << "ymat_sparse:\n";
   // ymat_sparse.Print(cout);

   // Array<int> cols;
   // Vector vals, cij(dim);;
   // int col_index;

   /* Get cij */
   // for (int i = 0; i < NDofs_L2; i++)
   // {
   //    for (int j = 0; j < NDofs_L2; j++)
   //    {
   //       cij = 0.;
   //       xmat_sparse.GetRow(i, cols, vals);
   //       col_index = cols.Find(j);
   //       if (col_index != -1) {
   //          cij[0] = vals[col_index];
   //       } else {
   //          // cout << "col not found\n";
   //          cij[0] = 0.;
   //       }

   //       ymat_sparse.GetRow(i, cols, vals);
   //       col_index = cols.Find(j);
   //       if (col_index != -1) {
   //          cij[1] = vals[col_index];
   //       } else {
   //          // cout << "col not found for i: " << i << ", j: " << j << "\n";
   //          cij[1] = 0.;
   //       }
   //       cout << "Cij[" << i << "][" << j << "] = ";
   //       cij.Print(cout);
   //    }
   // }

   /* Free memory */
   delete mesh;

   /* Check symmetry */
   if (xmat_sparse.IsSymmetric() < 1.E-10 && ymat_sparse.IsSymmetric() < 1.E-10)
   {
      cout << "Both matrices are symmetric.\n";
      return 0;
   }
   else
   {
      cout << "At least one matrix is not symmetric.\n";
      cout << "symm val x: " << xmat_sparse.IsSymmetric() << endl;
      cout << "symm val y: " << ymat_sparse.IsSymmetric() << endl;
      cout << "xmat_sparse:\n";
      xmat_sparse.Print(cout);
      cout << "ymat_sparse:\n";
      ymat_sparse.Print(cout);
      return 1;
   }

   return 0;
}

// In 1D, the basis functions are:
//   - f0(x) = 2x^2 - 3x + 1
//   - f1(x) = -4(x^2 - x)
//   - f2(x) = 2x^2 - x
// For higher dimensions, can take the tensor product of the 1D basis functions.
// I.e., in 2D the basis function for the following node is: 
//
//       p0(x,y) = (2x^2 -3x +1)(2y^2 -3y +1)
//   6--7--8
//   |  |  |
//   3--4--5
//   |  |  |
//   0--1--2
//
// This test will verify we compute c28 correctly. Note:
//    p2(x,y) = (2x^2 - x)(2y^2 - 3y + 1)
//    p8(x,y) = (2x^2 - 3x + 1)(2y^2 - 3y + 1)
// c28 = c28K - c28B
//    where
//   c28K = int_K p2(x,y) grad p8(x,y) dx dy
//   c28B = (1/2) int_DK p2(x,y) p8(x,y) nK dx dy
//    where nK = [1,0]^T.
// The correct value of c28 is:
//   c28 = [-1/120 -1/90]^T - [-1/120 0]^T = [0 -1/90]^T 
int ValidateOT2()
{
   cout << "========================================================================\n"
        << "\tValidate -ot 2.\n";
   // Create mesh
   Mesh *mesh = new Mesh(mesh_file_location, true, true);
   dim = mesh->Dimension();

   // Refine the mesh in serial to increase the resolution. In this example
   // we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
   // a command-line parameter. If the mesh is of NURBS type, we convert it
   // to a (piecewise-polynomial) high-order mesh.
   for (int lev = 0; lev < 1; lev++)
   {
      mesh->UniformRefinement();
   }

   // Define the parallel mesh by a partitioning of the serial mesh. Refine
   // this mesh further in parallel to increase the resolution. Once the
   // parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   ParMesh *smesh = new ParMesh(ParMesh::MakeRefined(*pmesh, 3, BasisType::ClosedUniform));;
   delete mesh;

   int NE = pmesh->GetNE(), ne_min, ne_max;
   MPI_Reduce(&NE, &ne_min, 1, MPI_INT, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&NE, &ne_max, 1, MPI_INT, MPI_MAX, 0, pmesh->GetComm());
   double hmin, hmax, kmin, kmax;
   pmesh->GetCharacteristics(hmin, hmax, kmin, kmax);

   // Set up problem
   ProblemBase * problem_class = new SodProblem(dim);

   // Define the parallel finite element spaces. We use:
   // - H1 (Q2, continuous) for mesh movement.
   // - L2 (Q0, discontinuous) for state variables
   // - CR/RT for mesh reconstruction at nodes
   H1_FECollection H1FEC(3, dim);
   H1_FECollection H1FEC_L(1, dim);
   L2_FECollection L2FEC(2, dim, BasisType::GaussLobatto);
   FiniteElementCollection * CRFEC;
   if (dim == 1)
   {
      CRFEC = new CrouzeixRaviartFECollection();
   }
   else
   {
      CRFEC = new RT_FECollection(0, dim);
   }

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace H1FESpace_L(pmesh, &H1FEC_L, dim);
   /* Finite element space solely constructed for continuous representation of density field */
   ParFiniteElementSpace H1cFESpace(pmesh, &H1FEC, 1);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);
   ParFiniteElementSpace CRFESpace(pmesh, CRFEC, dim);

   /* The monolithic BlockVector stores unknown fields as:
   *   - 0 -> position
   *   - 1 -> specific volume
   *   - 2 -> velocity (L2V)
   *   - 3 -> speific total energy
   */
   const int Vsize_l2 = L2FESpace.GetVSize();
   const int Vsize_l2v = L2VFESpace.GetVSize();
   const int Vsize_h1 = H1FESpace.GetVSize();
   Array<int> offset(5);
   offset[0] = 0;
   offset[1] = offset[0] + Vsize_h1;
   offset[2] = offset[1] + Vsize_l2;
   offset[3] = offset[2] + Vsize_l2v;
   offset[4] = offset[3] + Vsize_l2;
   BlockVector S(offset, Device::GetMemoryType());

   // Define GridFunction objects for the position, mesh velocity and specific
   // volume, velocity, and specific internal energy. At each step, each of 
   // these values will be updated.
   ParGridFunction x_gf, sv_gf, v_gf, ste_gf;
   ParGridFunction rho0_gf(&L2FESpace);
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   sv_gf.MakeRef(&L2FESpace, S, offset[1]);
   v_gf.MakeRef(&L2VFESpace, S, offset[2]);
   ste_gf.MakeRef(&L2FESpace, S, offset[3]);

   // Initialize x_gf using starting mesh positions
   pmesh->SetNodalGridFunction(&x_gf);
   // Sync the data location of x_gf with its base, S
   x_gf.SyncAliasMemory(S);

   // Change class variables into static std::functions since virtual static member functions are not an option
   // and Coefficient class requires std::function arguments
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> sv0_static = 
      std::bind(&ProblemBase::sv0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<void(const Vector &, const double, Vector &)> v0_static = 
      std::bind(&ProblemBase::v0, problem_class, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
   std::function<double(const Vector &,const double)> ste0_static = 
      std::bind(&ProblemBase::ste0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> rho0_static = 
      std::bind(&ProblemBase::rho0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> p0_static = 
      std::bind(&ProblemBase::p0, problem_class, std::placeholders::_1, std::placeholders::_2);

   // Initialize specific volume, velocity, and specific total energy
   FunctionCoefficient sv_coeff(sv0_static);
   sv_coeff.SetTime(t_init);
   sv_gf.ProjectCoefficient(sv_coeff);
   sv_gf.SyncAliasMemory(S);

   VectorFunctionCoefficient v_coeff(dim, v0_static);
   v_coeff.SetTime(t_init);
   v_gf.ProjectCoefficient(v_coeff);
   v_gf.SyncAliasMemory(S);

   // While the ste_coeff is not used for initialization in the Sedov case,
   // it is necessary for plotting the exact solution
   FunctionCoefficient ste_coeff(ste0_static);
   ste_coeff.SetTime(t_init);
   ste_gf.ProjectCoefficient(ste_coeff);
   ste_gf.SyncAliasMemory(S);

   // PLF to build mass vector
   FunctionCoefficient rho_coeff(rho0_static); 
   rho_coeff.SetTime(t_init);
   IntegrationRule ir = IntRules.Get(pmesh->GetElementBaseGeometry(0), 3 * H1FESpace.GetOrder(0) + L2FESpace.GetOrder(0) - 1);
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff, &ir));
   m->Assemble();

   /* Initialize rho0_gf */
   rho0_gf.ProjectCoefficient(rho_coeff);

   FunctionCoefficient p_coeff(p0_static);
   p_coeff.SetTime(t_init);

   /* Create Lagrangian Low Order Solver Object */
   LagrangianLOOperator hydro(dim, S.Size(), H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, rho_coeff, rho0_gf, m, ir, problem_class, offset, use_viscosity, elastic_eos, _mm, CFL);
   hydro.BuildCijMatrices();

   Vector c28(dim), c28_exact(dim), _t(dim);
   hydro.GetLocalCij(2, 8, c28);
   c28_exact[0] = 0., c28_exact[1] = -1./90.;
   subtract(c28, c28_exact, _t);
   if (_t.Norml2() > 1.e-12)
   {
      cout << "c28 computation does not match expected value!\n";
      cout << "Computed c28: ";
      c28.Print(cout);
      cout << "Expected c28: ";
      c28_exact.Print(cout);
      return 1;
   }

   cout << "Done\n"
        << "========================================================================\n\n";

   /* free memory */
   delete problem_class;
   delete CRFEC;
   delete m;
   delete pmesh;
   delete smesh;
   
   return 0;
}