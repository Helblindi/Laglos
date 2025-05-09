#include "mfem.hpp"
#include "laglos_solver.hpp"
#include "mfem/fem/coefficient.hpp"
#include "test_problems_include.h"
#include "var-config.h"
#include <cassert>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;
using namespace mfem;
using namespace hydroLO;

/* ---------------- Parameters to be used for tests ---------------- */
// Problem
const int dim = 2;
const int problem = 1;     // This value doesn't matter (1-8)
int order_mv = 2;          // Order of mesh movement approximation space
int order_u = 0;
double tol = 1e-12;
const char *mesh_file_location = "/data/ref-square.mesh";
std::string result = std::string(LAGLOS_DIR) + std::string(mesh_file_location);
const char* mesh_file = result.c_str();
int mesh_refinements = 1;
bool use_viscosity = true; // Doesn't matter
bool _mm = true;          // Doesn't matter
int mv_option = 2;
int fv_option = 2;
double CFL = 0.5;          // Doesn't matter
static int num_procs = 0;
static int myid = -1;
bool suppress_test_output = false;
/* ---------------- End Parameters ---------------- */

void velocity_exact(const Vector &x, const double & t, Vector &v);
int tester();
int test_elasticity();

int main(int argc, char *argv[])
{
   // Initialize MPI.
   Mpi::Init();
   num_procs = Mpi::WorldSize();
   myid = Mpi::WorldRank();

   // Set precision for all output values
   cout.precision(12);

   // Parse command line arguments
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&mesh_refinements, "-r", "--refinement", "Number of mesh refinements.");

   args.Parse();

   int d = 0;
   d += tester();
   d += test_elasticity();

   return d;
}

/*
Purpose: provide a place for quick tests of Laglos
*/
int tester()
{
   cout << "=== Running test space function ===\n";

   // Initialize mesh
   mfem::Mesh *mesh;
   mesh = new mfem::Mesh(mesh_file, true, true);

   // Refine the mesh
   for (int lev = 0; lev < mesh_refinements; lev++)
   {
      mesh->UniformRefinement();
   }

   // Construct parmesh object
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);            
   delete mesh;

   ProblemBase * problem_class = new ElasticShocktube();

   // Template FE stuff to construct hydro operator
   H1_FECollection H1FEC(order_mv, dim);
   H1_FECollection H1FEC_L(1, dim);
   L2_FECollection L2FEC(order_u, dim, BasisType::Positive);
   CrouzeixRaviartFECollection CRFEC;

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace H1FESpace_L(pmesh, &H1FEC_L, dim);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);
   ParFiniteElementSpace CRFESpace(pmesh, &CRFEC, dim);

   // Output information
   pmesh->ExchangeFaceNbrData();

   // Construct blockvector
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

   // Pair with corresponding gridfunctions
   // Only need position and velocity for testing
   ParGridFunction x_gf, sv_gf, v_gf, ste_gf;
   ParGridFunction rho0_gf(&L2FESpace);
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   sv_gf.MakeRef(&L2FESpace, S, offset[1]);
   v_gf.MakeRef(&L2VFESpace, S, offset[2]);
   ste_gf.MakeRef(&L2FESpace, S, offset[3]);

   pmesh->SetNodalGridFunction(&x_gf);
   x_gf.SyncAliasMemory(S);
   const ParGridFunction x0_gf = x_gf;

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
   sv_coeff.SetTime(0.);
   sv_gf.ProjectCoefficient(sv_coeff);

   VectorFunctionCoefficient v_coeff(dim, v0_static);
   v_coeff.SetTime(0.);
   v_gf.ProjectCoefficient(v_coeff);

   FunctionCoefficient ste_coeff(ste0_static);
   ste_coeff.SetTime(0.);
   ste_gf.ProjectCoefficient(ste_coeff);

   // Just leave templated for hydro construction
   FunctionCoefficient rho_coeff(rho0_static);
   rho_coeff.SetTime(0.);
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff));
   m->Assemble();
   rho0_gf.ProjectCoefficient(rho_coeff);

   // Ensure the sub-vectors x_gf, v_gf, and e_gf know the location of the
   // data in S. This operation simply updates the Memory validity flags of
   // the sub-vectors to match those of S.
   x_gf.SyncAliasMemory(S);
   sv_gf.SyncAliasMemory(S);
   v_gf.SyncAliasMemory(S);
   ste_gf.SyncAliasMemory(S);

   /* Create Lagrangian Low Order Solver Object */
   mfem::hydroLO::LagrangianLOOperator hydro(S.Size(), H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, rho0_gf, m, problem_class, offset, use_viscosity, _mm, CFL);
   hydro.SetMVOption(mv_option);
   hydro.SetFVOption(fv_option);
   hydro.SetElasticity(true);
   hydro.SetShearModulus(problem_class->get_shear_modulus());

   /* If the mesh hasn't moved, F should be the identity matrix */
   DenseMatrix F(dim), I(3), res(3);
   I = 0.; for (int i = 0; i < 3; i++) { I(i,i) = 1.; }
   hydro.elastic.ComputeAvgF(0, F);
   if (F.NumCols() != 3 && F.NumRows() != 3)
   {
      cout << "Improper dimensions of F.\n";
      return 1;
   }
   res = F;
   res -= I;
   if (res.FNorm() > tol)
   {
      cout << "F is not the identity matrix.\n";
      return 1;
   }

   /* If F is the identity, then es = 0 */
   DenseMatrix sigma(3);
   double e_sheer = -1.;
   hydro.elastic.ComputeS(0, rho0_gf[0], sigma);
   e_sheer = hydro.elastic.e_sheer(0);
   if (e_sheer > tol)
   {
      cout << "Sheer energy is not zero.\n";
      return 1;
   }

   /* Move with constant mesh velocity 1, F should still be identity */
   ParGridFunction dx_gf(&H1FESpace);
   dx_gf = 1.;
   const double dt = .01;
   x_gf.Add(dt, dx_gf);
   x_gf.SyncAliasMemory(S);
   pmesh->NewNodes(x_gf, false);

   hydro.elastic.ComputeAvgF(0, F);
   res = F;
   res -= I;
   if (res.FNorm() > tol)
   {
      cout << "Constant Velocity Field: F is not the identity matrix.\n";
      return 1;
   }
   // Reset x_gf
   x_gf = x0_gf;
   x_gf.SyncAliasMemory(S);
   pmesh->NewNodes(x_gf, false);

   


   return 0;
}

int test_elasticity()
{
   return 0;
}

