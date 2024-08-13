#include "mfem.hpp"
#include "laglos_solver.hpp"
#include "test_problems_include.h"
#include "var-config.h"
#include <cassert>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;
using namespace mfem;
using namespace hydrodynamics;

/* ------------- Problem Parameters ------------- */
const int dim = 2;
const int problem = 1;     // This value doesn't matter (1-8)
int order_mv = 2;          // Order of mesh movement approximation space
int order_u = 0;
const string flag = "testing";
double tol = 1e-12;
const char *mesh_file_location = "/tests/test-square.mesh";
std::string result = std::string(LAGLOS_DIR) + std::string(mesh_file_location);
const char* mesh_file = result.c_str();
int mesh_refinements = 1;
bool use_viscosity = true; // Doesn't matter
bool _mm = true;          // Doesn't matter
int mv_option = 3;
int fv_option = 2;
double CFL = 0.5;          // Doesn't matter
static int num_procs = 0;
static int myid = -1;
bool suppress_test_output = false;
bool create_table = false;
int lower_refinement = 2;
int upper_refinement = 7;
/* ---------------- End Parameters ---------------- */

void velocity_exact(const Vector &x, const double & t, Vector &v);
/* Define specific test functions here */
int TestWeightedCellAverageVelocityAtNode();

int main(int argc, char *argv[])
{
   // Initialize MPI.
   Mpi::Init();
   num_procs = Mpi::WorldSize();
   myid = Mpi::WorldRank();

   // Set precision for all output values
   cout.precision(12);

   int d = 0;

   d += TestWeightedCellAverageVelocityAtNode();

   // Run individual tests
   // Should any test fail, a 1 is returned.
   // If d is nonzero by the end of main, the whole test fails

   return d;
}

void velocity_exact(const Vector &x, const double & t, Vector &v)
{
   v = x;
}

int TestWeightedCellAverageVelocityAtNode()
{
   cout << "Test Weighted Cell Average Velocity At Node\n";

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
   Array<int> offset(6);
   offset[0] = 0;
   offset[1] = offset[0] + Vsize_h1;
   offset[2] = offset[1] + Vsize_h1;
   offset[3] = offset[2] + Vsize_l2;
   offset[4] = offset[3] + Vsize_l2v;
   offset[5] = offset[4] + Vsize_l2;
   BlockVector S(offset, Device::GetMemoryType());

   // Pair with corresponding gridfunctions
   // Only need position and velocity for testing
   ParGridFunction x_gf, mv_gf, sv_gf, v_gf, ste_gf;
   ParGridFunction v_cr_gf(&CRFESpace);
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   mv_gf.MakeRef(&H1FESpace, S, offset[1]);
   sv_gf.MakeRef(&L2FESpace, S, offset[2]);
   v_gf.MakeRef(&L2VFESpace, S, offset[3]);
   ste_gf.MakeRef(&L2FESpace, S, offset[4]);

   pmesh->SetNodalGridFunction(&x_gf);
   x_gf.SyncAliasMemory(S);

   Vector one(dim);
   one = 1.;
   Vector zero(dim);
   zero = 0.;

   VectorConstantCoefficient zero_coeff(zero);
   mv_gf.ProjectCoefficient(zero_coeff);
   mv_gf.SyncAliasMemory(S);

   // Initialize specific volume, velocity, and specific total energy
   ConstantCoefficient one_const_coeff(1.0);
   sv_gf.ProjectCoefficient(one_const_coeff);
   sv_gf.SyncAliasMemory(S);

   VectorFunctionCoefficient v_exact_coeff(dim, &velocity_exact);
   v_gf.ProjectCoefficient(v_exact_coeff);
   v_gf.SyncAliasMemory(S);

   ste_gf.ProjectCoefficient(one_const_coeff);
   ste_gf.SyncAliasMemory(S);

   // Just leave templated for hydro construction
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(one_const_coeff));
   m->Assemble();

   ProblemBase<dim> * problem_class = new SodProblem<dim>();

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);
   hydro.SetMVOption(mv_option);
   hydro.SetFVOption(fv_option);
   Vector vbar(dim), pos(dim);

   for (int node = 0; node < H1FESpace.GetNVDofs(); node++)
   {
      cout << "node: " << node << ", coord: ";
      hydro.GetNodePositionFromBV(S, node, pos);
      pos.Print(cout);
      cout << "vel: ";
      hydro.ComputeWeightedCellAverageVelocityAtNode(S, node, vbar);
      vbar.Print(cout);
   }

   Vector U(dim+2), cell_v(dim);
   for (int cell = 0; cell < L2FESpace.GetNDofs(); cell++)
   {
      hydro.GetCellStateVector(S, cell, U);
      problem_class->velocity(U, cell_v);
      cout << "cell v for cell " << cell << ": ";
      cell_v.Print(cout);
   }

   cout << "S: ";
   S.Print(cout);
}