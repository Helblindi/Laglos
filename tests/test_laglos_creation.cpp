#include "mfem.hpp"
#include "laglos_solver.hpp"
#include <cassert>
#include <fstream>
#include <sstream>
#include <cmath>
#include "matplotlibcpp.h"

using namespace std;
using namespace mfem;
namespace plt = matplotlibcpp;

/* ---------------- Parameters to be used for tests ---------------- */
// Linear velocity field
const double a = 1.,
             b = 1.,
             c = 1.,
             d = 1.,
             e = -1.,
             f = 1.;

// Problem
const int dim = 2;
const int problem = 8;     // This value doesn't matter (1-8)
int order_mv = 2;          // Order of mesh movement approximation space
int order_u = 0;
const string flag = "testing";
double tol = 1e-12;
const char *mesh_file = "../data/ref-square.mesh";
int mesh_refinements = 1;
bool use_viscosity = true; // Doesn't matter
bool _mm = false;          // Doesn't matter
double CFL = 0.5;          // Doesn't matter
static int num_procs = 0;
static int myid = -1;
bool suppress_test_output = false;
bool create_table = false;
int lower_refinement = 2;
int upper_refinement = 7;
/* ---------------- End Parameters ---------------- */

// Default reference values
static double _error_def = 0.;
static int _num_cells_def = 0;
void velocity_exact(const Vector &x, const double & t, Vector &v);
int test_mesh_initiation();

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
   args.AddOption(&create_table, "-ct", "--create-table", 
                  "-no-ct", "--no-create-table", 
                  "Option to use to create a table.");
   args.AddOption(&lower_refinement, "-lr", "--lower-refinement", 
                  "Lower bound for refinements for table creation.");
   args.AddOption(&upper_refinement, "-ur", "--upper-refinement", 
                  "Upper bound for refinements for table creation.");
   args.Parse();

   test_mesh_initiation();
}

void velocity_exact(const Vector &x, const double & t, Vector &v)
{
   if (dim != 2)
   {
      MFEM_ABORT("Dimension must be equal to 2 for mesh motion test.\n");
   }
   v[0] = a * x[0] + b * x[1] + c;
   v[1] = d * x[0] + e * x[1] + f;
}

int test_mesh_initiation()
{
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
   L2_FECollection L2FEC(order_u, dim, BasisType::Positive);

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);

   // Output information
   pmesh->ExchangeFaceNbrData();
   // if (!suppress_test_output)
   // {
   //    cout << "num cells: " << pmesh->GetNE() << endl;
   //    cout << "num interior faces: " << pmesh->GetNFbyType(FaceType::Interior) << endl;
   //    cout << "num boundary faces: " << pmesh->GetNFbyType(FaceType::Boundary) << endl;
   //    cout << "num total faces: " << pmesh->GetNumFaces() << endl;
   // }
   
   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-TestInitiation." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

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
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   mv_gf.MakeRef(&H1FESpace, S, offset[1]);
   sv_gf.MakeRef(&L2FESpace, S, offset[2]);
   v_gf.MakeRef(&L2VFESpace, S, offset[3]);
   ste_gf.MakeRef(&L2FESpace, S, offset[4]);

   pmesh->SetNodalGridFunction(&x_gf);
   x_gf.SyncAliasMemory(S);

   Vector mv0(dim);
   mv0[0] = 1.1;
   mv0[1] = 1.2;
   VectorConstantCoefficient mv_coeff(mv0);
   mv_gf.ProjectCoefficient(mv_coeff);
   mv_gf.SyncAliasMemory(S);

   // Initialize specific volume, velocity, and specific total energy
   ConstantCoefficient n_one(-1.0);
   sv_gf.ProjectCoefficient(n_one);
   sv_gf.SyncAliasMemory(S);

   Vector v0(dim);
   v0[0] = -1.1;
   v0[1] = -1.2;
   VectorFunctionCoefficient v_coeff(dim, velocity_exact);
   v_gf.ProjectCoefficient(v_coeff);
   v_gf.SyncAliasMemory(S);

   ConstantCoefficient n_two(-2.0);
   ste_gf.ProjectCoefficient(n_two);
   ste_gf.SyncAliasMemory(S);

   // BlockVector instantiated
   cout << "S: \n";
   S.Print(cout);

   // Just leave templated for hydro construction
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   
   mfem::hydrodynamics::LagrangianLOOperator<dim, problem> hydro(H1FESpace, L2FESpace, L2VFESpace, m, use_viscosity, _mm, CFL);

   const double t = 0;
   const double dt = 1.;

   // TODO: Test Getter and Setter of state variables
   Vector U(dim+2), new_vals(dim+2);
   int index = 1;

   cout << "\n-------------------------\n";
   cout << "U before changing values:\n";
   hydro.GetCellStateVector(S, index, U);
   U.Print(cout);

   new_vals = 10.;
   hydro.SetCellStateVector(S, index, new_vals);
   
   cout << "U after changing values:\n";
   hydro.GetCellStateVector(S, index, U);
   U.Print(cout);
   cout << "-------------------------\n\n";

   cout << "S: \n";
   S.Print(cout);

   cout << "\n-------------------------\n";
   cout << "Testing retrieval of intermediate face velocities.\n";
   hydro.compute_intermediate_face_velocities(S, t, flag, &velocity_exact);
   int face = 1;
   Vector face_vel(dim);
   hydro.get_intermediate_face_velocity(face, face_vel);
   cout << "face vel for face: " << face << endl;
   face_vel.Print(cout);
   cout << "-------------------------\n\n";

   // Delete remaining pointers
   delete pmesh;
   delete m;

   // Point dangling pointers to NULL
   pmesh = NULL;
   m = NULL;

   return 0;
}
