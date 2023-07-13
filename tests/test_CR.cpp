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
void RT_vel(mfem::hydrodynamics::LagrangianLOOperator<dim, problem> & hydro,
            const ParFiniteElementSpace & CRFESpace,
            const ParMesh * pmesh,
            const Table & element_face, 
            const int & cell,
            const int & node, 
            Vector &vel);

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

   int d = 0;
   d += test_mesh_initiation();
   return d;
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

/******************************************************************************************
* This is a test case to validate the Rannacher-Turek velocity reconstruction function.
* It assumes that the given test_velocity function is defined as:
*        C      x  +   D
*     | 1 1  | |x|   | 1 |
*     |      | | | + |   |
*     | 1 -1 | |y|   | 1 |
* Using this velocity, a LagrangianLOOperator class is instantiated, face velocities are
* computed exactly using this velocity, and then the velocity is computed at a given node.
******************************************************************************************/
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
   CrouzeixRaviartFECollection CRFEC;

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);
   ParFiniteElementSpace CRFESpace(pmesh, &CRFEC, dim);

   // Output information
   pmesh->ExchangeFaceNbrData();
   
   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-TestCR." << setfill('0') << setw(6) << myid;
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
   VectorConstantCoefficient v0_coeff(v0);
   v_gf.ProjectCoefficient(v0_coeff);
   v_gf.SyncAliasMemory(S);

   ConstantCoefficient n_two(-2.0);
   ste_gf.ProjectCoefficient(n_two);
   ste_gf.SyncAliasMemory(S);

   // Just leave templated for hydro construction
   ConstantCoefficient one(1.0);
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(one));
   m->Assemble();
   mfem::hydrodynamics::LagrangianLOOperator<dim, problem> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, use_viscosity, _mm, CFL);

   Vector _vel(dim), true_vel(dim);
   double t = 0., dt = 0.5;

   hydro.compute_intermediate_face_velocities(S, t, "testing", &velocity_exact);

   hydro.compute_node_velocities(S,t,dt);
   hydro.fill_face_velocities_with_average(S);
   hydro.fill_center_velocities_with_average(S);

   cout << "Printing S:\n";
   S.Print(cout);

   add(x_gf, dt, mv_gf, x_gf);
   pmesh->NewNodes(x_gf, false);

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-TestCR-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   // Delete remaining pointers
   delete pmesh;
   delete m;

   // Point dangling pointers to NULL
   pmesh = NULL;
   m = NULL;

   // if (true_vel[0] == _vel[0] && true_vel[1] == _vel[1]) { 
   //    return 0; 
   // }
   // else { 
   //    return 1; 
   // }
   // TODO: Properly implemment test
   return 0;
}