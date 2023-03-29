#include "mfem.hpp"
#include "laglos_solver.hpp"
#include <cassert>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;
using namespace mfem;

/* ---------------- Parameters to be used for test ---------------- */
// Linear velocity field
const double a = 1.,
             b = 1.,
             c = 1.,
             d = 1.,
             e = 1.,
             f = 1.;

// Problem
const int dim = 2;
const int problem = 8;     // This value doesn't matter (1-8)
int order_mv = 2;          // Order of mesh movement approximation space
int order_u = 0;
const string flag = "testing";
double tol = 1e-12;
const char *mesh_file = "../data/ref-square.mesh";
int mesh_refinements = 6;
bool use_viscosity = true; // Doesn't matter
bool _mm = true;           // Doesn't matter
double CFL = 0.5;          // Doesn't matter
static int num_procs = 0;
static int myid = -1;
bool suppress_test_output = false;
bool create_table = false;
int lower_refinement = 2;
int upper_refinement = 6;
/* ---------------- End Parameters ---------------- */

// Default reference values
static double _error_def = 0.;
static int _num_cells_def = 0;
int test_mesh_movement(double & _error = _error_def, int & _num_cells = _num_cells_def);
void create_error_convergence_table();

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

   // Handle either CMake test or error table production
   if (create_table)
   {
      // We will be outputting a tabulate table to console,
      // which is why we want to suppress the test's output.
      suppress_test_output = true;
      create_error_convergence_table();
      return 1;
   }
   else
   {
      // This is the call for the CMake test
      return test_mesh_movement();
   }
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

int test_mesh_movement(double & _error, int & _num_cells)
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
   if (!suppress_test_output)
   {
      cout << "num cells: " << pmesh->GetNE() << endl;
      cout << "num interior faces: " << pmesh->GetNFbyType(FaceType::Interior) << endl;
      cout << "num boundary faces: " << pmesh->GetNFbyType(FaceType::Boundary) << endl;
      cout << "num total faces: " << pmesh->GetNumFaces() << endl;
   }
   

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-TestMM." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   // Construct blockvector
   const int Vsize_l2 = L2FESpace.GetVSize();
   const int Vsize_l2v = L2VFESpace.GetVSize();
   const int Vsize_h1 = H1FESpace.GetVSize();
   Array<int> offset(2);
   offset[0] = 0;
   offset[1] = offset[0] + Vsize_h1;
   BlockVector S(offset, Device::GetMemoryType());

   // Pair with corresponding gridfunctions
   // Only need position and velocity for testing
   ParGridFunction x_gf, mv_gf;
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   mv_gf.MakeRef(&H1FESpace, S, offset[1]);

   pmesh->SetNodalGridFunction(&x_gf);
   x_gf.SyncAliasMemory(S);
   mv_gf.SyncAliasMemory(S);

   // Just leave templated for hydro construction
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   
   mfem::hydrodynamics::LagrangianLOOperator<dim, problem> hydro(H1FESpace, L2FESpace, L2VFESpace, m, use_viscosity, _mm, CFL);

   const double t = 0;
   const double dt = 0.1;
   
   hydro.compute_node_velocity_LS(S, t, dt, flag, &velocity_exact);
   if (!suppress_test_output)
   {
      cout << "Done computing node velocity.\n";
   }
   

   // Compute error
   double num = 0, denom = 0;
   Vector vel_approx(dim), vel_exact(dim), vel_error(dim), vertex_x(dim);
   int num_int_vertices = 0;
   // Iterate over vertices to compute error
   for (int vertex = 0; vertex < pmesh->GetNV(); vertex++)
   {
      if (!hydro.IsBdrVertex(vertex))
      {
         // Interior vertex
         vel_approx = 0.;
         vel_exact = 0.;
         num_int_vertices += 1;

         // Get Values
         hydro.get_node_position(S, vertex, vertex_x);
         hydro.get_node_velocity(S, vertex, vel_approx);
         velocity_exact(vertex_x, t, vel_exact);
         subtract(vel_approx, vel_exact, vel_error);

         // Add contribution to error computation
         denom += vel_exact.Norml2();
         num += vel_error.Norml2();

         /* Optionally output information */
         
         if (vel_error.Norml2() > 1e-9 && !suppress_test_output)
         {
            cout << "------------------------------\n";
            cout << "Vertex: " << vertex << endl;
            cout << "Vertex position:\n";
            vertex_x.Print(cout);
            cout << "Exact velocity:\n";
            vel_exact.Print(cout);
            cout << "Approx velocity:\n";
            vel_approx.Print(cout);
            
            cout << "\tvel_Error: \n";
            vel_error.Print(cout);
            cout << "error: " << vel_error.Norml2() << endl;
         }
      } 
   } // End vertex iterator

   _error = num / denom;
   _num_cells = pmesh->GetNE();
   if (!suppress_test_output)
   {
      cout << "Error: " << _error << endl;
      // cout << "Denom: " << denom << ", Num interior vertices: " << num_int_vertices << endl;
      // cout << "Num: " << num << endl;
   }

   // Delete remaining pointers
   delete pmesh, m;

   // Point dangling pointers to NULL
   pmesh = NULL;
   m = NULL;

   if (_error < tol)
   {
      return 0; // Test passed
   }
   else 
   {
      return 1; // Test failed
   }
}

/* 
* Function: create_error_convergence_table
* 
* Purpose: The purpose of this function is to iterate over different refinement
* levels, run the test on the mesh movement, and retrieve the corresponding error
* to create a table.  We should expect to see the error reducing as we refine the mesh.
*/
void create_error_convergence_table()
{
   double error = 0;
   int _num_cells = 0;
   int ref_initial = 2, ref_final = 7; // Change these

   // Table information
   cout << "\\begin{center}\n"
        << "\\begin{tabular}{\n"
        << "l\n"
        << "S[round-mode=places, round-precision=3, table-format=2.3e-2]\n"
        << "}\n"
        << "\\hline\n"
        << "{\\# cells} & {Error} \\\\\n"
        << "\\hline\n";

   // Loop over range of
   for (int _ref = lower_refinement; _ref <= upper_refinement; _ref++)
   {
      // Reset error
      error = 0;
      _num_cells = 0;

      // Change global constant representing mesh refinements
      mesh_refinements = _ref;

      // Run test
      test_mesh_movement(error, _num_cells);
      // cout << "Testing mesh movement at " << _ref << " refinements.\n";
      // cout << "The total number of cells in our mesh is: " << _num_cells << endl;
      // cout << "The corresponding error is: " << error << endl;
      cout << "\t" << _num_cells << " & " << error << " \\\\\n";
   }
   cout << "\\hline\n"
        << "\\end{tabular}\n"
        << "\\end{center}\n";
   
}

// \begin{tabular}{
//   l
//   S[round-mode=places, round-precision=3, table-format=2.3]
//   S[round-mode=places, round-precision=2, table-format=2.2]
//   S[round-mode=places, round-precision=3, table-format=2.3]
//   S[round-mode=places, round-precision=2, table-format=2.2]}
// \hline
//    {\# dof} & {L1 Error} & {Rate}     &   {L2 Error} & {Rate} \\
// \hline
//       16 &   1.43065  & {---}      &   1.34817  & {---}  \\
//       64 &   1.1216   & 0.351112 &   0.968713 & 0.476861\\
//      256 &   0.826564 & 0.440358 &   0.732866 & 0.402521  \\
//     1024 &   0.491592 & 0.749666 &   0.495077 & 0.565896   \\
//     4096 &   0.318258 & 0.627263 &   0.367941 & 0.428179  \\
//    16384 &   0.192916 & 0.722223 &   0.274031 & 0.425134 \\
// \hline
// \end{tabular}
