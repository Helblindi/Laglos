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
const char *mesh_file_location = "/data/ref-square.mesh";
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

   // Parse command line arguments
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&mesh_refinements, "-r", "--refinement", "Number of mesh refinements.");

   args.Parse();

   int d = 0;

   // Run individual tests
   // Should any test fail, a 1 is returned.
   // If d is nonzero by the end of main, the whole test fails

   return d;
}