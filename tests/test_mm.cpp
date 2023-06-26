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
const double a = 0.,
             b = -1.,
             c = 0.,
             d = 1.,
             e = 0.,
             f = 0.;

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
int test_mesh_movement(double & _error = _error_def, int & _num_cells = _num_cells_def);
int test_mass_conservation(double & _error = _error_def);
int test_area_conservation(double & _error = _error_def, int & _num_cells = _num_cells_def);
void create_error_convergence_table();
void move_corner_nodes_with_corrected_velocity(Vector &S, 
                                               const double & t,
                                               const double &dt,
                                               const string flag="NA", 
                                               void (*test_vel)(const Vector&, const double&, Vector&)=NULL); 

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
      // return test_mesh_movement();
      return test_area_conservation();
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
   const double dt = 1.;
   
   hydro.compute_node_velocities(S, t, dt, flag, &velocity_exact);
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
   delete pmesh;
   delete m;

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
      // test_mesh_movement(error, _num_cells);
      test_area_conservation(error, _num_cells);
      // cout << "Testing mesh movement at " << _ref << " refinements.\n";
      // cout << "The total number of cells in our mesh is: " << _num_cells << endl;
      // cout << "The corresponding error is: " << error << endl;
      cout << "\t" << _num_cells << " & " << error << " \\\\\n";
   }
   cout << "\\hline\n"
        << "\\end{tabular}\n"
        << "\\end{center}\n";
   
}


/* Just moving nodes, averaging faces */
int test_area_conservation(double & _error, int & _num_cells)
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
   _num_cells = pmesh->GetNE();
   if (!suppress_test_output)
   {
      cout << "num cells: " << _num_cells << endl;
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
   mv_gf.SyncAliasMemory(S);   

   // zero coefficient for sv and ste
   ConstantCoefficient zero(0.0);
   sv_gf.ProjectCoefficient(zero);
   sv_gf.SyncAliasMemory(S);
   ste_gf.ProjectCoefficient(zero);
   ste_gf.SyncAliasMemory(S);

   // Prescribe initial velocity profile on cells
   VectorFunctionCoefficient v_coeff(dim, velocity_exact);
   v_gf.ProjectCoefficient(v_coeff);
   v_gf.SyncAliasMemory(S);

   // Just leave templated for hydro construction
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   
   mfem::hydrodynamics::LagrangianLOOperator<dim, problem> hydro(H1FESpace, L2FESpace, L2VFESpace, m, use_viscosity, _mm, CFL);

   double t = 0;
   const double dt = .9;
   Array<double> cell_areas(L2FESpace.GetNE());
   Array<double> bmn_sum_cells(L2FESpace.GetNE());

   // Store initial mass of cells
   for (int ci = 0; ci < L2FESpace.GetNE(); ci++)
   {
      double k = pmesh->GetElementVolume(ci);

      cell_areas[ci] = k;
      assert(k >= 0.);
   }

   /* ************************
   Compute nodal velocities
   *************************** */ 
   Vector vertex_x(dim), vertex_v(dim), v_pre_update(dim), v_post_update(dim), res(dim);
   DenseMatrix C(dim);
   Vector D(dim);

   hydro.compute_node_velocities(S, t, dt, flag, &velocity_exact);
   /* Fill faces with average velocities */
   // hydro.fill_face_velocities_with_average(S, flag, &velocity_exact);
   hydro.compute_corrective_face_velocities(S, t, dt, flag, &velocity_exact);
   hydro.fill_center_velocities_with_average(S, flag, &velocity_exact);

   /* ************************
   Compute cell bmn sum (RHS)
   *************************** */ 
   double num = 0.;
   _error = 0.;
   double error_denom = 0.;
   double bmn_sum = 0.;

   Array<int> fids, oris, row;
   mfem::Mesh::FaceInformation FI;
   Vector secant(dim), n_vec(dim), n_int(dim), vdof1_x(dim), vdof2_x(dim), face_x(dim), face_v(dim);
   int face_vdof1, face_vdof2, face_dof;
   int c, cp;
   Vector Uc(dim+2), Ucp(dim+2), c_vec(dim), Vf(dim);

   // vector for plotting the error of each cell
   std::vector<double> cells_py(L2FESpace.GetNE());
   std::vector<double> cell_errors_py(L2FESpace.GetNE());

   for (int ci = 0; ci < L2FESpace.GetNE(); ci++)
   {
      // cout << "===== Computing cell bmn sum for cell: " << ci << " =====\n";
      // Compute sum of bmn on cell
      pmesh->GetElementEdges(ci, fids, oris);

      // iterate over faces
      bmn_sum = 0.;
      for (int j = 0; j < fids.Size(); j++)
      {
         Vf = 0.;
         int face = fids[j];
         FI = pmesh->GetFaceInformation(face);
         // cout << "face: " << face << endl;

         // Compute intermediate face velocity on fly
         c = FI.element[0].index;
         cp = FI.element[1].index;
         // cout << "cell c: " << c << ", cell cp: " << cp << endl;
         hydro.GetCellStateVector(S, c, Uc);

         H1FESpace.GetFaceDofs(fids[j], row);
         int face_vdof1 = row[1], face_vdof2 = row[0], face_dof = row[2];
         hydro.get_node_position(S, face_vdof1, vdof1_x);
         hydro.get_node_position(S, face_vdof2, vdof2_x);
         hydro.get_node_position(S, face_dof, face_x);

         // if (FI.IsInterior())
         // {
         //    // this vector is only needed for interior faces
         //    hydro.GetCellStateVector(S, cp, Ucp);
         //    // Get normal, d, and |F|
         //    hydro.CalcOutwardNormalInt(S, c, face, c_vec);
         //    double c_norm = c_vec.Norml2();
         //    c_vec *= 2.; // Needed to accomodate for the 1/2 in 4.7 that isn't present in 5.7a
         //    double F = c_vec.Norml2();
         //    n_vec = c_vec;
         //    n_vec /= F;

         //    double _d = hydro.compute_lambda_max(Uc, Ucp, n_vec) * c_norm;

         //    Vf = hydro.velocity(Uc);
         //    Vf += hydro.velocity(Ucp);
         //    Vf *= 0.5;
         //    double coeff = _d * (Ucp[0] - Uc[0]) / F;
         //    Vf.Add(coeff, n_vec);
         // }
         // else
         // {
         //    assert(FI.IsBoundary());
         //    Vf = hydro.velocity(Uc);
         // }
         
         velocity_exact(face_x, t, Vf);

         // if (!suppress_test_output)
         // {
         //    cout << "face: " << fids[j] << ", face x: " << endl;
         //    face_x.Print(cout);
         // }

         // subtract(vdof2_x, vdof1_x, secant);
         // double face_length = secant.Norml2();
         // secant /= face_length;
         // Vector sec_norm(dim);
         // sec_norm = secant;
         // hydro.Orthogonal(sec_norm);
         // n_vec = secant;
         // hydro.Orthogonal(n_vec);
         hydro.CalcOutwardNormalInt(S, ci, face, n_int);
         n_vec = n_int;
         double F = n_vec.Norml2();
         n_vec /= F;

         // hydro.get_node_velocity(S, face_dof, face_v);
         // hydro.get_intermediate_face_velocity(fids[j], Vf);
         // cout << "face_x:\n";
         // face_x.Print(cout);
         // cout << "Vf:\n";
         // Vf.Print(cout);
         // cout << "outward normal vector:\n";
         // n_vec.Print(cout);
         // cout << "outward normal vector secant:\n";
         // sec_norm.Print(cout);

         // double bmn = face_v * n_vec;
         double bmn = Vf * n_vec;
         bmn *= F;
         // bmn *= n_vec.Norml2();
         // cout << "bmn: " << bmn << endl << endl;
         bmn_sum += bmn;

         // if (!suppress_test_output)
         // {
         //    cout << "Face velocity:\n";
         //    face_v.Print(cout);
         //    cout << "n_vec:\n";
         //    n_vec.Print(cout);
         //    cout << "face length: " << face_length << endl;
         //    cout << "bmn: " << bmn << endl;
         // }
         
      } // End face iterator
      bmn_sum_cells[ci] = bmn_sum;
   }

   /* ************************
   Move the mesh according to 
   previously computed nodal 
   velocities
   *************************** */ 
   add(x_gf, dt, mv_gf, x_gf);
   pmesh->NewNodes(x_gf, false);
   t += dt;
   if (!suppress_test_output)
   {
      cout << "Done moving mesh.\n";
   }

   cout << "S:\n";
   S.Print(cout);
   /* ************************
   Compute cell area error
   *************************** */ 
   for (int ci = 0; ci < L2FESpace.GetNE(); ci++)
   {
      if (!suppress_test_output)
      {
         cout << "Computing cell area error for cell: " << ci << endl;
      }
      // Compute cell area side
      double k = pmesh->GetElementVolume(ci);
      double _bmn_approx = (k - cell_areas[ci]) / dt;
      cout << "k: " << k << endl;
      cout << "cell area: " << cell_areas[ci] << endl;
      cout << "dt: " << dt << endl;

      cout << "Cell wise comparison:\n";
      cout << "cell area difference: " << _bmn_approx << ", bmn_sum: " << bmn_sum_cells[ci] << endl;

      double val = abs(_bmn_approx - bmn_sum_cells[ci]);
      num += val;
      error_denom += abs(bmn_sum_cells[ci]);
      cells_py[ci] = ci;
      cell_errors_py[ci] = val;
   }

   if (error_denom < tol)
   {
      _error = 0.;
   }
   else
   {
      _error = abs(num / error_denom);
   }

   if (!suppress_test_output)
   {
      cout << "Total area error is: " << _error << endl;
   }

   /* ************************
   Visualization
   *************************** */ 
   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-TestMM-post-move." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   // Plots cell errors using Python's matplotlib library
   if (!suppress_test_output)
   {
      plt::figure_size(1200, 780);
      plt::title("Error by cell");
      plt::plot(cells_py, cell_errors_py, "og");
      plt::xlabel("cell");
      plt::ylabel("error");
      plt::show();
   }

   if (_error < tol)
   {
      return 0; // Test passed
   }
   else 
   {
      return 1; // Test failed
   }
}
