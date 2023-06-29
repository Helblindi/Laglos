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
void velocity_exact(const Vector &x, const double & t, Vector &v);
int test_mesh_initiation();
void RT_vel(mfem::hydrodynamics::LagrangianLOOperator<dim, problem> & hydro,
            const ParFiniteElementSpace & CRFESpace,
            const ParMesh * pmesh,
            const Table & cell_face, 
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
   CrouzeixRaviartFECollection CRFEC;

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);
   ParFiniteElementSpace CRFESpace(pmesh, &CRFEC, dim);

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
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   mfem::hydrodynamics::LagrangianLOOperator<dim, problem> hydro(H1FESpace, L2FESpace, L2VFESpace, m, use_viscosity, _mm, CFL);

   // Output information on CRFESpace
   cout << "Crouzeix-Raviart info:\n";
   cout << "CRFESpace.GetVSize(): " << CRFESpace.GetVSize() << endl;
   cout << "CRFESpace.TrueVSize(): " << CRFESpace.TrueVSize() << endl;
   cout << "CRFESpace.GlobalTrueVSize(): " << CRFESpace.GlobalTrueVSize() << endl;
   cout << "CRFESpace.GetNDofs(): " << CRFESpace.GetNDofs() << endl;
   cout << "CRFESpace.GetNVDofs(): " << CRFESpace.GetNVDofs() << endl;
   cout << "CFFESpace.GetNE(): " << CRFESpace.GetNE() << endl;

   IntegrationPoint test_ip;
   test_ip.Set2(1.,0.);
   Vector shape(4);
   DenseMatrix dshape(4,2);
   dshape = 0.;
   shape = 0.;

   const FiniteElement * fe = CRFESpace.GetFE(2);
   fe->CalcShape(test_ip, shape);
   fe->CalcDShape(test_ip, dshape);
   IntegrationRule ir = fe->GetNodes();
   cout << "====== ir: " << ir.GetNPoints() << endl;
   // ir.Print(cout);

   cout << "shape:\n";
   shape.Print(cout);

   cout << "dshape:\n";
   dshape.Print(cout);

   // element to node table
   Table cell_face = CRFESpace.GetElementToDofTable();
   Array<int> cell_face_row;
   int row_length = 0;
   for (int cell = 0; cell < CRFESpace.GetNE(); cell++)
   {  
      cell_face.GetRow(cell, cell_face_row);
      row_length = cell_face_row.Size();
      cout << "cell: " << cell << ", row:\n";
      cell_face_row.Print(cout);
   }

   int node = 1;
   Vector _vel(dim);
   RT_vel(hydro, CRFESpace, pmesh, cell_face, 1, node, _vel);

   const double t = 0;
   const double dt = 1.;

   // Delete remaining pointers
   delete pmesh;
   delete m;

   // Point dangling pointers to NULL
   pmesh = NULL;
   m = NULL;

   return 0;
}

void RT_vel(mfem::hydrodynamics::LagrangianLOOperator<dim, problem> & hydro,
            const ParFiniteElementSpace & CRFESpace,
            const ParMesh * pmesh,
            const Table & cell_face, 
            const int & cell,
            const int & node, 
            Vector &vel)
{
   // Get DoFs corresponding to cell
   Array<int> cell_face_row;
   cell_face.GetRow(cell, cell_face_row);
   int row_length = cell_face_row.Size();
   cout << "faces on cell:\n";
   for (int i = 0; i < row_length; i++)
   {
      cout << "face: " << cell_face_row[i] << endl;
   }

   // TODO: Get integration point corresponding to node location
   Array<int> verts;
   pmesh->GetElementVertices(cell, verts);
   // cout << "verts for element 0: \n";
   // for (int i = 0; i < verts.Size(); i++)
   // {
   //    cout << "vert: " << verts[i] << endl;
   // }

   IntegrationPoint test_ip;
   if (node == verts[0]) {
      test_ip.Set2(0., 0.);
      // cout << "ip is (0,0)\n";
   } else if (node == verts[1]) {
      test_ip.Set2(1., 0.);
      // cout << "ip is (1,0)\n";
   } else if (node == verts[2]) {
      test_ip.Set2(1., 1.);
      // cout << "ip is (1,1)\n";
   } else if (node == verts[3]) {
      test_ip.Set2(0., 1.);
      // cout << "ip is (0,1)\n";
   } else { 
      // Invalid if node is not contained in cell
      MFEM_ABORT("Incorrect node provided.\n"); 
   }

   // Evaluate reference shape functions at integration point
   Vector shapes(4);
   const FiniteElement * fe = CRFESpace.GetFE(cell);
   fe->CalcShape(test_ip, shapes);
   // cout << "shapes: \n";
   // shapes.Print(cout);

   // Sum over faces evaluating local velocity at vertex
   Vector Ci(dim);
   Ci = 0.;
   Vector face_vel(dim);
   Vector v_face_intermediate(CRFESpace.GetVSize());
   v_face_intermediate[0] = 1.;
   v_face_intermediate[12] = 1.;
   v_face_intermediate[1] = 2.;
   v_face_intermediate[13] = 2.;
   v_face_intermediate[2] = 3.;
   v_face_intermediate[14] = 3.;
   v_face_intermediate[3] = 4.;
   v_face_intermediate[15] = 4.;

   for (int face_it = 0; face_it < row_length; face_it++)
   {
      int face = cell_face_row[face_it];

      // Get intermediate face velocities corresponding to faces
      // hydro.get_intermediate_face_velocity(face, face_vel);
      face_vel[0] = v_face_intermediate[face];
      face_vel[1] = v_face_intermediate[face + CRFESpace.GetNDofs()];

      // cout << "face it: " << face_it << ", face: " << face << ", shape: " << shapes[face_it] << endl;
      // cout << "face_vel:\n";
      // face_vel.Print(cout);
      Ci.Add(shapes[face_it], face_vel);
   }

   Ci.Print(cout);
}