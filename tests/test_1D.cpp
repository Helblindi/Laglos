#include "mfem.hpp"
#include "var-config.h"
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
double a = 1.,
       b = 1.,
       c = 1.,
       d = 1.,
       e = 1.,
       f = 1.;

// Problem
const int dim = 1;
const int problem = 8;     // This value doesn't matter (1-8)
int order_mv = 2;          // Order of mesh movement approximation space
int order_u = 0;
const string flag = "testing";
double tol = 1e-12;
const char *mesh_file_location = "/data/ref-segment.mesh";
std::string result = std::string(LAGLOS_DIR) + std::string(mesh_file_location);
const char* mesh_file = result.c_str();
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
int test_1d_mesh();
int test_vel_field_1();

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


   int d = test_1d_mesh();

   test_vel_field_1();

   cout << "Must return 0 for test to pass.  d = " << d << endl;

   // Must return 0 for test to pass
   return d;

}

void velocity_exact(const Vector &x, const double & t, Vector &v)
{
   switch(dim)
   {
      case 1:
      {
         v[0] = a * x[0] + c;
         break;
      }
      case 2:
      {
         v[0] = a * x[0] + b * x[1] + c;
         v[1] = d * x[0] + e * x[1] + f;
         break;
      }
      case 3:
      {
         MFEM_ABORT("3D not implemented\n");
         break;
      }
      default:
      {
         MFEM_ABORT("Invalid dimension supplied\n");
      }
   }

   
}

int test_1d_mesh()
{  
   double _error = 0;
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

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-Test1D-basic." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

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
   mv0 = 1.;

   VectorFunctionCoefficient v_exact_coeff(dim, &velocity_exact);
   mv_gf.ProjectCoefficient(v_exact_coeff);
   mv_gf.SyncAliasMemory(S);

   // Initialize specific volume, velocity, and specific total energy
   ConstantCoefficient one(1.0);
   sv_gf.ProjectCoefficient(one);
   sv_gf.SyncAliasMemory(S);

   v_gf.ProjectCoefficient(v_exact_coeff);
   v_gf.SyncAliasMemory(S);

   ste_gf.ProjectCoefficient(one);
   ste_gf.SyncAliasMemory(S);

   // Just leave templated for hydro construction
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(one));
   m->Assemble();

   mfem::hydrodynamics::LagrangianLOOperator<dim, problem> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, use_viscosity, _mm, CFL);

   cout << "x_gf:\n";
   x_gf.Print(cout);

   const double t = 0;
   const double dt = 0.;
   Array<double> cell_areas(L2FESpace.GetNE());
   Array<double> bmn_sum_cells(L2FESpace.GetNE());
   // hydro.compute_intermediate_face_velocities(S, t,"NA", &velocity_exact);

   // Store initial mass of cells
   for (int ci = 0; ci < L2FESpace.GetNE(); ci++)
   {
      double k = pmesh->GetElementVolume(ci);

      cell_areas[ci] = k;
      assert(k >= 0.);
   }

   hydrodynamics::LagrangianLOOperator<dim, problem>::DofEntity _dof_entity;
   int edof;
   for (int dof = 0; dof < H1FESpace.GetNDofs(); dof++)
   {      
      hydro.GetEntityDof(dof, _dof_entity, edof);
      cout << "dof: " << dof << ", entity: " << _dof_entity << ", edof: " << edof << endl;
   }

   return 0;
}

/*
Purpose:
   To test Remark 5.3.
   This function tests a the RT computation of a continuous velocity field from a prescribed 
   velocity field that is linear in the first component and constant in the second.

   v_exact = | 1 0 | |x| + |1|
             | 0 0 | |y|   |1|
*/
int test_vel_field_1()
{
   a = 1., b = 0., c = 1., d = 0., e = 0., f = 1.;

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

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-Test1D-linear." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

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

   VectorFunctionCoefficient v_exact_coeff(dim, &velocity_exact);
   ParGridFunction v_exact_gf(&H1FESpace);
   v_exact_gf.ProjectCoefficient(v_exact_coeff);

   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   mv_gf.MakeRef(&H1FESpace, S, offset[1]);
   sv_gf.MakeRef(&L2FESpace, S, offset[2]);
   v_gf.MakeRef(&L2VFESpace, S, offset[3]);
   ste_gf.MakeRef(&L2FESpace, S, offset[4]);

   pmesh->SetNodalGridFunction(&x_gf);
   x_gf.SyncAliasMemory(S);

   Vector one(dim);
   one = 1.;
   VectorConstantCoefficient one_coeff(one);
   mv_gf.ProjectCoefficient(one_coeff);
   mv_gf.SyncAliasMemory(S);

   // Initialize specific volume, velocity, and specific total energy
   ConstantCoefficient one_const_coeff(1.0);
   sv_gf.ProjectCoefficient(one_const_coeff);
   sv_gf.SyncAliasMemory(S);

   v_gf.ProjectCoefficient(one_coeff);
   v_gf.SyncAliasMemory(S);

   ste_gf.ProjectCoefficient(one_const_coeff);
   ste_gf.SyncAliasMemory(S);

   // Just leave templated for hydro construction
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(one_const_coeff));
   m->Assemble();

   mfem::hydrodynamics::LagrangianLOOperator<dim, problem> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, use_viscosity, _mm, CFL);

   Vector _vel(dim), vec_res(dim);
   double t = 0., dt = 1.;
   DenseMatrix dm(dim);

   hydro.compute_intermediate_face_velocities(S, t, "testing", &velocity_exact);

   for (int node_it = 0; node_it < H1FESpace.GetNDofs() - L2FESpace.GetNDofs(); node_it++)
   {
      hydro.compute_node_velocity_RT(node_it, dt, vec_res);
      hydro.update_node_velocity(S, node_it, vec_res);
   }

   // hydro.compute_corrective_face_velocities(S, t, dt);
   hydro.fill_center_velocities_with_average(S);

   /* ************************
   Displace Velocities
   *************************** */ 
   socketstream vis_vh, vis_vexact;
   char vishost[] = "localhost";
   int visport = 19916;

   MPI_Barrier(pmesh->GetComm());

   int Wx = 0, Wy = 0;
   const int Ww = 350, Wh = 350;
   int offx = Ww+10, offy = Wh+45;


   hydrodynamics::VisualizeField(vis_vh, vishost, visport, mv_gf, "Mesh Velocity", Wx, Wy, Ww, Wh);

   Wx += offx;

   hydrodynamics::VisualizeField(vis_vexact, vishost, visport, v_exact_gf, "Exact Velocity", Wx, Wy, Ww, Wh);


   /* ************************
   Move the mesh according to previously computed nodal velocities
   *************************** */ 
   add(x_gf, dt, mv_gf, x_gf);
   pmesh->NewNodes(x_gf, false);
   t += dt;
   if (!suppress_test_output)
   {
      cout << "Done moving mesh.\n";
   }


   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-Test1D-linear-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   // Since we prescribe the exact velocity on the faces, mv_gf should be the same as v_exact
   Vector vel_err(mv_gf.Size());
   subtract(mv_gf, v_exact_gf, vel_err);
   double _error = vel_err.Norml2();

   if (abs(_error) < tol)
   {
      return 0; // Test passed
   }
   else 
   {
      cout << "error: " << _error << endl;
      return 1; // Test failed
   }
}

