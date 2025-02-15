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
using namespace hydroLO;

/* ---------------- Parameters to be used for tests ---------------- */
// Linear velocity field
double a = 1.,
       b = 1.,
       c = 1.,
       d = 1.,
       e = -0.5,
       f = 1.,
       aq = 0.,
       bq = 0.,
       dq = 0.,
       eq = 0.;

// Problem
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

void velocity_exact(const Vector &x, const double & t, Vector &v);
int tester();
int test_flux();
int test_vel_field_1();
int test_CSV_getter_setter(); // LagrangianLOOperator::GetCellStateVector tester
int test_sod_hydro();

int test_smooth_hydro();
void plot_mv_smooth();

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
   // d += test_flux();
   // d += test_vel_field_1();
   // d += test_CSV_getter_setter();
   // d += test_sod_hydro();
   // d += test_smooth_hydro();

   // plot_mv_smooth();

   return d;
}


void velocity_exact(const Vector &x, const double & t, Vector &v)
{
   switch(dim)
   {
      case 1:
      {
         v[0] = aq * pow(x[0],2) + a * x[0] + c;
         break;
      }
      case 2:
      {
         v[0] = aq * pow(x[0],2) + bq * pow(x[1],2) + a * x[0] + b * x[1] + c;
         v[1] = dq * pow(x[0],2) + eq * pow(x[1],2) + d * x[0] + e * x[1] + f;
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

   VectorFunctionCoefficient v_exact_coeff(dim, &velocity_exact);
   mv_gf.ProjectCoefficient(v_exact_coeff);
   mv_gf.SyncAliasMemory(S);

   // Initialize specific volume, velocity, and specific total energy
   ConstantCoefficient one_const_coeff(1.0);
   sv_gf.ProjectCoefficient(one_const_coeff);
   sv_gf.SyncAliasMemory(S);

   VectorConstantCoefficient one_coeff(one);
   v_gf.ProjectCoefficient(one_coeff);
   v_gf.SyncAliasMemory(S);

   ste_gf.ProjectCoefficient(one_const_coeff);
   ste_gf.SyncAliasMemory(S);

   // Just leave templated for hydro construction
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(one_const_coeff));
   m->Assemble();

   ProblemBase<dim> * problem_class = new SodProblem<dim>();

   mfem::hydroLO::LagrangianLOOperator<dim> hydro(H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);
   hydro.SetMVOption(mv_option);
   hydro.SetFVOption(fv_option);

   Vector _vel(dim), vec_res(dim);
   double t = 0., dt = .0000001;
   DenseMatrix dm(dim);

   hydro.CalculateTimestep(S);
   dt = hydro.GetTimestep();
   cout << "dt: " << dt << endl;

   hydro.ComputeMeshVelocities(S, S, t, dt); // For now just use S as S_old. Will need to replace
   // cout << "Done computing mesh velocities\n";

   // DenseMatrix A(3,2);
   // A(0,0) = 1.;
   // A(0,1) = 2.;
   // A(1,0) = 3.;
   // A(1,1) = 4.;
   // A(2,0) = 5.;
   // A(2,1) = 6.;

   // DenseMatrix AT = A;
   // AT.Transpose();

   // cout << "mat A:\n";
   // A.Print(cout);
   // cout << "mat AT:\n";
   // AT.Print(cout);

   // DenseMatrix ATA(2,2);
   // Mult(AT,A,ATA);
   // cout << "mat ATA:\n";
   // ATA.Print(cout);


   return 0;
}

/*
Purpose:
   To test the flux function in 1D to validate that the returned DenseMatrix is correct.

Input:
   dim = 2
   problem = 1
       | 2./5. |          | -2.  -3. |
   U = |  .3   | , F(U) = | .75   0  |
       |  .4   |          |  0   .75 |
       |  1.   | ,        | 2.25  3. | 
*/
int test_flux()
{
   cout << "test_flux\n";
   
   a = 0., b = 0., c = .3, d = 0., e = 0., f = .4;
   aq = 0., bq = 0., dq = 0., eq = 0.;
   
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

   VectorFunctionCoefficient v_exact_coeff(dim, &velocity_exact);
   mv_gf.ProjectCoefficient(v_exact_coeff);
   mv_gf.SyncAliasMemory(S);

   // Initialize specific volume, velocity, and specific total energy
   ConstantCoefficient one_const_coeff(1.0);
   sv_gf.ProjectCoefficient(one_const_coeff);
   sv_gf.SyncAliasMemory(S);

   VectorConstantCoefficient one_coeff(one);
   v_gf.ProjectCoefficient(one_coeff);
   v_gf.SyncAliasMemory(S);

   ste_gf.ProjectCoefficient(one_const_coeff);
   ste_gf.SyncAliasMemory(S);

   // Just leave templated for hydro construction
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(one_const_coeff));
   m->Assemble();

   ProblemBase<dim> * problem_class = new SodProblem<dim>();

   mfem::hydroLO::LagrangianLOOperator<dim> hydro(H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);
   hydro.SetMVOption(mv_option);
   hydro.SetFVOption(fv_option);

   double _error = 0.;

   Vector U(dim+2);
   U[0] = 2./5.;
   U[1] = .3;
   U[2] = .4;
   U[3] = 1.;

   DenseMatrix dm(dim+2, dim), dm_exact(dim+2, dim), dm_error(dim+2, dim);
   dm = problem_class->flux(U);
   dm_exact(0,0) = -.3, dm_exact(0,1) = -.4;
   dm_exact(1,0) = .875, dm_exact(1,1) = 0.;
   dm_exact(2,0) = 0., dm_exact(2,1) = .875;
   dm_exact(3,0) = .2625, dm_exact(3,1) = .35;

   Add(dm, dm_exact, -1., dm_error);
   _error += dm_error.FNorm();


   if (abs(_error) < tol)
   {
      return 0; // Test passed
   }
   else 
   {
      Vector test_vel(2);
      test_vel(0) = .3;
      test_vel(1) = .4;
      cout << "norm: " << test_vel.Norml2() << endl;
      cout << "pressure: " << problem_class->pressure(U) << endl;
      cout << "sie: " << problem_class->specific_internal_energy(U) << endl;
      cout << "For the flux we have:\n";
      dm.Print(cout);
      cout << "We should have:\n";
      dm_exact.Print(cout);
      cout << "error: " << _error << endl;
      cout << "failed flux check.\n";
      return 1; // Test failed
   }
}

/*
Purpose:
   To test Remark 5.3.
   This function tests a the RT computation of a continuous velocity field from a prescribed 
   velocity field that is linear in the first component and 0 in the second.

   v_exact = | 1 0 | |x| + |0|
             | 0 0 | |y|   |0|
*/
int test_vel_field_1()
{
   cout << "Test2D::test_vel_field_1\n";
   a = 1., b = 0., c = 0., d = 0., e = 0., f = 0.;

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
      mesh_name << "../results/mesh-Test2D-linear." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

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

   Vector zero(dim);
   zero = 0.;
   VectorConstantCoefficient zero_coeff(zero);
   mv_gf.ProjectCoefficient(zero_coeff);
   mv_gf.SyncAliasMemory(S);

   // Initialize specific volume, velocity, and specific total energy
   ConstantCoefficient one_const_coeff(1.0);
   sv_gf.ProjectCoefficient(one_const_coeff);
   sv_gf.SyncAliasMemory(S);

   v_gf.ProjectCoefficient(v_exact_coeff);
   v_gf.SyncAliasMemory(S);

   ste_gf.ProjectCoefficient(one_const_coeff);
   ste_gf.SyncAliasMemory(S);

   // Just leave templated for hydro construction
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(one_const_coeff));
   m->Assemble();

   ProblemBase<dim> * problem_class = new SodProblem<dim>();

   mfem::hydroLO::LagrangianLOOperator<dim> hydro(H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   cout << "Done constructing hydro op\n";

   Vector _vel(dim), vec_res(dim);
   double t = 0., dt = .0000001;
   DenseMatrix dm(dim);

   hydro.CalculateTimestep(S);
   dt = hydro.GetTimestep();
   cout << "dt: " << dt << endl;

   // hydro.ComputeMeshVelocities(S, t, dt, flag, &velocity_exact);
   hydro.ComputeMeshVelocities(S, S, t, dt); // For now just use S as S_old. Will need to replace
   cout << "Done computing mesh velocities\n";

   /* ************************
   Displace Velocities
   *************************** */ 
   socketstream vis_vh, vis_vexact, vis_vcrgf, vis_ifv_exact;
   char vishost[] = "localhost";
   int visport = 19916;

   MPI_Barrier(pmesh->GetComm());

   int Wx = 0, Wy = 0;
   const int Ww = 350, Wh = 350;
   int offx = Ww+10, offy = Wh+45;

   hydroLO::VisualizeField(vis_vh, vishost, visport, mv_gf, "2D: Reconstructed Mesh Velocity", Wx, Wy, Ww, Wh);

   Wx += offx;

   hydroLO::VisualizeField(vis_vexact, vishost, visport, v_exact_gf, "2D: Exact Velocity", Wx, Wy, Ww, Wh);

   Wx += offx;
   ParGridFunction v_CR_gf(&CRFESpace);
   hydro.GetVCRgf(v_CR_gf);

   hydroLO::VisualizeField(vis_vcrgf, vishost, visport, v_CR_gf, "2D: v_CR_gf", Wx, Wy, Ww, Wh);

   Wx += offx;

   ParGridFunction ifv_exact(&CRFESpace);
   ifv_exact.ProjectCoefficient(v_exact_coeff);
   hydroLO::VisualizeField(vis_ifv_exact, vishost, visport, ifv_exact, "2D: ifv_exact", Wx, Wy, Ww, Wh);


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
      mesh_name << "../results/mesh-Test2D-linear-moved." << setfill('0') << setw(6) << myid;
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
      cout << "error: " << _error << endl;
      cout << "passed Test2D::test_vel_field_1\n";
      return 0; // Test passed
   }
   else 
   {
      cout << "error: " << _error << endl;
      cout << "Failed test_vel_field_1\n";
      cout << "error gf:\n";
      vel_err.Print(cout);
      return 1; // Test failed
   }
}

int test_CSV_getter_setter()
{
   a = 1., b = 1., c = 1., d = 1., e = -1., f = 1.;

   // Initialize mesh
   mfem::Mesh *mesh;
   mesh = new mfem::Mesh(mesh_file, true, true);

   // Refine the mesh
   for (int lev = 0; lev < 3; lev++)
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
      mesh_name << "../results/mesh-Test2D-CSV." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

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

   VectorFunctionCoefficient v_exact_coeff(dim, &velocity_exact);

   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   mv_gf.MakeRef(&H1FESpace, S, offset[1]);
   sv_gf.MakeRef(&L2FESpace, S, offset[2]);
   v_gf.MakeRef(&L2VFESpace, S, offset[3]);
   ste_gf.MakeRef(&L2FESpace, S, offset[4]);

   pmesh->SetNodalGridFunction(&x_gf);
   x_gf.SyncAliasMemory(S);

   Vector n_one(dim);
   n_one = -1.;
   VectorConstantCoefficient n_one_coeff(n_one);
   mv_gf.ProjectCoefficient(n_one_coeff);
   mv_gf.SyncAliasMemory(S);

   // Initialize specific volume, velocity, and specific total energy
   ConstantCoefficient one_const_coeff(1.0);
   sv_gf.ProjectCoefficient(one_const_coeff);
   sv_gf.SyncAliasMemory(S);

   Vector vel(dim);
   vel[0] = 1.1;
   vel[1] = 1.2;
   VectorConstantCoefficient vel_coeff(vel);
   v_gf.ProjectCoefficient(vel_coeff);
   v_gf.SyncAliasMemory(S);

   ConstantCoefficient one_three_coeff(1.3);
   ste_gf.ProjectCoefficient(one_three_coeff);
   ste_gf.SyncAliasMemory(S);

   // Just leave templated for hydro construction
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(one_const_coeff));
   m->Assemble();

   ProblemBase<dim> * problem_class = new SodProblem<dim>();

   mfem::hydroLO::LagrangianLOOperator<dim> hydro(H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   cout << "S:\n";
   S.Print(cout);
   Vector _vel(dim), vec_res(dim);
   double t = 0., dt = 1.;

   int cell = 5;
   Vector res(dim+2), new_res(dim+2);
   new_res = 0.;
   res[0] = 0.;
   res[1] = 10.;
   res[2] = 20.;
   res[3] = 30.;
   hydro.SetCellStateVector(S, cell, res);

   cout << "S:\n";
   S.Print(cout);

   hydro.GetCellStateVector(S, cell, new_res);

   new_res.Add(-1., res);
   double _error = new_res.Norml2();

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


/***
 * Purpose:
 * Following a discussion with Valdimir on 2023-08-16, I should be checking the hydrodynamics
 * between the 1D and the 2D run of the sod problem, as well as the mesh velocities.  The corner 
 * nodes in 2D should be equivalent to their face counterparts in 1D, otherwise the 2d mesh
 * velocity calculation will require face correction to be locally mass conservative, which is
 * something that should not happen.
 * 
 * For validation, I simply looked at the meshes and check that the cells by column on the
 * 2D mesh are the same as their respective column on the 1d mesh.  The indices for columns
 * is hardcoded.
*/
int test_sod_hydro()
{
   cout << "Test2D::test_sod_hydro\n";
   int _mesh_refinements = 2; // This is enough refinements to check, don't modify or test will crash
   int ms = 10;
   double t = 0., dt_1d = .001, dt_2d = 0.001;

   /*************************************
    * RUN 1D first
   **************************************/

   // Initialize 1D mesh
   const int dim_1d = 1;
   std::string result_1d = std::string(LAGLOS_DIR) + "/data/ref-segment.mesh";
   const char* mesh_file_1d = result_1d.c_str();
   mfem::Mesh *mesh_1d = new mfem::Mesh(mesh_file_1d, true, true);

   // Refine the mesh
   for (int lev = 0; lev < _mesh_refinements; lev++)
   {
      mesh_1d->UniformRefinement();
   }

   // Construct parmesh object
   ParMesh *pmesh_1d = new ParMesh(MPI_COMM_WORLD, *mesh_1d);            
   delete mesh_1d;

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-Test2D-sod_hydro_1d." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh_1d->Print(omesh);
   }

   // Template FE stuff to construct hydro operator
   H1_FECollection H1FEC_1d(order_mv, dim_1d);
   H1_FECollection H1FEC_1d_L(1, dim_1d);
   L2_FECollection L2FEC_1d(order_u, dim_1d, BasisType::Positive);
   CrouzeixRaviartFECollection CRFEC_1d;

   ParFiniteElementSpace H1FESpace_1d(pmesh_1d, &H1FEC_1d, dim_1d);
   ParFiniteElementSpace H1FESpace_1d_L(pmesh_1d, &H1FEC_1d_L, dim_1d);
   ParFiniteElementSpace L2FESpace_1d(pmesh_1d, &L2FEC_1d);
   ParFiniteElementSpace L2VFESpace_1d(pmesh_1d, &L2FEC_1d, dim_1d);
   ParFiniteElementSpace CRFESpace_1d(pmesh_1d, &CRFEC_1d, dim_1d);

   // Output information
   pmesh_1d->ExchangeFaceNbrData();

   // Construct blockvector
   const int Vsize_l2_1d = L2FESpace_1d.GetVSize();
   const int Vsize_l2v_1d = L2VFESpace_1d.GetVSize();
   const int Vsize_h1_1d = H1FESpace_1d.GetVSize();
   Array<int> offset_1d(6);
   offset_1d[0] = 0;
   offset_1d[1] = offset_1d[0] + Vsize_h1_1d;
   offset_1d[2] = offset_1d[1] + Vsize_h1_1d;
   offset_1d[3] = offset_1d[2] + Vsize_l2_1d;
   offset_1d[4] = offset_1d[3] + Vsize_l2v_1d;
   offset_1d[5] = offset_1d[4] + Vsize_l2_1d;
   BlockVector S_1d(offset_1d, Device::GetMemoryType());

   // Pair with corresponding gridfunctions
   // Only need position and velocity for testing
   ParGridFunction x_gf_1d, mv_gf_1d, sv_gf_1d, v_gf_1d, ste_gf_1d;
   x_gf_1d.MakeRef(&H1FESpace_1d, S_1d, offset_1d[0]);
   mv_gf_1d.MakeRef(&H1FESpace_1d, S_1d, offset_1d[1]);
   sv_gf_1d.MakeRef(&L2FESpace_1d, S_1d, offset_1d[2]);
   v_gf_1d.MakeRef(&L2VFESpace_1d, S_1d, offset_1d[3]);
   ste_gf_1d.MakeRef(&L2FESpace_1d, S_1d, offset_1d[4]);

   Vector zero_1d(dim_1d);
   zero_1d = 0.;
   VectorConstantCoefficient zero_coeff_1d(zero_1d);
   mv_gf_1d.ProjectCoefficient(zero_coeff_1d);
   mv_gf_1d.SyncAliasMemory(S_1d);

   pmesh_1d->SetNodalGridFunction(&x_gf_1d);
   x_gf_1d.SyncAliasMemory(S_1d);

   // Initialize specific volume, velocity, and specific total energy
   // Change class variables into static std::functions since virtual static member functions are not an option
   // and Coefficient class requires std::function arguments
   ProblemBase<dim_1d> * problem_class_1d = new SodProblem<dim_1d>();
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> sv0_static_1d = 
      std::bind(&ProblemBase<dim_1d>::sv0, problem_class_1d, std::placeholders::_1, std::placeholders::_2);
   std::function<void(const Vector &, const double, Vector &)> v0_static_1d = 
      std::bind(&ProblemBase<dim_1d>::v0, problem_class_1d, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
   std::function<double(const Vector &,const double)> ste0_static_1d = 
      std::bind(&ProblemBase<dim_1d>::ste0, problem_class_1d, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> rho0_static_1d = 
      std::bind(&ProblemBase<dim_1d>::rho0, problem_class_1d, std::placeholders::_1, std::placeholders::_2);

   VectorFunctionCoefficient v_exact_coeff_1d(dim_1d, v0_static_1d);

   // Initialize specific volume, velocity, and specific total energy
   FunctionCoefficient sv_coeff_1d(sv0_static_1d);
   sv_coeff_1d.SetTime(t);
   sv_gf_1d.ProjectCoefficient(sv_coeff_1d);
   sv_gf_1d.SyncAliasMemory(S_1d);

   v_exact_coeff_1d.SetTime(t);
   v_gf_1d.ProjectCoefficient(v_exact_coeff_1d);
   // Sync the data location of v_gf with its base, S
   v_gf_1d.SyncAliasMemory(S_1d);

   FunctionCoefficient ste_coeff_1d(ste0_static_1d);
   ste_coeff_1d.SetTime(t);
   ste_gf_1d.ProjectCoefficient(ste_coeff_1d);
   ste_gf_1d.SyncAliasMemory(S_1d);

   // PLF to build mass vector
   FunctionCoefficient rho_coeff_1d(rho0_static_1d); 
   rho_coeff_1d.SetTime(t);
   ParLinearForm *m_1d = new ParLinearForm(&L2FESpace_1d);
   m_1d->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff_1d));
   m_1d->Assemble();

   cout << "gridFunctions initiated.\n";

   mfem::hydroLO::LagrangianLOOperator<dim_1d> hydro_1d(H1FESpace_1d, H1FESpace_1d_L, L2FESpace_1d, L2VFESpace_1d, CRFESpace_1d, m_1d, problem_class_1d, use_viscosity, _mm, CFL);

   cout << "Done constructing hydro op\n";

   bool is_collapsed = false;
   for (int ts = 1; ts < ms; ts++)
   {
      // Verify hardcoded timestep is small enough
      hydro_1d.BuildDijMatrix(S_1d);
      hydro_1d.CalculateTimestep(S_1d);
      if (dt_1d > hydro_1d.GetTimestep()) 
      {
         MFEM_ABORT("Sod hydro :: 1D Too big timestep\n");
      }

      hydro_1d.MakeTimeStep(S_1d, t, dt_1d, is_collapsed);

      // Ensure the sub-vectors x_gf, v_gf, and ste_gf know the location of the
      // data in S. This operation simply updates the Memory validity flags of
      // the sub-vectors to match those of S.
      x_gf_1d.SyncAliasMemory(S_1d);
      mv_gf_1d.SyncAliasMemory(S_1d);
      sv_gf_1d.SyncAliasMemory(S_1d);
      v_gf_1d.SyncAliasMemory(S_1d);
      ste_gf_1d.SyncAliasMemory(S_1d);

      // Make sure that the mesh corresponds to the new solution state. This is
      // needed, because some time integrators use different S-type vectors
      // and the oper object might have redirected the mesh positions to those.
      pmesh_1d->NewNodes(x_gf_1d, false);

      ParGridFunction mc_gf_1d(&L2FESpace_1d);
      hydro_1d.CheckMassConservation(S_1d, mc_gf_1d);
   }

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-Test2D-sod_hydro_1d-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh_1d->Print(omesh);
   }


   /*************************************
    * Run 2D now
   **************************************/
   const int dim_2d = 2;
   std::string result_2d = std::string(LAGLOS_DIR) + "/data/ref-square.mesh";
   const char* mesh_file_2d = result_2d.c_str();
   mfem::Mesh *mesh_2d = new mfem::Mesh(mesh_file_2d, true, true);

   // Refine the mesh
   for (int lev = 0; lev < _mesh_refinements; lev++)
   {
      mesh_2d->UniformRefinement();
   }

   // Construct parmesh object
   ParMesh *pmesh_2d = new ParMesh(MPI_COMM_WORLD, *mesh_2d);            
   delete mesh_2d;

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-Test2D-sod_hydro_2d." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh_2d->Print(omesh);
   }

   // Template FE stuff to construct hydro operator
   H1_FECollection H1FEC_2d(order_mv, dim_2d);
   H1_FECollection H1FEC_2d_L(1, dim_2d);
   L2_FECollection L2FEC_2d(order_u, dim_2d, BasisType::Positive);
   CrouzeixRaviartFECollection CRFEC_2d;

   ParFiniteElementSpace H1FESpace_2d(pmesh_2d, &H1FEC_2d, dim_2d);
   ParFiniteElementSpace H1FESpace_2d_L(pmesh_2d, &H1FEC_2d_L, dim_2d);
   ParFiniteElementSpace L2FESpace_2d(pmesh_2d, &L2FEC_2d);
   ParFiniteElementSpace L2VFESpace_2d(pmesh_2d, &L2FEC_2d, dim_2d);
   ParFiniteElementSpace CRFESpace_2d(pmesh_2d, &CRFEC_2d, dim_2d);

   // Output information
   pmesh_2d->ExchangeFaceNbrData();

   // Construct blockvector
   const int Vsize_l2_2d = L2FESpace_2d.GetVSize();
   const int Vsize_l2v_2d = L2VFESpace_2d.GetVSize();
   const int Vsize_h1_2d = H1FESpace_2d.GetVSize();
   Array<int> offset_2d(6);
   offset_2d[0] = 0;
   offset_2d[1] = offset_2d[0] + Vsize_h1_2d;
   offset_2d[2] = offset_2d[1] + Vsize_h1_2d;
   offset_2d[3] = offset_2d[2] + Vsize_l2_2d;
   offset_2d[4] = offset_2d[3] + Vsize_l2v_2d;
   offset_2d[5] = offset_2d[4] + Vsize_l2_2d;
   BlockVector S_2d(offset_2d, Device::GetMemoryType());

   // Pair with corresponding gridfunctions
   // Only need position and velocity for testing
   ParGridFunction x_gf_2d, mv_gf_2d, sv_gf_2d, v_gf_2d, ste_gf_2d;
   x_gf_2d.MakeRef(&H1FESpace_2d, S_2d, offset_2d[0]);
   mv_gf_2d.MakeRef(&H1FESpace_2d, S_2d, offset_2d[1]);
   sv_gf_2d.MakeRef(&L2FESpace_2d, S_2d, offset_2d[2]);
   v_gf_2d.MakeRef(&L2VFESpace_2d, S_2d, offset_2d[3]);
   ste_gf_2d.MakeRef(&L2FESpace_2d, S_2d, offset_2d[4]);

   Vector zero_2d(dim_2d);
   zero_2d = 0.;
   VectorConstantCoefficient zero_coeff_2d(zero_2d);
   mv_gf_2d.ProjectCoefficient(zero_coeff_2d);
   mv_gf_2d.SyncAliasMemory(S_2d);

   pmesh_2d->SetNodalGridFunction(&x_gf_2d);
   x_gf_2d.SyncAliasMemory(S_2d);

   t = 0.;

   // Initialize specific volume, velocity, and specific total energy
   // Change class variables into static std::functions since virtual static member functions are not an option
   // and Coefficient class requires std::function arguments
   ProblemBase<dim_2d> * problem_class_2d = new SodProblem<dim_2d>();
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> sv0_static_2d = 
      std::bind(&ProblemBase<dim_2d>::sv0, problem_class_2d, std::placeholders::_1, std::placeholders::_2);
   std::function<void(const Vector &, const double, Vector &)> v0_static_2d = 
      std::bind(&ProblemBase<dim_2d>::v0, problem_class_2d, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
   std::function<double(const Vector &,const double)> ste0_static_2d = 
      std::bind(&ProblemBase<dim_2d>::ste0, problem_class_2d, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> rho0_static_2d = 
      std::bind(&ProblemBase<dim_2d>::rho0, problem_class_2d, std::placeholders::_1, std::placeholders::_2);

   VectorFunctionCoefficient v_exact_coeff_2d(dim_2d, v0_static_2d);

   // Initialize specific volume, velocity, and specific total energy
   FunctionCoefficient sv_coeff_2d(sv0_static_2d);
   sv_coeff_2d.SetTime(t);
   sv_gf_2d.ProjectCoefficient(sv_coeff_2d);
   sv_gf_2d.SyncAliasMemory(S_2d);

   v_exact_coeff_2d.SetTime(t);
   v_gf_2d.ProjectCoefficient(v_exact_coeff_2d);
   // Sync the data location of v_gf with its base, S
   v_gf_2d.SyncAliasMemory(S_2d);

   FunctionCoefficient ste_coeff_2d(ste0_static_2d);
   ste_coeff_2d.SetTime(t);
   ste_gf_2d.ProjectCoefficient(ste_coeff_2d);
   ste_gf_2d.SyncAliasMemory(S_2d);

   // PLF to build mass vector
   FunctionCoefficient rho_coeff_2d(rho0_static_2d); 
   rho_coeff_2d.SetTime(t);
   ParLinearForm *m_2d = new ParLinearForm(&L2FESpace_2d);
   m_2d->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff_2d));
   m_2d->Assemble();

   cout << "gridFunctions initiated.\n";

   mfem::hydroLO::LagrangianLOOperator<dim_2d> hydro_2d(H1FESpace_2d, H1FESpace_2d_L, L2FESpace_2d, L2VFESpace_2d, CRFESpace_2d, m_2d, problem_class_2d, use_viscosity, _mm, CFL);

   cout << "Done constructing hydro op\n";

   is_collapsed = false;

   for (int ts = 1; ts < ms; ts++)
   {
      // Verify hardcoded timestep is small enough
      hydro_2d.BuildDijMatrix(S_2d);
      hydro_2d.CalculateTimestep(S_2d);
      if (dt_2d > hydro_2d.GetTimestep()) 
      {
         MFEM_ABORT("Sod hydro :: 1D Too big timestep\n");
      }

      hydro_2d.MakeTimeStep(S_2d, t, dt_2d, is_collapsed);

      // Ensure the sub-vectors x_gf, v_gf, and ste_gf know the location of the
      // data in S. This operation simply updates the Memory validity flags of
      // the sub-vectors to match those of S.
      x_gf_2d.SyncAliasMemory(S_2d);
      mv_gf_2d.SyncAliasMemory(S_2d);
      sv_gf_2d.SyncAliasMemory(S_2d);
      v_gf_2d.SyncAliasMemory(S_2d);
      ste_gf_2d.SyncAliasMemory(S_2d);

      // Make sure that the mesh corresponds to the new solution state. This is
      // needed, because some time integrators use different S-type vectors
      // and the oper object might have redirected the mesh positions to those.
      pmesh_2d->NewNodes(x_gf_2d, false);

      ParGridFunction mc_gf_2d(&L2FESpace_2d);
      hydro_2d.CheckMassConservation(S_2d, mc_gf_2d);
   }

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-Test2D-sod_hydro_2d-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh_2d->Print(omesh);
   }

   cout << "mv_1d:\n";
   mv_gf_1d.Print(cout);

   cout << "mv_2d:\n";
   mv_gf_2d.Print(cout);

   // hydro verification
   // Verify that columns of 2d sim match columns of 1d sim
   int a0[4] = {0, 3, 12, 15};
   int a1[4] = {4, 7, 8, 11};
   int a2[4] = {1, 2, 13, 14};
   int a3[4] = {5, 6, 9, 10};

   double _sv0 = sv_gf_1d[0], _sv1 = sv_gf_1d[1], _sv2 = sv_gf_1d[2], _sv3 = sv_gf_1d[3];
   double _v0 = v_gf_1d[0], _v1 = v_gf_1d[1], _v2 = v_gf_1d[2], _v3 = v_gf_1d[3];
   double _ste0 = ste_gf_1d[0], _ste1 = ste_gf_1d[1], _ste2 = ste_gf_1d[2], _ste3 = ste_gf_1d[3];

   for (int i = 0; i < 4; i++)
   {
      // check sv
      if (_sv0 != sv_gf_2d[a0[i]] || _sv1 != sv_gf_2d[a1[i]] ||
          _sv2 != sv_gf_2d[a2[i]] || _sv3 != sv_gf_2d[a3[i]])
      {
         cout << "Test2D::test_sod_hydro failed on sv.\n";
         return 1;
      }
      
      // check v
      if (abs(_v0 - v_gf_2d[a0[i]]) > tol || abs(_v1 - v_gf_2d[a1[i]]) > tol ||
          abs(_v2 - v_gf_2d[a2[i]]) > tol || abs(_v3 - v_gf_2d[a3[i]]) > tol) 
      {
         cout << "Test2D::test_sod_hydro failed on v.\n"; 

         cout << "i: " << i << endl;
         cout << _v0 << ", " << v_gf_2d[a0[i]] << endl;
         cout << _v1 << ", " << v_gf_2d[a1[i]] << endl;
         cout << _v2 << ", " << v_gf_2d[a2[i]] << endl;
         cout << _v3 << ", " << v_gf_2d[a3[i]] << endl;

         return 1; 
      }

      // check ste
      if (_ste0 != ste_gf_2d[a0[i]] || _ste1 != ste_gf_2d[a1[i]] ||
          _ste2 != ste_gf_2d[a2[i]] || _ste3 != ste_gf_2d[a3[i]]) 
      {
         cout << "Test2D::test_sod_hydro failed on ste.\n"; 
         return 1; 
      }
   }

   // mesh velocity verification
   // Verify that columns of 2d sim match columns of 1d sim
   int b0[5] = {0, 3, 7, 12, 20};
   int b1[5] = {1, 2, 5, 14, 16};
   int b2[5] = {4, 6, 8, 10, 18};
   int b3[5] = {9, 11, 19, 21, 21};
   int b4[5] = {13, 15, 17, 22, 23};
   // double _sv0 = sv_gf_1d[0], _sv1 = sv_gf_1d[1], _sv2 = sv_gf_1d[2], _sv3 = sv_gf_1d[3];
   double _mv0 = mv_gf_1d[0], _mv1 = mv_gf_1d[1], _mv2 = mv_gf_1d[2], _mv3 = mv_gf_1d[3], _mv4 = mv_gf_1d[4];

   for (int i = 0; i < 5; i++)
   {
      // check v
      if (abs(_mv0 - mv_gf_2d[b0[i]]) > tol || abs(_mv1 - mv_gf_2d[b1[i]]) > tol ||
          abs(_mv2 - mv_gf_2d[b2[i]]) > tol || abs(_mv3 - mv_gf_2d[b3[i]]) > tol ||
          abs(_mv4 - mv_gf_2d[b4[i]]) > tol) 
      {
         cout << "Test2D::test_sod_hydro failed on mv.\n"; 

         cout << "i: " << i << endl;
         cout << _mv0 << ", " << mv_gf_2d[b0[i]] << endl;
         cout << _mv1 << ", " << mv_gf_2d[b1[i]] << endl;
         cout << _mv2 << ", " << mv_gf_2d[b2[i]] << endl;
         cout << _mv3 << ", " << mv_gf_2d[b3[i]] << endl;
         cout << _mv4 << ", " << mv_gf_2d[b4[i]] << endl;

         return 1; 
      }
   }

   return 0;
}

/**
 * This function provides the same functionality as the above Sod problem, but with the smooth 
 * wave.
*/
int test_smooth_hydro()
{
   cout << "Test2D::test_smooth_hydro\n";
   int _mesh_refinements = 2; // This is enough refinements to check, don't modify or test will crash
   double t = 0., dt_1d = .001, dt_2d = 0.001;

   /*************************************
    * RUN 1D first
   **************************************/

   // Initialize 1D mesh
   const int dim_1d = 1;
   std::string result_1d = std::string(LAGLOS_DIR) + "/data/ref-segment.mesh";
   const char* mesh_file_1d = result_1d.c_str();
   mfem::Mesh *mesh_1d = new mfem::Mesh(mesh_file_1d, true, true);

   // Refine the mesh
   for (int lev = 0; lev < _mesh_refinements; lev++)
   {
      mesh_1d->UniformRefinement();
   }

   // Construct parmesh object
   ParMesh *pmesh_1d = new ParMesh(MPI_COMM_WORLD, *mesh_1d);            
   delete mesh_1d;

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-Test2D-smooth_1d." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh_1d->Print(omesh);
   }

   // Template FE stuff to construct hydro operator
   H1_FECollection H1FEC_1d(order_mv, dim_1d);
   H1_FECollection H1FEC_1d_L(1, dim_1d);
   L2_FECollection L2FEC_1d(order_u, dim_1d, BasisType::Positive);
   CrouzeixRaviartFECollection CRFEC_1d;

   ParFiniteElementSpace H1FESpace_1d(pmesh_1d, &H1FEC_1d, dim_1d);
   ParFiniteElementSpace H1FESpace_1d_L(pmesh_1d, &H1FEC_1d_L, dim_1d);
   ParFiniteElementSpace L2FESpace_1d(pmesh_1d, &L2FEC_1d);
   ParFiniteElementSpace L2VFESpace_1d(pmesh_1d, &L2FEC_1d, dim_1d);
   ParFiniteElementSpace CRFESpace_1d(pmesh_1d, &CRFEC_1d, dim_1d);

   // Output information
   pmesh_1d->ExchangeFaceNbrData();

   // Construct blockvector
   const int Vsize_l2_1d = L2FESpace_1d.GetVSize();
   const int Vsize_l2v_1d = L2VFESpace_1d.GetVSize();
   const int Vsize_h1_1d = H1FESpace_1d.GetVSize();
   Array<int> offset_1d(6);
   offset_1d[0] = 0;
   offset_1d[1] = offset_1d[0] + Vsize_h1_1d;
   offset_1d[2] = offset_1d[1] + Vsize_h1_1d;
   offset_1d[3] = offset_1d[2] + Vsize_l2_1d;
   offset_1d[4] = offset_1d[3] + Vsize_l2v_1d;
   offset_1d[5] = offset_1d[4] + Vsize_l2_1d;
   BlockVector S_1d(offset_1d, Device::GetMemoryType());

   // Pair with corresponding gridfunctions
   // Only need position and velocity for testing
   ParGridFunction x_gf_1d, mv_gf_1d, sv_gf_1d, v_gf_1d, ste_gf_1d;
   x_gf_1d.MakeRef(&H1FESpace_1d, S_1d, offset_1d[0]);
   mv_gf_1d.MakeRef(&H1FESpace_1d, S_1d, offset_1d[1]);
   sv_gf_1d.MakeRef(&L2FESpace_1d, S_1d, offset_1d[2]);
   v_gf_1d.MakeRef(&L2VFESpace_1d, S_1d, offset_1d[3]);
   ste_gf_1d.MakeRef(&L2FESpace_1d, S_1d, offset_1d[4]);

   Vector zero_1d(dim_1d);
   zero_1d = 0.;
   VectorConstantCoefficient zero_coeff_1d(zero_1d);
   mv_gf_1d.ProjectCoefficient(zero_coeff_1d);
   mv_gf_1d.SyncAliasMemory(S_1d);

   pmesh_1d->SetNodalGridFunction(&x_gf_1d);
   x_gf_1d.SyncAliasMemory(S_1d);

   // Initialize specific volume, velocity, and specific total energy
   // Change class variables into static std::functions since virtual static member functions are not an option
   // and Coefficient class requires std::function arguments
   ProblemBase<dim_1d> * problem_class_1d = new SmoothWave<dim_1d>();
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> sv0_static_1d = 
      std::bind(&ProblemBase<dim_1d>::sv0, problem_class_1d, std::placeholders::_1, std::placeholders::_2);
   std::function<void(const Vector &, const double, Vector &)> v0_static_1d = 
      std::bind(&ProblemBase<dim_1d>::v0, problem_class_1d, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
   std::function<double(const Vector &,const double)> ste0_static_1d = 
      std::bind(&ProblemBase<dim_1d>::ste0, problem_class_1d, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> rho0_static_1d = 
      std::bind(&ProblemBase<dim_1d>::rho0, problem_class_1d, std::placeholders::_1, std::placeholders::_2);

   VectorFunctionCoefficient v_exact_coeff_1d(dim_1d, v0_static_1d);

   // Initialize specific volume, velocity, and specific total energy
   FunctionCoefficient sv_coeff_1d(sv0_static_1d);
   sv_coeff_1d.SetTime(t);
   sv_gf_1d.ProjectCoefficient(sv_coeff_1d);
   sv_gf_1d.SyncAliasMemory(S_1d);

   v_exact_coeff_1d.SetTime(t);
   v_gf_1d.ProjectCoefficient(v_exact_coeff_1d);
   // Sync the data location of v_gf with its base, S
   v_gf_1d.SyncAliasMemory(S_1d);

   FunctionCoefficient ste_coeff_1d(ste0_static_1d);
   ste_coeff_1d.SetTime(t);
   ste_gf_1d.ProjectCoefficient(ste_coeff_1d);
   ste_gf_1d.SyncAliasMemory(S_1d);

   // PLF to build mass vector
   FunctionCoefficient rho_coeff_1d(rho0_static_1d); 
   rho_coeff_1d.SetTime(t);
   ParLinearForm *m_1d = new ParLinearForm(&L2FESpace_1d);
   m_1d->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff_1d));
   m_1d->Assemble();

   cout << "gridFunctions initiated.\n";

   mfem::hydroLO::LagrangianLOOperator<dim_1d> hydro_1d(H1FESpace_1d, H1FESpace_1d_L, L2FESpace_1d, L2VFESpace_1d, CRFESpace_1d, m_1d, problem_class_1d, use_viscosity, _mm, CFL);

   cout << "Done constructing hydro op\n";

   // Verify hardcoded timestep is small enough
   hydro_1d.BuildDijMatrix(S_1d);
   hydro_1d.CalculateTimestep(S_1d);
   if (dt_1d > hydro_1d.GetTimestep()) 
   {
      MFEM_ABORT("Sod hydro :: 1D Too big timestep\n");
   }

   bool is_collapsed = false;
   hydro_1d.MakeTimeStep(S_1d, t, dt_1d, is_collapsed);

   // Ensure the sub-vectors x_gf, v_gf, and ste_gf know the location of the
   // data in S. This operation simply updates the Memory validity flags of
   // the sub-vectors to match those of S.
   x_gf_1d.SyncAliasMemory(S_1d);
   mv_gf_1d.SyncAliasMemory(S_1d);
   sv_gf_1d.SyncAliasMemory(S_1d);
   v_gf_1d.SyncAliasMemory(S_1d);
   ste_gf_1d.SyncAliasMemory(S_1d);

   // Make sure that the mesh corresponds to the new solution state. This is
   // needed, because some time integrators use different S-type vectors
   // and the oper object might have redirected the mesh positions to those.
   pmesh_1d->NewNodes(x_gf_1d, false);

   ParGridFunction mc_gf_1d(&L2FESpace_1d);
   hydro_1d.CheckMassConservation(S_1d, mc_gf_1d);

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-Test2D-smooth_1d-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh_1d->Print(omesh);
   }


   /*************************************
    * Run 2D now
   **************************************/
   const int dim_2d = 2;
   std::string result_2d = std::string(LAGLOS_DIR) + "/data/ref-square.mesh";
   const char* mesh_file_2d = result_2d.c_str();
   mfem::Mesh *mesh_2d = new mfem::Mesh(mesh_file_2d, true, true);

   // Refine the mesh
   for (int lev = 0; lev < _mesh_refinements; lev++)
   {
      mesh_2d->UniformRefinement();
   }

   // Construct parmesh object
   ParMesh *pmesh_2d = new ParMesh(MPI_COMM_WORLD, *mesh_2d);            
   delete mesh_2d;

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-Test2D-smooth_2d." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh_2d->Print(omesh);
   }

   // Template FE stuff to construct hydro operator
   H1_FECollection H1FEC_2d(order_mv, dim_2d);
   H1_FECollection H1FEC_2d_L(1, dim_2d);
   L2_FECollection L2FEC_2d(order_u, dim_2d, BasisType::Positive);
   CrouzeixRaviartFECollection CRFEC_2d;

   ParFiniteElementSpace H1FESpace_2d(pmesh_2d, &H1FEC_2d, dim_2d);
   ParFiniteElementSpace H1FESpace_2d_L(pmesh_2d, &H1FEC_2d_L, dim_2d);
   ParFiniteElementSpace L2FESpace_2d(pmesh_2d, &L2FEC_2d);
   ParFiniteElementSpace L2VFESpace_2d(pmesh_2d, &L2FEC_2d, dim_2d);
   ParFiniteElementSpace CRFESpace_2d(pmesh_2d, &CRFEC_2d, dim_2d);

   // Output information
   pmesh_2d->ExchangeFaceNbrData();

   // Construct blockvector
   const int Vsize_l2_2d = L2FESpace_2d.GetVSize();
   const int Vsize_l2v_2d = L2VFESpace_2d.GetVSize();
   const int Vsize_h1_2d = H1FESpace_2d.GetVSize();
   Array<int> offset_2d(6);
   offset_2d[0] = 0;
   offset_2d[1] = offset_2d[0] + Vsize_h1_2d;
   offset_2d[2] = offset_2d[1] + Vsize_h1_2d;
   offset_2d[3] = offset_2d[2] + Vsize_l2_2d;
   offset_2d[4] = offset_2d[3] + Vsize_l2v_2d;
   offset_2d[5] = offset_2d[4] + Vsize_l2_2d;
   BlockVector S_2d(offset_2d, Device::GetMemoryType());

   // Pair with corresponding gridfunctions
   // Only need position and velocity for testing
   ParGridFunction x_gf_2d, mv_gf_2d, sv_gf_2d, v_gf_2d, ste_gf_2d;
   x_gf_2d.MakeRef(&H1FESpace_2d, S_2d, offset_2d[0]);
   mv_gf_2d.MakeRef(&H1FESpace_2d, S_2d, offset_2d[1]);
   sv_gf_2d.MakeRef(&L2FESpace_2d, S_2d, offset_2d[2]);
   v_gf_2d.MakeRef(&L2VFESpace_2d, S_2d, offset_2d[3]);
   ste_gf_2d.MakeRef(&L2FESpace_2d, S_2d, offset_2d[4]);

   Vector zero_2d(dim_2d);
   zero_2d = 0.;
   VectorConstantCoefficient zero_coeff_2d(zero_2d);
   mv_gf_2d.ProjectCoefficient(zero_coeff_2d);
   mv_gf_2d.SyncAliasMemory(S_2d);

   pmesh_2d->SetNodalGridFunction(&x_gf_2d);
   x_gf_2d.SyncAliasMemory(S_2d);

   t = 0.;

   // Initialize specific volume, velocity, and specific total energy
   // Change class variables into static std::functions since virtual static member functions are not an option
   // and Coefficient class requires std::function arguments
   ProblemBase<dim_2d> * problem_class_2d = new SmoothWave<dim_2d>();
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> sv0_static_2d = 
      std::bind(&ProblemBase<dim_2d>::sv0, problem_class_2d, std::placeholders::_1, std::placeholders::_2);
   std::function<void(const Vector &, const double, Vector &)> v0_static_2d = 
      std::bind(&ProblemBase<dim_2d>::v0, problem_class_2d, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
   std::function<double(const Vector &,const double)> ste0_static_2d = 
      std::bind(&ProblemBase<dim_2d>::ste0, problem_class_2d, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> rho0_static_2d = 
      std::bind(&ProblemBase<dim_2d>::rho0, problem_class_2d, std::placeholders::_1, std::placeholders::_2);

   VectorFunctionCoefficient v_exact_coeff_2d(dim_2d, v0_static_2d);

   // Initialize specific volume, velocity, and specific total energy
   FunctionCoefficient sv_coeff_2d(sv0_static_2d);
   sv_coeff_2d.SetTime(t);
   sv_gf_2d.ProjectCoefficient(sv_coeff_2d);
   sv_gf_2d.SyncAliasMemory(S_2d);

   v_exact_coeff_2d.SetTime(t);
   v_gf_2d.ProjectCoefficient(v_exact_coeff_2d);
   // Sync the data location of v_gf with its base, S
   v_gf_2d.SyncAliasMemory(S_2d);

   FunctionCoefficient ste_coeff_2d(ste0_static_2d);
   ste_coeff_2d.SetTime(t);
   ste_gf_2d.ProjectCoefficient(ste_coeff_2d);
   ste_gf_2d.SyncAliasMemory(S_2d);

   // PLF to build mass vector
   FunctionCoefficient rho_coeff_2d(rho0_static_2d); 
   rho_coeff_2d.SetTime(t);
   ParLinearForm *m_2d = new ParLinearForm(&L2FESpace_2d);
   m_2d->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff_2d));
   m_2d->Assemble();

   cout << "gridFunctions initiated.\n";

   mfem::hydroLO::LagrangianLOOperator<dim_2d> hydro_2d(H1FESpace_2d, H1FESpace_2d_L, L2FESpace_2d, L2VFESpace_2d, CRFESpace_2d, m_2d, problem_class_2d, use_viscosity, _mm, CFL);

   cout << "Done constructing hydro op\n";

   // Verify hardcoded timestep is small enough
   hydro_2d.BuildDijMatrix(S_2d);
   hydro_2d.CalculateTimestep(S_2d);
   if (dt_2d > hydro_2d.GetTimestep()) 
   {
      MFEM_ABORT("Sod hydro :: 1D Too big timestep\n");
   }

   is_collapsed = false;
   hydro_2d.MakeTimeStep(S_2d, t, dt_2d, is_collapsed);

   // Ensure the sub-vectors x_gf, v_gf, and ste_gf know the location of the
   // data in S. This operation simply updates the Memory validity flags of
   // the sub-vectors to match those of S.
   x_gf_2d.SyncAliasMemory(S_2d);
   mv_gf_2d.SyncAliasMemory(S_2d);
   sv_gf_2d.SyncAliasMemory(S_2d);
   v_gf_2d.SyncAliasMemory(S_2d);
   ste_gf_2d.SyncAliasMemory(S_2d);

   // Make sure that the mesh corresponds to the new solution state. This is
   // needed, because some time integrators use different S-type vectors
   // and the oper object might have redirected the mesh positions to those.
   pmesh_2d->NewNodes(x_gf_2d, false);

   ParGridFunction mc_gf_2d(&L2FESpace_2d);
   hydro_2d.CheckMassConservation(S_2d, mc_gf_2d);

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-Test2D-smooth_2d-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh_2d->Print(omesh);
   }

   cout << "mv_1d:\n";
   mv_gf_1d.Print(cout);

   cout << "mv_2d:\n";
   mv_gf_2d.Print(cout);

   // Verify that columns of 2d sim match columns of 1d sim
   int a0[4] = {0, 3, 12, 15};
   int a1[4] = {4, 7, 8, 11};
   int a2[4] = {1, 2, 13, 14};
   int a3[4] = {5, 6, 9, 10};

   double _sv0 = sv_gf_1d[0], _sv1 = sv_gf_1d[1], _sv2 = sv_gf_1d[2], _sv3 = sv_gf_1d[3];
   double _v0 = v_gf_1d[0], _v1 = v_gf_1d[1], _v2 = v_gf_1d[2], _v3 = v_gf_1d[3];
   double _ste0 = ste_gf_1d[0], _ste1 = ste_gf_1d[1], _ste2 = ste_gf_1d[2], _ste3 = ste_gf_1d[3];

   for (int i = 0; i < 4; i++)
   {
      // check sv
      if (abs(_sv0 - sv_gf_2d[a0[i]]) > tol || abs(_sv1 - sv_gf_2d[a1[i]]) > tol ||
          abs(_sv2 - sv_gf_2d[a2[i]]) > tol || abs(_sv3 - sv_gf_2d[a3[i]]) > tol)
      {
         cout << "Test2D::test_smooth_hydro failed on sv.\n";
         return 1;
      }
      
      // check v
      if (abs(_v0 - v_gf_2d[a0[i]]) > tol || abs(_v1 - v_gf_2d[a1[i]]) > tol ||
          abs(_v2 - v_gf_2d[a2[i]]) > tol || abs(_v3 - v_gf_2d[a3[i]]) > tol) 
      {
         cout << "Test2D::test_smooth_hydro failed on v.\n"; 

         cout << "i: " << i << endl;
         cout << _v0 << ", " << v_gf_2d[a0[i]] << endl;
         cout << _v1 << ", " << v_gf_2d[a1[i]] << endl;
         cout << _v2 << ", " << v_gf_2d[a2[i]] << endl;
         cout << _v3 << ", " << v_gf_2d[a3[i]] << endl;

         return 1; 
      }

      // check ste
      if (_ste0 != ste_gf_2d[a0[i]] || _ste1 != ste_gf_2d[a1[i]] ||
          _ste2 != ste_gf_2d[a2[i]] || _ste3 != ste_gf_2d[a3[i]]) 
      {
         cout << "Test2D::test_smooth_hydro failed on ste.\n"; 
         return 1; 
      }
   }

   return 0;
}

/*
Purpose:
   Plot various velocity fields from the smooth problem.
*/
void plot_mv_smooth()
{
   cout << "Test2D::plot_mv_smooth\n";
   a = 1., b = 0., c = 0., d = 0., e = 0., f = 0.;

   double t = 0., dt = .001;

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
      mesh_name << "../results/mesh-Test2D-smooth." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

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

   ParGridFunction v_exact_gf(&H1FESpace);

   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   mv_gf.MakeRef(&H1FESpace, S, offset[1]);
   sv_gf.MakeRef(&L2FESpace, S, offset[2]);
   v_gf.MakeRef(&L2VFESpace, S, offset[3]);
   ste_gf.MakeRef(&L2FESpace, S, offset[4]);

   Vector zero(dim);
   zero = 0.;
   VectorConstantCoefficient zero_coeff(zero);
   mv_gf.ProjectCoefficient(zero_coeff);
   mv_gf.SyncAliasMemory(S);

   pmesh->SetNodalGridFunction(&x_gf);
   x_gf.SyncAliasMemory(S);

   // Initialize specific volume, velocity, and specific total energy
   // Change class variables into static std::functions since virtual static member functions are not an option
   // and Coefficient class requires std::function arguments
   ProblemBase<dim> * problem_class = new SmoothWave<dim>();
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> sv0_static = 
      std::bind(&ProblemBase<dim>::sv0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<void(const Vector &, const double, Vector &)> v0_static = 
      std::bind(&ProblemBase<dim>::v0, problem_class, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
   std::function<double(const Vector &,const double)> ste0_static = 
      std::bind(&ProblemBase<dim>::ste0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> rho0_static = 
      std::bind(&ProblemBase<dim>::rho0, problem_class, std::placeholders::_1, std::placeholders::_2);

   VectorFunctionCoefficient v_exact_coeff(dim, v0_static);

   // Initialize specific volume, velocity, and specific total energy
   FunctionCoefficient sv_coeff(sv0_static);
   sv_coeff.SetTime(t);
   sv_gf.ProjectCoefficient(sv_coeff);
   sv_gf.SyncAliasMemory(S);

   v_exact_coeff.SetTime(t);
   v_gf.ProjectCoefficient(v_exact_coeff);
   // Sync the data location of v_gf with its base, S
   v_gf.SyncAliasMemory(S);

   v_exact_gf.ProjectCoefficient(v_exact_coeff);

   FunctionCoefficient ste_coeff(ste0_static);
   ste_coeff.SetTime(t);
   ste_gf.ProjectCoefficient(ste_coeff);
   ste_gf.SyncAliasMemory(S);

   // PLF to build mass vector
   FunctionCoefficient rho_coeff(rho0_static); 
   rho_coeff.SetTime(t);
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff));
   m->Assemble();

   cout << "gridFunctions initiated.\n";

   mfem::hydroLO::LagrangianLOOperator<dim> hydro(H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   cout << "Done constructing hydro op\n";

   Vector _vel(dim), vec_res(dim);
   DenseMatrix dm(dim);

   hydro.CalculateTimestep(S);
   dt = hydro.GetTimestep();
   cout << "dt: " << dt << endl;

   hydro.ComputeMeshVelocities(S, S, t, dt); // For now just use S as S_old. Will need to replace
   cout << "Done computing mesh velocities\n";

   /* ************************
   Displace Velocities
   *************************** */ 
   socketstream vis_vh, vis_vexact, vis_vcrgf, vis_ifv_exact;
   char vishost[] = "localhost";
   int visport = 19916;

   MPI_Barrier(pmesh->GetComm());

   int Wx = 0, Wy = 0;
   const int Ww = 350, Wh = 350;
   int offx = Ww+10, offy = Wh+45;

   hydroLO::VisualizeField(vis_vh, vishost, visport, mv_gf, "2D Smooth: Reconstructed Mesh Velocity", Wx, Wy, Ww, Wh);

   Wx += offx;

   hydroLO::VisualizeField(vis_vexact, vishost, visport, v_exact_gf, "2D Smooth: Exact Velocity", Wx, Wy, Ww, Wh);

   Wx += offx;
   ParGridFunction v_CR_gf(&CRFESpace);
   hydro.GetVCRgf(v_CR_gf);

   hydroLO::VisualizeField(vis_vcrgf, vishost, visport, v_CR_gf, "2D Smooth: v_CR_gf", Wx, Wy, Ww, Wh);

   Wx += offx;

   ParGridFunction ifv_exact(&CRFESpace);
   ifv_exact.ProjectCoefficient(v_exact_coeff);
   hydroLO::VisualizeField(vis_ifv_exact, vishost, visport, ifv_exact, "2D Smooth: ifv_exact", Wx, Wy, Ww, Wh);

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
      mesh_name << "../results/mesh-Test2D-smooth-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   return;
}