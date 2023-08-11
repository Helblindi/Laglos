#include "mfem.hpp"
#include "laglos_solver.hpp"
#include "sod.h"
#include "var-config.h"
#include <cassert>
#include <fstream>
#include <sstream>
#include <cmath>
#include "matplotlibcpp.h"

using namespace std;
using namespace mfem;
using namespace hydrodynamics;
namespace plt = matplotlibcpp;

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
bool _mm = false;          // Doesn't matter
double CFL = 0.5;          // Doesn't matter
static int num_procs = 0;
static int myid = -1;
bool suppress_test_output = false;
bool create_table = false;
int lower_refinement = 2;
int upper_refinement = 7;
/* ---------------- End Parameters ---------------- */

// void RT_int_grad(ParGridFunction & CR_v_gf, 
//                  ParMesh * pmesh,   
//                  const IntegrationRule * ir, 
//                  const int cell, 
//                  DenseMatrix & res);
void velocity_exact(const Vector &x, const double & t, Vector &v);
int test_flux();
int test_vel_field_1();
int test_CSV_getter_setter(); // LagrangianLOOperator::GetCellStateVector tester

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
   d += test_flux();
   d += test_vel_field_1();
   d += test_CSV_getter_setter();
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

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

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
   velocity field that is linear in the first component and constant in the second.

   v_exact = | 1 0 | |x| + |1|
             | 0 0 | |y|   |1|
*/
int test_vel_field_1()
{
   cout << "Test2D::test_vel_field_1\n";
   a = 1., b = 0., c = 0., d = 0., e = 1., f = 0.;

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

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);


   cout << "Done constructing hydro op\n";

   Vector _vel(dim), vec_res(dim);
   double t = 0., dt = 1;
   DenseMatrix dm(dim);

   // hydro.ComputeMeshVelocities(S, t, dt, flag, &velocity_exact);
   hydro.ComputeMeshVelocities(S, t, dt);
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

   hydrodynamics::VisualizeField(vis_vh, vishost, visport, mv_gf, "2D: Reconstructed Mesh Velocity", Wx, Wy, Ww, Wh);

   Wx += offx;

   hydrodynamics::VisualizeField(vis_vexact, vishost, visport, v_exact_gf, "2D: Exact Velocity", Wx, Wy, Ww, Wh);

   Wx += offx;
   ParGridFunction v_CR_gf(&CRFESpace);
   hydro.get_vcrgf(v_CR_gf);

   hydrodynamics::VisualizeField(vis_vcrgf, vishost, visport, v_CR_gf, "2D: v_CR_gf", Wx, Wy, Ww, Wh);

   Wx += offx;

   ParGridFunction ifv_exact(&CRFESpace);
   ifv_exact.ProjectCoefficient(v_exact_coeff);
   hydrodynamics::VisualizeField(vis_ifv_exact, vishost, visport, ifv_exact, "2D: ifv_exact", Wx, Wy, Ww, Wh);


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
   mesh_refinements = 3;

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
      mesh_name << "../results/mesh-Test2D-CSV." << setfill('0') << setw(6) << myid;
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

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

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