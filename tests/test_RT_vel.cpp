#include "mfem.hpp"
#include "laglos_solver.hpp"
#include "problem_template.h"
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
       b = 0.,
       c = 0.,
       d = 0.,
       e = 0,
       f = 0.,
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
int mesh_refinements = 2;
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
int test_Vi_geo();
int test_Ci_geo();
int test_RT_vel();
int test_RT_int_grad();
int test_RT_int_grad2();
int test_RT_int_grad_quadratic();
int test_determinant();
void plot_velocities();
void test_vel_field_1();
void test_vel_field_2();

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
   d += test_Vi_geo();
   d += test_Ci_geo();
   d += test_RT_vel();
   d += test_RT_int_grad();
   d += test_RT_int_grad2();
   d += test_RT_int_grad_quadratic();
   d += test_determinant();
   // plot_velocities();
   // test_vel_field_1();
   // test_vel_field_2();

   return d;
}

void velocity_exact(const Vector &x, const double & t, Vector &v)
{
   if (dim != 2)
   {
      MFEM_ABORT("Dimension must be equal to 2 for mesh motion test.\n");
   }
   v[0] = aq * pow(x[0],2) + bq * pow(x[1],2) + a * x[0] + b * x[1] + c;
   v[1] = dq * pow(x[0],2) + eq * pow(x[1],2) + d * x[0] + e * x[1] + f;
}

// void test_VDLFGI()
// {
//    // Initialize mesh
//    mfem::Mesh *mesh;
//    mesh = new mfem::Mesh(mesh_file, true, true);

//    // Refine the mesh
//    for (int lev = 0; lev < mesh_refinements; lev++)
//    {
//       mesh->UniformRefinement();
//    }

//    // Construct parmesh object
//    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);            
//    delete mesh;

//    // Output mesh to be visualized
//    // Can be visualized with glvis -np # -m mesh-test-moved
//    {
//       int precision = 12;
//       ostringstream mesh_name;
//       mesh_name << "../results/mesh-TestVDLFGI." << setfill('0') << setw(6) << myid;
//       ofstream omesh(mesh_name.str().c_str());
//       omesh.precision(precision);
//       pmesh->Print(omesh);
//    }

//    // Test stuff
//    Vector ones(2);
//    ones = 1.0;
//    VectorConstantCoefficient ones_coeff(ones);

//    CrouzeixRaviartFECollection CRFEC;
//    ParFiniteElementSpace CRFESpace(pmesh, &CRFEC, dim);

//    // Set int rule
//    IntegrationRules IntRulesLo(0, Quadrature1D::GaussLobatto);
//    const IntegrationRule gir = IntRules.Get(CRFESpace.GetFE(0)->GetGeomType(), 1);

//    VectorDomainLFGradIntegrator * vdlfgi = new VectorDomainLFGradIntegrator(ones_coeff);
//    vdlfgi->SetIntRule(&gir);

//    // Output integration rule information
//    cout << "Order of integration rule: " << gir.GetOrder() << endl;
//    cout << "number of integration points: " << gir.GetNPoints() << endl;
//    cout << "weights:\n";
//    const Array<double> weights = gir.GetWeights();
//    for (int i = 0; i < gir.GetNPoints(); i++)
//    {
//       cout << "weight at i = " << i << ": " << weights[i] << endl;
//       IntegrationPoint point = gir.IntPoint(i);
//       cout << "x: " << point.x << ", y: " << point.y << endl;
//    }

//    // Print element vector for each cell
//    Vector elvec;
//    Vector cent;
//    DenseMatrix dshape(4,2);
//    for (int el_index = 0; el_index < pmesh->GetNE(); el_index++)
//    {
//       const FiniteElement * fe = CRFESpace.GetFE(el_index);
//       ElementTransformation * trans = pmesh->GetElementTransformation(el_index);

//       const int dim = fe->GetDim();
//       const int dof = fe->GetDof();
//       const int sdim = trans->GetSpaceDim();
//       const int order = fe->GetOrder();
//       cout << "el_index: " << el_index << ", dim: " << dim << ", dof: " << dof << ", sdim: " << sdim << ", order: " << order << endl;
//       trans->Transform(Geometries.GetCenter(pmesh->GetElementBaseGeometry(el_index)),cent);
//       vdlfgi->AssembleRHSElementVect(*fe, *trans, elvec);
//       // elvec.Print(cout);
//       cout << "cent: \n";
//       cent.Print(cout);
//       for (int i = 0; i < gir.GetNPoints(); i++)
//       {
//          const IntegrationPoint &ip = gir.IntPoint(i);
//          trans->SetIntPoint(&ip);

//          fe->CalcPhysDShape(*trans, dshape);
//          cout << "dshape for el " << el_index << " at integration point " << i << endl;
//          dshape.Print(cout);
//       }
//    }
   
//    ParLinearForm lf(&CRFESpace);
//    lf.AddDomainIntegrator(vdlfgi);
//    lf.Assemble();
//    lf.Print(cout);

//    /* Print out diagnostics */
//    cout << "num cells: " << pmesh->GetNE() << endl;
//    cout << "num interior faces: " << pmesh->GetNFbyType(FaceType::Interior) << endl;
//    cout << "num boundary faces: " << pmesh->GetNFbyType(FaceType::Boundary) << endl;
//    cout << "num total faces: " << pmesh->GetNumFaces() << endl;
//    cout << "CR num DoFs: " << CRFESpace.GetNDofs() << endl;
//    cout << "CRFESpace.GetNE(): " << CRFESpace.GetNE() << endl;
//    cout << "CR true dofs: " << CRFESpace.GetTrueVSize() << endl;

// }

/*
Purpose:
   In a linear velocity field, Vgeo should be the same as the exact velocity at the given node.  This 
   iterates over all nodes and validates this.
*/
int test_Vi_geo()
{
   cout << "Testing Vi_geo computation\n";
   aq = 0., bq =0., dq = 0., eq = 0.;
   a = 2., b = 4., c = 6., d = 8., e = 9., f = 10.;

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

   ProblemBase<dim> * problem_class = new ProblemTemplate<dim>();

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   // Necessary parameters
   Vector Vgeo(dim), V_ex(dim), V_error(dim);
   double _error = 0.;
   double dt = 1., t = 0.;

   // Must compute the intermediate face velocities before computing geometric velocity
   hydro.compute_intermediate_face_velocities(S, t, "testing", &velocity_exact);

   // Iterate across all nodes and verify exact velocity and Vigeo is the same at all the 
   // geometric vertices.
   for (int node_it = 0; node_it < H1FESpace.GetNDofs() - L2FESpace.GetNDofs(); node_it++)
   {
      hydro.compute_geo_V(node_it, Vgeo);
      hydro.get_node_velocity(S, node_it, V_ex);

      cout << "--- Velocities on node: " << node_it << " ---\n";
      cout << "geo V: ";
      Vgeo.Print(cout);
      cout << "V_ex: ";
      V_ex.Print(cout);

      subtract(V_ex, Vgeo, V_error);
      _error += V_error.Norml2();
   }

   if (abs(_error) < tol)
   {
      return 0; // Test passed
   }
   else 
   {
      cout << "failed test_Vi_geo.\n";
      cout << "error: " << _error << endl;
      return 1; // Test failed
   }
}

/*
Purpose:
   In a linear velocity field, Cgeo should be the same at every node.  This 
   iterates over all nodes and validates this.
*/
int test_Ci_geo()
{
   cout << "test_Ci_geo\n";
   aq = 0., bq = 0., dq = 0., eq = 0.;
   a = 2., b = 4., c = 6., d = 8., e = 9., f = 10.;
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

   ProblemBase<dim> * problem_class = new ProblemTemplate<dim>();

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   // Necessary parameters
   DenseMatrix Cgeo(dim), Cgeo_0(dim), C_error(dim);
   double _error = 0.;
   double dt = 1., t = 0.;

   // Must compute the intermediate face velocities before computing geometric velocity
   hydro.compute_intermediate_face_velocities(S, t, "testing", &velocity_exact);

   // Iterate across all nodes and verify exact velocity and Vigeo is the same at all the 
   // geometric vertices.
   hydro.compute_geo_C(0, Cgeo_0);
   for (int node_it = 1; node_it < H1FESpace.GetNDofs() - L2FESpace.GetNDofs(); node_it++)
   {
      hydro.compute_geo_C(node_it, Cgeo);
      Add(Cgeo, Cgeo_0, -1., C_error);
      _error += C_error.FNorm();

      cout << "--- Ci on node: " << node_it << " ---\n";
      Cgeo.Print(cout);
   }

   if (abs(_error) < 1e-10)
   {
      return 0; // Test passed
   }
   else 
   {
      cout << "failed test_Ci_geo\n";
      cout << "error: " << _error << endl;
      return 1; // Test failed
   }
}

/*
Purpose:
   The purpose of this function is to verify our RT velocity field satisfies equation (5.10), 
   i.e. that

    1  /
   --- | v_h^v ds = V_F.
   |F| /F   
*/
int test_RT_vel()
{
   cout << "test_RT_vel\n";
   aq = 0., bq = 0., dq = 0., eq = 0.;
   a = 1., b = 2., c = 3., d = -4., e = -5., f = -6.;
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
      mesh_name << "../results/mesh-TestRTvel." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   // Template FE stuff to construct hydro operator
   H1_FECollection H1FEC(order_mv, dim);
   L2_FECollection L2FEC(order_u, dim, BasisType::Positive);
   CrouzeixRaviartFECollection CRFEC; // defaults to order 1

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

   ProblemBase<dim> * problem_class = new SodProblem<dim>();

   Vector one(dim);
   one = 1.;

   VectorFunctionCoefficient v_exact_coeff(dim, &velocity_exact);
   mv_gf.ProjectCoefficient(v_exact_coeff);
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

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   // Compute intermediate face velocities
   double t = 0., error = 0.;
   int num_exact = 0;

   hydro.compute_intermediate_face_velocities(S, t, "testing", &velocity_exact);
   // hydro.compute_intermediate_face_velocities(S, t);
   ParGridFunction v_CR_gf(&CRFESpace), v_CR_gf_int_el1(&CRFESpace);
   hydro.get_vcrgf(v_CR_gf);

   // Set integration rule for face
   cout << "setting integration rule\n";
   IntegrationRules IntRulesLo(0, Quadrature1D::GaussLobatto);
   IntegrationRule IR_face = IntRules.Get(CRFESpace.GetFaceElement(0)->GetGeomType(), 2);
   cout << "face geometry type: " << CRFESpace.GetFaceElement(0)->GetGeomType() << endl;
   cout << "num Q points: " << IR_face.GetNPoints() << endl;
   cout << "finished setting integration rule\n";
   mfem::Mesh::FaceInformation FI;
   Vector n_int(dim), n_vec(dim), Vf(dim), val(dim), int_vel(dim), el1_int_vel(dim), el2_int_vel(dim);
   Vector el1_error(dim), el2_error(dim);

   // Iterate over faces
   for (int face = 0; face < L2FESpace.GetNF(); face++) // face iterator
   {
      cout << "---face--- " << face << endl;
      el1_int_vel = 0.;
      el2_int_vel = 0.;
      FI = pmesh->GetFaceInformation(face);
      hydro.CalcOutwardNormalInt(S, FI.element[0].index, face, n_int);
      n_vec = n_int;
      double F = n_vec.Norml2();

      hydro.get_intermediate_face_velocity(face, Vf);
      FaceElementTransformations * FEtrans = pmesh->GetFaceElementTransformations(face);
      ElementTransformation & trans_el1 = FEtrans->GetElement1Transformation();

      if (FI.IsInterior()) { 
         cout << "interior face\n"; 
         
         ElementTransformation & trans_el2 = FEtrans->GetElement2Transformation();

         // Iterate over quadrature
         for (int j = 0; j < IR_face.GetNPoints(); j++)
         {
            const IntegrationPoint &ip_f = IR_face.IntPoint(j);
            FEtrans->SetAllIntPoints(&ip_f);

            const IntegrationPoint &ip_el1 = FEtrans->GetElement1IntPoint();
            const IntegrationPoint &ip_el2 = FEtrans->GetElement2IntPoint();

            v_CR_gf.GetVectorValue(trans_el1, ip_el1, val);
            el1_int_vel.Add(ip_f.weight, val);

            v_CR_gf.GetVectorValue(trans_el2, ip_el2, val);
            el2_int_vel.Add(ip_f.weight, val);            
         }
         
         // Compute face contribution to error
         subtract(el1_int_vel, Vf, el1_error);
         subtract(el1_int_vel, Vf, el2_error);
         double e1 = el1_error.Norml2();
         double e2 = el2_error.Norml2();

         if (e1 < tol && e2 < tol)
         {
            num_exact += 1;
         }
         else
         {
            cout << "Unequal velocities:\n";
            cout << "el1_int_vel corresponding to cell " << FI.element[0].index << ": ";
            el1_int_vel.Print(cout);
            cout << "el2_int_vel corresponding to cell " << FI.element[1].index << ": ";
            el2_int_vel.Print(cout);
            cout << "Vf: ";
            Vf.Print(cout);
         }
         error += e1;
         error += e2;
      } 
      else { 
         cout << "boundary face\n"; 

         // Iterate over quadrature
         for (int j = 0; j < IR_face.GetNPoints(); j++)
         {
            const IntegrationPoint &ip_f = IR_face.IntPoint(j);
            FEtrans->SetAllIntPoints(&ip_f);

            const IntegrationPoint &ip_el1 = FEtrans->GetElement1IntPoint();

            v_CR_gf.GetVectorValue(trans_el1, ip_el1, val);
            el1_int_vel.Add(ip_f.weight, val);
         }
         
         // Compute face contribution to error
         subtract(el1_int_vel, Vf, el1_error);
         double e1 = el1_error.Norml2();
         
         if (e1 < tol)
         {
            num_exact += 1;
         }
         else
         {
            cout << "Unequal velocities:\n";
            cout << "cell_int_velocities: ";
            el1_int_vel.Print(cout);
            cout << "Vf: ";
            Vf.Print(cout);
         }

         error += e1;
      }
   }

   double percent_exact = (double)num_exact / (double)L2FESpace.GetNF();
   cout << "norm of error: " << error << endl;
   cout << "percent exact: " << percent_exact << endl;
   if (percent_exact < 1.)
   {
      return 1;
   }
   else
   {
      return 0;
   }
}

/*
Purpose:
   Test the integral of the gradient of the RT velocity function given an integration rule.
   This test case uses a linear prescribed velocity.   
*/
int test_RT_int_grad()
{
   cout << "Testing RT_int_grad function\n";
   // Ensure there is no quadratic component to the velocity field
   aq = 0., bq = 0., dq = 0., eq = 0.;
   a = 1., b = 2., c = 3., d = 4., e = 5., f = 6.;

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

   ProblemBase<dim> * problem_class = new ProblemTemplate<dim>();

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   // Necessary parameters
   DenseMatrix dm(dim), true_grad(dim), grad_error(dim);
   true_grad(0,0) = a;
   true_grad(0,1) = b;
   true_grad(1,0) = d;
   true_grad(1,1) = e;;
   double _error = 0.;
   double dt = 1., t = 0.;

   // Output exact gradient
   cout << "The true gradient matrix should be: \n";
   true_grad.Print(cout);

   // Must compute the intermediate face velocities before computing geometric velocity
   hydro.compute_intermediate_face_velocities(S, t, "testing", &velocity_exact);

   for (int cell_it = 0; cell_it < L2FESpace.GetNE(); cell_it++)
   {
      hydro.RT_int_grad(cell_it, dm);
      cout << "Dense Matrix for cell: " << cell_it << endl;
      dm.Print(cout);
      // assert(false);

      Add(dm, true_grad, -1., grad_error);
      _error += grad_error.FNorm();
   }
   
   if (abs(_error) < 1e-10)
   {
      return 0; // Test passed
   }
   else 
   {
      cout << "failed test_RT_int_grad\n";
      cout << "error: " << _error << endl;
      return 1; // Test failed
   }
}

/*
Purpose:
   Test the integral of the gradient of the RT velocity function given an integration rule.
   This test case uses a linear prescribed velocity.   
*/
int test_RT_int_grad2()
{
   cout << "Testing RT_int_grad function\n";
   // Ensure there is no quadratic component to the velocity field
   aq = 0., bq = 0., dq = 0., eq = 0.;
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

   ProblemBase<dim> * problem_class = new ProblemTemplate<dim>();

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   // Necessary parameters
   DenseMatrix dm(dim), true_grad(dim), grad_error(dim);
   true_grad(0,0) = a;
   true_grad(0,1) = b;
   true_grad(1,0) = d;
   true_grad(1,1) = e;;
   double _error = 0.;
   double dt = 1., t = 0.;

   // Output exact gradient
   cout << "The true gradient matrix should be: \n";
   true_grad.Print(cout);

   // Must compute the intermediate face velocities before computing geometric velocity
   hydro.compute_intermediate_face_velocities(S, t, "testing", &velocity_exact);

   for (int cell_it = 0; cell_it < L2FESpace.GetNE(); cell_it++)
   {
      hydro.RT_int_grad(cell_it, dm);
      cout << "Dense Matrix for cell: " << cell_it << endl;
      dm.Print(cout);
      // assert(false);

      Add(dm, true_grad, -1., grad_error);
      _error += grad_error.FNorm();
   }
   
   if (abs(_error) < 1e-10)
   {
      return 0; // Test passed
   }
   else 
   {
      cout << "failed test_RT_int_grad2\n";
      cout << "error: " << _error << endl;
      return 1; // Test failed
   }
}

/*
Purpose:
   Test the integral of the gradient of the RT velocity function given an integration rule.
   This test case uses a quadratic prescribed velocity. 

   True grad = | x+2   3y+4  |
               | 5x+8  7y+9 |
*/
int test_RT_int_grad_quadratic()
{
   cout << "Testing RT_int_grad_quadratic function\n";
   // Ensure there is no quadratic component to the velocity field
   aq = 1./2., bq = 3./2., dq = 5./2., eq = 7./2.;
   a = 2., b = 4., c = 6., d = 8., e = 9., f = 10.;

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

   ProblemBase<dim> * problem_class = new ProblemTemplate<dim>();

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   // Necessary parameters
   DenseMatrix dm(dim), true_grad(dim), grad_error(dim);
   double _error = 0.;
   double dt = 1., t = 0.;

   // Must compute the intermediate face velocities before computing geometric velocity
   hydro.compute_intermediate_face_velocities(S, t, "testing", &velocity_exact);
   Vector cell_x(dim);
   for (int cell_it = 0; cell_it < L2FESpace.GetNE(); cell_it++)
   {
      // Compute true gradient as it will not be constant across mesh
      int cell_vdof = cell_it + H1FESpace.GetNVDofs() + L2FESpace.GetNF();
      hydro.get_node_position(S, cell_vdof, cell_x);
      true_grad(0,0) = cell_x[0] + 2.;
      true_grad(0,1) = 3.*cell_x[1] + 4.;
      true_grad(1,0) = 5.*cell_x[0] + 8.;
      true_grad(1,1) = 7.*cell_x[1] + 9.;

      // Compute approximate gradient
      hydro.RT_int_grad(cell_it, dm);
      cout << "Dense Matrix for cell: " << cell_it << endl;
      dm.Print(cout);
      // Output exact gradient
      cout << "The true gradient matrix should be: \n";
      true_grad.Print(cout);
      // assert(false);

      Add(dm, true_grad, -1., grad_error);
      _error += grad_error.FNorm();
   }
   
   if (abs(_error) < 1e-10)
   {
      return 0; // Test passed
   }
   else 
   {
      cout << "failed test_RT_int_grad_quadratic\n";
      cout << "error: " << _error << endl;
      return 1; // Test failed
   }
}


int test_determinant()
{
   cout << "test_determinant\n";

   // Since compute_determinant is not a static member function, I will instantiate a baseline hydro operator
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

   ProblemBase<dim> * problem_class = new ProblemTemplate<dim>();

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   // Bread and butter of test
   // First determinant test
   DenseMatrix C(2);
   C = 1.;
   double dt = 1.;
   double d = 0.;

   hydro.compute_determinant(C, dt, d);
   if (d != 2.)
   {
      cout << "Failed 1st determinant test.\n";
      return 1;
   }
   else { cout << "Passed 1st derminant test.\n"; }

   // Second determinant test
   C = 0.;
   C(0,0) = 1.;
   d = 0.;
   hydro.compute_determinant(C, dt, d);
   if (d != 1.5)
   {
      cout << "Failed 2nd determinant test.\n";
      return 1;
   }
   else { cout << "Passed 2nd derminant test.\n"; }

   // Third determinant test
   C = -1.;
   C(0,1) = 5.;
   C(1,0) = 1.;
   d = 0.;
   hydro.compute_determinant(C, dt, d);
   if (d != 1.)
   {
      cout << "Failed 3rd determinant test.\n";
      return 1;
   }
   else { cout << "Passed 3rd derminant test.\n"; }

   // Fourth determinant test
   C = 0.;
   C(0,0) = 4.;
   C(0,1) = -8.;
   C(1,0) = 1.;
   C(1,1) = 4.;

   d = 0.;
   hydro.compute_determinant(C, dt, d);
   if (d != 3.)
   {
      cout << "Failed 4th determinant test.\n";
      return 1;
   }
   else { cout << "Passed 4th derminant test.\n"; }

   return 0;
}

void plot_velocities()
{
   a = 1., b = 0., c = 0., d = 0., e = -0.5, f = 0.;
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

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-TestVDLFGI." << setfill('0') << setw(6) << myid;
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

   ProblemBase<dim> * problem_class = new ProblemTemplate<dim>();

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   Vector _vel(dim), true_vel(dim);
   double t = 0.;

   hydro.compute_intermediate_face_velocities(S, t, "testing", &velocity_exact);

   VectorFunctionCoefficient v_CR_coeff(dim, &velocity_exact);
   v_cr_gf.ProjectCoefficient(v_CR_coeff);
   ParGridFunction v_exact_gf(&H1FESpace);
   v_exact_gf.ProjectCoefficient(v_CR_coeff);

   cout << "Printing v_cr_gf:\n";
   v_cr_gf.Print(cout);

   // Set int rule
   IntegrationRules IntRulesLo(0, Quadrature1D::GaussLobatto);
   const IntegrationRule gir = IntRules.Get(CRFESpace.GetFE(0)->GetGeomType(), 1);

   // Output integration rule information
   cout << "Order of integration rule: " << gir.GetOrder() << endl;
   cout << "number of integration points: " << gir.GetNPoints() << endl;
   cout << "weights:\n";
   const Array<double> weights = gir.GetWeights();
   for (int i = 0; i < gir.GetNPoints(); i++)
   {
      cout << "weight at i = " << i << ": " << weights[i] << endl;
      IntegrationPoint point = gir.IntPoint(i);
      cout << "x: " << point.x << ", y: " << point.y << endl;
   }
   ParFiniteElementSpace CRc(CRFESpace.GetParMesh(), CRFESpace.FEColl(), 1);
   const int size = CRc.GetVSize();
   cout << "Size: " << size << endl;
   ParGridFunction vgeo_gf(&H1FESpace);

   cout << "Diagnostics:\n";
   cout << "CRFESpace.GetVSize(): " << CRFESpace.GetVSize() << endl;
   cout << "CRFESpace.TrueVSize(): " << CRFESpace.TrueVSize() << endl;
   cout << "CRFESpace.GetNDofs(): " << CRFESpace.GetNDofs() << endl;

   cout << "CRc.GetVSize(): " << CRc.GetVSize() << endl;
   cout << "CRc.TrueVSize(): " << CRc.TrueVSize() << endl;
   cout << "CRc.GetNDofs(): " << CRc.GetNDofs() << endl;

   Array<double> cell_areas(L2FESpace.GetNE());
   Array<double> bmn_sum_cells(L2FESpace.GetNE());

   // Store initial mass of cells
   for (int ci = 0; ci < L2FESpace.GetNE(); ci++)
   {
      double k = pmesh->GetElementVolume(ci);

      cell_areas[ci] = k;
      assert(k >= 0.);
   }

   DenseMatrix res(dim);

   Vector x(dim), vec_res(dim);
   double dt = 1.;
   bool is_dt_changed = false;
   for (int node_it = 0; node_it < H1FESpace.GetNDofs() - L2FESpace.GetNDofs(); node_it++)
   {
      // hydro.get_node_position(S, node_it, x);
      // cout << "node position for node: " << node << endl;
      // x.Print(cout);

      hydro.compute_node_velocity_RT(node_it, dt, vec_res, is_dt_changed);
      hydro.update_node_velocity(S, node_it, vec_res);

      // Also store for plotting simple the averaged geometric V
      hydro.compute_geo_V(node_it, vec_res);
      for (int i = 0; i < dim; i++)
      {
         int index = node_it + i*H1FESpace.GetNDofs();
         vgeo_gf[index] = vec_res[i];
      }

      // restart nodal velocity computation if the timestep has been restricted
      if (is_dt_changed)
      {
         is_dt_changed = false;
         node_it = -1;
      }
   }

   // Fill center with average over corner vertices
   Array<int> verts;
   Vector vel_center(dim);
   for (int ci = 0; ci < L2FESpace.GetNDofs(); ci++)
   {
      double vel_center_x = 0., vel_center_y = 0.;
      pmesh->GetElementVertices(ci, verts);
      for (int j = 0; j < verts.Size(); j++)
      {
         int index = verts[j];
         vel_center_x += vgeo_gf[index];
         index = verts[j] + H1FESpace.GetNDofs();
         vel_center_y += vgeo_gf[index];
      }
      vel_center_x *= 1. / verts.Size();
      vel_center_y *= 1. / verts.Size();

      vgeo_gf[H1FESpace.GetNDofs() - L2FESpace.GetNDofs() + ci] = vel_center_x;
      vgeo_gf[2 * H1FESpace.GetNDofs() - L2FESpace.GetNDofs() + ci] = vel_center_y;
   }

   // hydro.compute_corrective_face_velocities(S, t, dt);
   hydro.fill_face_velocities_with_average(S);
   hydro.fill_center_velocities_with_average(S);

   Array<int> fids, oris, row;
   mfem::Mesh::FaceInformation FI;
   double bmn_sum;
   int c, cp;
   Vector Uc(dim+2), Ucp(dim+2), c_vec(dim), Vf(dim), face_x(dim), vdof1_x(dim), vdof2_x(dim), n_int(dim), n_vec(dim);

   // for (int ci = 0; ci < L2FESpace.GetNE(); ci++)
   // {
   //    // Compute sum of bmn on cell
   //    pmesh->GetElementEdges(ci, fids, oris);

   //    // iterate over faces
   //    bmn_sum = 0.;
   //    for (int j = 0; j < fids.Size(); j++)
   //    {
   //       Vf = 0.;
   //       int face = fids[j];
   //       FI = pmesh->GetFaceInformation(face);

   //       // Compute intermediate face velocity on fly
   //       c = FI.element[0].index;
   //       cp = FI.element[1].index;
   //       hydro.GetCellStateVector(S, c, Uc);

   //       H1FESpace.GetFaceDofs(fids[j], row);
   //       int face_vdof1 = row[1], face_vdof2 = row[0], face_dof = row[2];
   //       hydro.get_node_position(S, face_vdof1, vdof1_x);
   //       hydro.get_node_position(S, face_vdof2, vdof2_x);
   //       hydro.get_node_position(S, face_dof, face_x);
         
   //       velocity_exact(face_x, t, Vf);

   //       hydro.CalcOutwardNormalInt(S, ci, face, n_int);
   //       n_vec = n_int;
   //       double F = n_vec.Norml2();
   //       n_vec /= F;

   //       double bmn = Vf * n_vec;
   //       bmn *= F;
   //       bmn_sum += bmn;
         
   //    } // End face iterator
   //    bmn_sum_cells[ci] = bmn_sum;
   // }

   /* ************************
   Move the mesh according to previously computed nodal velocities
   *************************** */ 
   cout << "x_gf pre move:\n";
   x_gf.Print(cout);
   add(x_gf, dt, mv_gf, x_gf);
   pmesh->NewNodes(x_gf, false);
   cout << "x_gf post move:\n";
   x_gf.Print(cout);
   t += dt;
   if (!suppress_test_output)
   {
      cout << "Done moving mesh.\n";
   }

   double num = 0., error_denom = 0., _error = 0.;
   /* ************************
   Compute cell area error
   *************************** */ 
   std::vector<double> cells_py(L2FESpace.GetNE());
   std::vector<double> cell_errors_py(L2FESpace.GetNE());
   // for (int ci = 0; ci < L2FESpace.GetNE(); ci++)
   // {
   //    if (!suppress_test_output)
   //    {
   //       cout << "Computing cell area error for cell: " << ci << endl;
   //    }
   //    // Compute cell area side
   //    double k = pmesh->GetElementVolume(ci);
   //    double _bmn_approx = (k - cell_areas[ci]) / dt;
   //    cout << "cell area pre move: " << cell_areas[ci] << endl;
   //    cout << "cell area post move: " << k << endl;
   //    cout << "dt: " << dt << endl;

   //    cout << "Cell wise comparison:\n";
   //    cout << "cell area difference: " << _bmn_approx << ", bmn_sum: " << bmn_sum_cells[ci] << endl;

   //    double val = abs(_bmn_approx - bmn_sum_cells[ci]);
   //    num += val;
   //    error_denom += abs(bmn_sum_cells[ci]);
   //    cells_py[ci] = ci;
   //    cell_errors_py[ci] = val;
   // }

   // if (error_denom < tol)
   // {
   //    _error = 0.;
   // }
   // else
   // {
   //    _error = abs(num / error_denom);
   // }

   /* ************************
   Visualization
   *************************** */ 
   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-TestVDLFGI-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   socketstream vis_vh, vis_vgeo, vis_vexact;
   char vishost[] = "localhost";
   int visport = 19916;

   MPI_Barrier(pmesh->GetComm());
   vis_vgeo.precision(8);

   int Wx = 0, Wy = 0;
   const int Ww = 350, Wh = 350;
   int offx = Ww+10, offy = Wh+45;

   hydrodynamics::VisualizeField(vis_vgeo, vishost, visport, vgeo_gf, "Geometric velocity", Wx, Wy, Ww, Wh);

   Wx += offx;

   hydrodynamics::VisualizeField(vis_vh, vishost, visport, mv_gf, "Mesh Velocity", Wx, Wy, Ww, Wh);

   Wx += offx;

   hydrodynamics::VisualizeField(vis_vexact, vishost, visport, v_exact_gf, "Exact Velocity", Wx, Wy, Ww, Wh);

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
   
}

/*
Purpose:
   This function tests a the RT computation of a continuous velocity field from a prescribed 
   velocity field that is linear in the first component and constant in the second.

   v_exact = | 5 0 | |x| + |1|
             | 0 0 | |y|   |1|
*/
void test_vel_field_1()
{
   a = 5., b = 0., c = 1., d = 0., e = 3., f = 1.;
   aq = 0., bq = 0., dq = 0., eq = 0.;

   double t = 0., dt = 0.0001;

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
      mesh_name << "../results/mesh-TestVel1." << setfill('0') << setw(6) << myid;
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
   ParGridFunction vgeo_gf(&H1FESpace);
   ParGridFunction v_cr_gf(&CRFESpace);

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

   ProblemBase<dim> * problem_class = new ProblemTemplate<dim>();

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   Vector _vel(dim), vec_res(dim);
   DenseMatrix dm(dim);

   hydro.compute_intermediate_face_velocities(S, t, "testing", &velocity_exact);
   bool is_dt_changed = false;
   for (int node_it = 0; node_it < H1FESpace.GetNDofs() - L2FESpace.GetNDofs(); node_it++)
   {
      cout << "\nnode: " << node_it << endl;
      hydro.compute_geo_C(node_it, dm);
      cout << "Ci:\n";
      dm.Print(cout);
      hydro.compute_node_velocity_RT(node_it, dt, vec_res, is_dt_changed);
      // hydro.get_intermediate_face_velocity(node_it, vec_res);
      cout << "RT v: ";
      vec_res.Print(cout);
      hydro.update_node_velocity(S, node_it, vec_res);

      // Also store for plotting simple the averaged geometric V
      hydro.compute_geo_V(node_it, vec_res);
      cout << "Vgeo: ";
      vec_res.Print(cout);
      for (int i = 0; i < dim; i++)
      {
         int index = node_it + i*H1FESpace.GetNDofs();
         vgeo_gf[index] = vec_res[i];
      }

      // restart nodal velocity computation if the timestep has been restricted
      if (is_dt_changed)
      {
         is_dt_changed = false;
         node_it = -1;
      }
   }

   // Fill center with average over corner vertices
   Array<int> verts;
   Vector vel_center(dim);
   for (int ci = 0; ci < L2FESpace.GetNDofs(); ci++)
   {
      double vel_center_x = 0., vel_center_y = 0.;
      pmesh->GetElementVertices(ci, verts);
      for (int j = 0; j < verts.Size(); j++)
      {
         int index = verts[j];
         vel_center_x += vgeo_gf[index];
         index = verts[j] + H1FESpace.GetNDofs();
         vel_center_y += vgeo_gf[index];
      }
      vel_center_x *= 1. / verts.Size();
      vel_center_y *= 1. / verts.Size();

      vgeo_gf[H1FESpace.GetNDofs() - L2FESpace.GetNDofs() + ci] = vel_center_x;
      vgeo_gf[2 * H1FESpace.GetNDofs() - L2FESpace.GetNDofs() + ci] = vel_center_y;
   }

   // hydro.compute_corrective_face_velocities(S, t, dt);
   hydro.fill_face_velocities_with_average(S);
   hydro.fill_center_velocities_with_average(S);

   /* ************************
   Displace Velocities
   *************************** */ 
   socketstream vis_vh, vis_vgeo, vis_vexact;
   char vishost[] = "localhost";
   int visport = 19916;

   MPI_Barrier(pmesh->GetComm());
   vis_vgeo.precision(8);

   int Wx = 0, Wy = 0;
   const int Ww = 350, Wh = 350;
   int offx = Ww+10, offy = Wh+45;

   hydrodynamics::VisualizeField(vis_vgeo, vishost, visport, vgeo_gf, "Geometric velocity", Wx, Wy, Ww, Wh);

   Wx += offx;

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
      mesh_name << "../results/mesh-TestVel1-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }
}

/*
Purpose:
   This function tests a the RT computation of a continuous velocity field from a prescribed 
   velocity field that is linear in the first component and constant in the second.

   v_exact = | 1/2  3/2 | |x^2| + | 2 4 | |x| + |6 |
             | 5/2  7/2 | |y^2|   | 8 9 | |y|   |10|
*/
void test_vel_field_2()
{
   cout << "Testing quadratic velocity field\n";
   aq = 1./2., bq = 3./2., dq = 5./2., eq = 7./2.;
   a = 2., b = 4., c = 6., d = 8., e = 9., f = 10.;

   double t = 0., dt = 0.0001;

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
      mesh_name << "../results/mesh-TestVel1." << setfill('0') << setw(6) << myid;
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
   ParGridFunction vgeo_gf(&H1FESpace);
   ParGridFunction v_cr_gf(&CRFESpace);

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

   ProblemBase<dim> * problem_class = new ProblemTemplate<dim>();

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   Vector _vel(dim), vec_res(dim);
   DenseMatrix dm(dim);

   hydro.compute_intermediate_face_velocities(S, t, "testing", &velocity_exact);
   bool is_dt_changed = false;
   for (int node_it = 0; node_it < H1FESpace.GetNDofs() - L2FESpace.GetNDofs(); node_it++)
   {
      cout << "\nnode: " << node_it << endl;
      hydro.compute_geo_C(node_it, dm);
      cout << "Ci:\n";
      dm.Print(cout);
      hydro.compute_node_velocity_RT(node_it, dt, vec_res, is_dt_changed);
      // hydro.get_intermediate_face_velocity(node_it, vec_res);
      cout << "RT v: ";
      vec_res.Print(cout);
      hydro.update_node_velocity(S, node_it, vec_res);

      // Also store for plotting simple the averaged geometric V
      hydro.compute_geo_V(node_it, vec_res);
      cout << "Vgeo: ";
      vec_res.Print(cout);
      for (int i = 0; i < dim; i++)
      {
         int index = node_it + i*H1FESpace.GetNDofs();
         vgeo_gf[index] = vec_res[i];
      }

      // restart nodal velocity computation if the timestep has been restricted
      if (is_dt_changed)
      {
         is_dt_changed = false;
         node_it = -1;
      }
   }

   // Fill center with average over corner vertices
   Array<int> verts;
   Vector vel_center(dim);
   for (int ci = 0; ci < L2FESpace.GetNDofs(); ci++)
   {
      double vel_center_x = 0., vel_center_y = 0.;
      pmesh->GetElementVertices(ci, verts);
      for (int j = 0; j < verts.Size(); j++)
      {
         int index = verts[j];
         vel_center_x += vgeo_gf[index];
         index = verts[j] + H1FESpace.GetNDofs();
         vel_center_y += vgeo_gf[index];
      }
      vel_center_x *= 1. / verts.Size();
      vel_center_y *= 1. / verts.Size();

      vgeo_gf[H1FESpace.GetNDofs() - L2FESpace.GetNDofs() + ci] = vel_center_x;
      vgeo_gf[2 * H1FESpace.GetNDofs() - L2FESpace.GetNDofs() + ci] = vel_center_y;
   }

   // hydro.compute_corrective_face_velocities(S, t, dt);
   hydro.fill_face_velocities_with_average(S);
   hydro.fill_center_velocities_with_average(S);

   /* ************************
   Displace Velocities
   *************************** */ 
   socketstream vis_vh, vis_vgeo, vis_vexact;
   char vishost[] = "localhost";
   int visport = 19916;

   MPI_Barrier(pmesh->GetComm());
   vis_vgeo.precision(8);

   int Wx = 0, Wy = 0;
   const int Ww = 350, Wh = 350;
   int offx = Ww+10, offy = Wh+45;

   hydrodynamics::VisualizeField(vis_vgeo, vishost, visport, vgeo_gf, "Geometric velocity", Wx, Wy, Ww, Wh);

   Wx += offx;

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
      mesh_name << "../results/mesh-TestVel1-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }
}