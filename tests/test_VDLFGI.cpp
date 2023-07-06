#include "mfem.hpp"
#include "laglos_solver.hpp"
#include <cassert>
#include <fstream>
#include <sstream>
#include <cmath>
#include "matplotlibcpp.h"

using namespace std;
using namespace mfem;

/* ---------------- Parameters to be used for tests ---------------- */
// Linear velocity field
const double a = 2.,
             b = 1.,
             c = 1.,
             d = 1.,
             e = -1.,
             f = 1.,
             aq = 1.,
             bq = 2.,
             dq = 3.,
             eq = 4.;

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

void velocity_exact(const Vector &x, const double & t, Vector &v);
void test_VDLFGI();
void test_integral();

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

   // test_VDLFGI();
   test_integral();
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

void test_VDLFGI()
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

   // Test stuff
   Vector ones(2);
   ones = 1.0;
   VectorConstantCoefficient ones_coeff(ones);

   CrouzeixRaviartFECollection CRFEC;
   ParFiniteElementSpace CRFESpace(pmesh, &CRFEC, dim);

   // Set int rule
   IntegrationRules IntRulesLo(0, Quadrature1D::GaussLobatto);
   const IntegrationRule gir = IntRules.Get(CRFESpace.GetFE(0)->GetGeomType(), 1);

   VectorDomainLFGradIntegrator * vdlfgi = new VectorDomainLFGradIntegrator(ones_coeff);
   vdlfgi->SetIntRule(&gir);

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

   // Print element vector for each cell
   Vector elvec;
   Vector cent;
   DenseMatrix dshape(4,2);
   for (int el_index = 0; el_index < pmesh->GetNE(); el_index++)
   {
      const FiniteElement * fe = CRFESpace.GetFE(el_index);
      ElementTransformation * trans = pmesh->GetElementTransformation(el_index);

      const int dim = fe->GetDim();
      const int dof = fe->GetDof();
      const int sdim = trans->GetSpaceDim();
      const int order = fe->GetOrder();
      cout << "el_index: " << el_index << ", dim: " << dim << ", dof: " << dof << ", sdim: " << sdim << ", order: " << order << endl;
      trans->Transform(Geometries.GetCenter(pmesh->GetElementBaseGeometry(el_index)),cent);
      vdlfgi->AssembleRHSElementVect(*fe, *trans, elvec);
      // elvec.Print(cout);
      cout << "cent: \n";
      cent.Print(cout);
      for (int i = 0; i < gir.GetNPoints(); i++)
      {
         const IntegrationPoint &ip = gir.IntPoint(i);
         trans->SetIntPoint(&ip);

         fe->CalcPhysDShape(*trans, dshape);
         cout << "dshape for el " << el_index << " at integration point " << i << endl;
         dshape.Print(cout);
      }
   }
   
   ParLinearForm lf(&CRFESpace);
   lf.AddDomainIntegrator(vdlfgi);
   lf.Assemble();
   lf.Print(cout);

   /* Print out diagnostics */
   cout << "num cells: " << pmesh->GetNE() << endl;
   cout << "num interior faces: " << pmesh->GetNFbyType(FaceType::Interior) << endl;
   cout << "num boundary faces: " << pmesh->GetNFbyType(FaceType::Boundary) << endl;
   cout << "num total faces: " << pmesh->GetNumFaces() << endl;
   cout << "CR num DoFs: " << CRFESpace.GetNDofs() << endl;
   cout << "CRFESpace.GetNE(): " << CRFESpace.GetNE() << endl;
   cout << "CR true dofs: " << CRFESpace.GetTrueVSize() << endl;

}



void test_integral()
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

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-TestIntegral." << setfill('0') << setw(6) << myid;
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
   
   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-TestIntegral." << setfill('0') << setw(6) << myid;
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
   ParGridFunction v_cr_gf(&CRFESpace);
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

   int node = 15;
   int cell = 1;
   Vector _vel(dim), true_vel(dim);
   double t = 0.;
   hydro.compute_intermediate_face_velocities(S, t, "testing", &velocity_exact);

   cout << "Printing v_cr_gf:\n";
   v_cr_gf.Print(cout);

   VectorFunctionCoefficient v_CR_coeff(dim, &velocity_exact);
   v_cr_gf.ProjectCoefficient(v_CR_coeff);

   cout << "Printing v_cr_gf:\n";
   v_cr_gf.Print(cout);




   // hydro.RT_corner_velocity(cell, node, _vel);
   // Set int rule
   IntegrationRules IntRulesLo(0, Quadrature1D::GaussLobatto);
   const IntegrationRule gir = IntRules.Get(CRFESpace.GetFE(0)->GetGeomType(), 2);

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
   ParGridFunction CRc_gf(&CRc);

   cout << "Diagnostics:\n";
   cout << "CRFESpace.GetVSize(): " << CRFESpace.GetVSize() << endl;
   cout << "CRFESpace.TrueVSize(): " << CRFESpace.TrueVSize() << endl;
   cout << "CRFESpace.GetNDofs(): " << CRFESpace.GetNDofs() << endl;

   cout << "CRc.GetVSize(): " << CRc.GetVSize() << endl;
   cout << "CRc.TrueVSize(): " << CRc.TrueVSize() << endl;
   cout << "CRc.GetNDofs(): " << CRc.GetNDofs() << endl;


   // Print element vector for each cell
   for (int c = 0; c < dim; c++)
   {
      cout << "--- Iterating on dimension: " << c << " ---\n";
      // Get dimensional component of pargridfunction
      CRc_gf.MakeRef(&CRc, v_cr_gf, c*size);
      cout << "CRc_gf:\n";
      CRc_gf.Print(cout);

      for (int el_index = 0; el_index < pmesh->GetNE(); el_index++)
      {
         const FiniteElement * fe = CRFESpace.GetFE(el_index);
         ElementTransformation * trans = pmesh->GetElementTransformation(el_index);

         const int dim = fe->GetDim();
         const int dof = fe->GetDof();
         const int sdim = trans->GetSpaceDim();
         const int order = fe->GetOrder();
         cout << "el_index: " << el_index << ", dim: " << dim << ", dof: " << dof << ", sdim: " << sdim << ", order: " << order << endl;
         // trans->Transform(Geometries.GetCenter(pmesh->GetElementBaseGeometry(el_index)),cent);

         mfem::Vector grad;
         for (int i = 0; i < gir.GetNPoints(); i++)
         {
            const IntegrationPoint &ip = gir.IntPoint(i);
            trans->SetIntPoint(&ip);

            cout << "el " << el_index << " at integration point " << i << endl;
            CRc_gf.GetGradient(*trans, grad);
            cout << "grad: \n";
            grad.Print(cout); 

         }
      }
   }
}