#include "mfem.hpp"
#include <cassert>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;
using namespace mfem;

// Problem
const int dim = 2;
int order_mv = 2;          // Order of mesh movement approximation space
int order_u = 0;
const string flag = "testing";
double tol = 1e-12;
const char *mesh_file = "../data/ref-square.mesh";
int mesh_refinements = 1;
static int num_procs = 0;
static int myid = -1;
/* ---------------- End Parameters ---------------- */

void test_VDLFGI();

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

   test_VDLFGI();
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

   // Output first integration point
   IntegrationPoint point = gir.IntPoint(0);

   cout << "x: " << point.x << ", y: " << point.y << endl;

   // Output integration rule information
   cout << "Order of integration rule: " << gir.GetOrder() << endl;
   cout << "number of integration points: " << gir.GetNPoints() << endl;
   cout << "weights:\n";
   const Array<double> weights = gir.GetWeights();
   for (int i = 0; i < gir.GetNPoints(); i++)
   {
      cout << "weight at i = " << i << ": " << weights[i] << endl;
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
      // trans->Transform(Geometries.GetCenter(pmesh->GetElementBaseGeometry(el_index)),cent);
      vdlfgi->AssembleRHSElementVect(*fe, *trans, elvec);
      // elvec.Print(cout);
      // cout << "cent: \n";
      // cent.Print(cout);
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