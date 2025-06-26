#include "mfem.hpp"
#include "laglos_solver.hpp"
#include "test_problems_include.h"
#include "var-config.h"
#include <cassert>

using namespace std;
using namespace mfem;
using namespace hydroLO;

/* Parameters */
int dim = 2;
int order_mv = 2;          // Order of mesh movement approximation space
int order_u = 0;
string result = string(LAGLOS_DIR) + "data/ref-square.mesh";
const char* mesh_file_location = result.c_str();

int mesh_refinements = 1;
bool use_viscosity = true; // Doesn't matter
const bool _mm = true;          // Doesn't matter
int mv_option = 2;
int fv_option = 2;
double CFL = 0.5;          // Doesn't matter
int rs_levels = 2;
double t_init = 0.0;
double t_final = 0.225;
const int elastic_eos = 0;

static int num_procs = 0;
static int myid = -1;
/* ---------------- End Parameters ---------------- */

int CompareCijComputation();
void TestU2();
int ValidateCijComputationOrder1();
int TestDGNormalIntegrator();

int main(int argc, char *argv[])
{
   // Initialize MPI.
   Mpi::Init();
   num_procs = Mpi::WorldSize();
   myid = Mpi::WorldRank();

   int d = 0;
   d += CompareCijComputation();
   d += ValidateCijComputationOrder1();
   d += TestDGNormalIntegrator();
   // TestU2();

   return d;
}

int CompareCijComputation()
{
   // On all processors, use the default builtin 1D/2D/3D mesh or read the
   // serial one given on the command line.
   Mesh *mesh = new Mesh(mesh_file_location, true, true);
   dim = mesh->Dimension();

   // Refine the mesh in serial to increase the resolution. In this example
   // we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
   // a command-line parameter. If the mesh is of NURBS type, we convert it
   // to a (piecewise-polynomial) high-order mesh.
   for (int lev = 0; lev < rs_levels; lev++)
   {
      mesh->UniformRefinement();
   }

   // Define the parallel mesh by a partitioning of the serial mesh. Refine
   // this mesh further in parallel to increase the resolution. Once the
   // parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;

   int NE = pmesh->GetNE(), ne_min, ne_max;
   MPI_Reduce(&NE, &ne_min, 1, MPI_INT, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&NE, &ne_max, 1, MPI_INT, MPI_MAX, 0, pmesh->GetComm());
   double hmin, hmax, kmin, kmax;
   pmesh->GetCharacteristics(hmin, hmax, kmin, kmax);

   // Set up problem
   ProblemBase * problem_class = new SodProblem(dim);

   // Define the parallel finite element spaces. We use:
   // - H1 (Q2, continuous) for mesh movement.
   // - L2 (Q0, discontinuous) for state variables
   // - CR/RT for mesh reconstruction at nodes
   H1_FECollection H1FEC(order_mv, dim);
   H1_FECollection H1FEC_L(1, dim);
   L2_FECollection L2FEC(order_u, dim, BasisType::Positive);
   FiniteElementCollection * CRFEC;
   if (dim == 1)
   {
      CRFEC = new CrouzeixRaviartFECollection();
   }
   else
   {
      CRFEC = new RT_FECollection(0, dim);
   }

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace H1FESpace_L(pmesh, &H1FEC_L, dim);
   /* Finite element space solely constructed for continuous representation of density field */
   ParFiniteElementSpace H1cFESpace(pmesh, &H1FEC, 1);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);
   ParFiniteElementSpace CRFESpace(pmesh, CRFEC, dim);

   /* The monolithic BlockVector stores unknown fields as:
   *   - 0 -> position
   *   - 1 -> specific volume
   *   - 2 -> velocity (L2V)
   *   - 3 -> speific total energy
   */
   const int Vsize_l2 = L2FESpace.GetVSize();
   const int Vsize_l2v = L2VFESpace.GetVSize();
   const int Vsize_h1 = H1FESpace.GetVSize();
   Array<int> offset(5);
   offset[0] = 0;
   offset[1] = offset[0] + Vsize_h1;
   offset[2] = offset[1] + Vsize_l2;
   offset[3] = offset[2] + Vsize_l2v;
   offset[4] = offset[3] + Vsize_l2;
   BlockVector S(offset, Device::GetMemoryType());

   // Define GridFunction objects for the position, mesh velocity and specific
   // volume, velocity, and specific internal energy. At each step, each of 
   // these values will be updated.
   ParGridFunction x_gf, sv_gf, v_gf, ste_gf;
   ParGridFunction rho0_gf(&L2FESpace);
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   sv_gf.MakeRef(&L2FESpace, S, offset[1]);
   v_gf.MakeRef(&L2VFESpace, S, offset[2]);
   ste_gf.MakeRef(&L2FESpace, S, offset[3]);

   // Initialize x_gf using starting mesh positions
   pmesh->SetNodalGridFunction(&x_gf);
   // Sync the data location of x_gf with its base, S
   x_gf.SyncAliasMemory(S);

   // Change class variables into static std::functions since virtual static member functions are not an option
   // and Coefficient class requires std::function arguments
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> sv0_static = 
      std::bind(&ProblemBase::sv0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<void(const Vector &, const double, Vector &)> v0_static = 
      std::bind(&ProblemBase::v0, problem_class, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
   std::function<double(const Vector &,const double)> ste0_static = 
      std::bind(&ProblemBase::ste0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> rho0_static = 
      std::bind(&ProblemBase::rho0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> p0_static = 
      std::bind(&ProblemBase::p0, problem_class, std::placeholders::_1, std::placeholders::_2);

   // Initialize specific volume, velocity, and specific total energy
   FunctionCoefficient sv_coeff(sv0_static);
   sv_coeff.SetTime(t_init);
   sv_gf.ProjectCoefficient(sv_coeff);
   sv_gf.SyncAliasMemory(S);

   VectorFunctionCoefficient v_coeff(dim, v0_static);
   v_coeff.SetTime(t_init);
   v_gf.ProjectCoefficient(v_coeff);
   v_gf.SyncAliasMemory(S);

   // While the ste_coeff is not used for initialization in the Sedov case,
   // it is necessary for plotting the exact solution
   FunctionCoefficient ste_coeff(ste0_static);
   ste_coeff.SetTime(t_init);
   ste_gf.ProjectCoefficient(ste_coeff);
   ste_gf.SyncAliasMemory(S);

   // PLF to build mass vector
   FunctionCoefficient rho_coeff(rho0_static); 
   rho_coeff.SetTime(t_init);
   IntegrationRule ir = IntRules.Get(pmesh->GetElementBaseGeometry(0), 3 * H1FESpace.GetOrder(0) + L2FESpace.GetOrder(0) - 1);
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff, &ir));
   m->Assemble();

   /* Initialize rho0_gf */
   rho0_gf.ProjectCoefficient(rho_coeff);

   FunctionCoefficient p_coeff(p0_static);
   p_coeff.SetTime(t_init);

   /* Create Lagrangian Low Order Solver Object */
   LagrangianLOOperator hydro(dim, S.Size(), H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, rho_coeff, rho0_gf, m, ir, problem_class, offset, use_viscosity, elastic_eos, _mm, CFL);
   hydro.BuildCijMatrices();
   // assert(false);
   Array<int> fids, oris;
   mfem::Mesh::FaceInformation FI;
   Vector n_int(dim), c(dim), cij(dim), _t(dim);
   int cj = 0;
   for (int ci = 0; ci < NE; ci++)
   {
      pmesh->GetElementEdges(ci, fids, oris);
      for (int j=0; j < fids.Size(); j++) // Face iterator
      {
         FI = pmesh->GetFaceInformation(fids[j]);
         if (FI.IsInterior())
         {
            // Get index information/state vector for second cell
            if (ci == FI.element[0].index) { 
               cj = FI.element[1].index; 
            }
            else { 
               assert(ci == FI.element[1].index);
               cj = FI.element[0].index; 
            }

            hydro.CalcOutwardNormalInt(S, ci, fids[j], n_int);
            c = n_int;
            c /= 2.;

            hydro.GetLocalCij(ci, cj, cij);
            
            subtract(c, cij, _t);
            if (_t.Norml2() > 1.e-10)
            {
               cout << "cij computation does not match outward normal!\n";
               cout << "i: " << ci << ", j: " << cj << endl;
               cout << "Outward normal: ";
               c.Print(cout);
               cout << "cij: ";
               cij.Print(cout);
               return 1;
            }

         }
      }
   }

   /* free memory */
   delete problem_class;
   delete CRFEC;
   delete m;
   delete pmesh;

   return 0;
}


void TestU2()
{
   cout << "Testing U2 problem..." << endl;
   // On all processors, use the default builtin 1D/2D/3D mesh or read the
   // serial one given on the command line.
   Mesh *mesh = new Mesh(mesh_file_location, true, true);
   dim = mesh->Dimension();

   // Refine the mesh in serial to increase the resolution. In this example
   // we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
   // a command-line parameter. If the mesh is of NURBS type, we convert it
   // to a (piecewise-polynomial) high-order mesh.
   for (int lev = 0; lev < 1; lev++)
   {
      mesh->UniformRefinement();
   }

   // Define the parallel mesh by a partitioning of the serial mesh. Refine
   // this mesh further in parallel to increase the resolution. Once the
   // parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;

   int NE = pmesh->GetNE(), ne_min, ne_max;
   MPI_Reduce(&NE, &ne_min, 1, MPI_INT, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&NE, &ne_max, 1, MPI_INT, MPI_MAX, 0, pmesh->GetComm());
   double hmin, hmax, kmin, kmax;
   pmesh->GetCharacteristics(hmin, hmax, kmin, kmax);

   // Set up problem
   ProblemBase * problem_class = new SodProblem(dim);

   // Define the parallel finite element spaces. We use:
   // - H1 (Q2, continuous) for mesh movement.
   // - L2 (Q0, discontinuous) for state variables
   // - CR/RT for mesh reconstruction at nodes
   H1_FECollection H1FEC(order_mv, dim);
   H1_FECollection H1FEC_L(1, dim);
   L2_FECollection L2FEC(1, dim, BasisType::Positive);
   FiniteElementCollection * CRFEC;
   if (dim == 1)
   {
      CRFEC = new CrouzeixRaviartFECollection();
   }
   else
   {
      CRFEC = new RT_FECollection(0, dim);
   }

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace H1FESpace_L(pmesh, &H1FEC_L, dim);
   /* Finite element space solely constructed for continuous representation of density field */
   ParFiniteElementSpace H1cFESpace(pmesh, &H1FEC, 1);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);
   ParFiniteElementSpace CRFESpace(pmesh, CRFEC, dim);

   /* The monolithic BlockVector stores unknown fields as:
   *   - 0 -> position
   *   - 1 -> specific volume
   *   - 2 -> velocity (L2V)
   *   - 3 -> speific total energy
   */
   const int Vsize_l2 = L2FESpace.GetVSize();
   const int Vsize_l2v = L2VFESpace.GetVSize();
   const int Vsize_h1 = H1FESpace.GetVSize();
   Array<int> offset(5);
   offset[0] = 0;
   offset[1] = offset[0] + Vsize_h1;
   offset[2] = offset[1] + Vsize_l2;
   offset[3] = offset[2] + Vsize_l2v;
   offset[4] = offset[3] + Vsize_l2;
   BlockVector S(offset, Device::GetMemoryType());

   // Define GridFunction objects for the position, mesh velocity and specific
   // volume, velocity, and specific internal energy. At each step, each of 
   // these values will be updated.
   ParGridFunction x_gf, sv_gf, v_gf, ste_gf;
   ParGridFunction rho0_gf(&L2FESpace);
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   sv_gf.MakeRef(&L2FESpace, S, offset[1]);
   v_gf.MakeRef(&L2VFESpace, S, offset[2]);
   ste_gf.MakeRef(&L2FESpace, S, offset[3]);

   // Initialize x_gf using starting mesh positions
   pmesh->SetNodalGridFunction(&x_gf);
   // Sync the data location of x_gf with its base, S
   x_gf.SyncAliasMemory(S);

   // Change class variables into static std::functions since virtual static member functions are not an option
   // and Coefficient class requires std::function arguments
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> sv0_static = 
      std::bind(&ProblemBase::sv0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<void(const Vector &, const double, Vector &)> v0_static = 
      std::bind(&ProblemBase::v0, problem_class, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
   std::function<double(const Vector &,const double)> ste0_static = 
      std::bind(&ProblemBase::ste0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> rho0_static = 
      std::bind(&ProblemBase::rho0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> p0_static = 
      std::bind(&ProblemBase::p0, problem_class, std::placeholders::_1, std::placeholders::_2);

   // Initialize specific volume, velocity, and specific total energy
   FunctionCoefficient sv_coeff(sv0_static);
   sv_coeff.SetTime(t_init);
   sv_gf.ProjectCoefficient(sv_coeff);

   VectorFunctionCoefficient v_coeff(dim, v0_static);
   v_coeff.SetTime(t_init);
   v_gf.ProjectCoefficient(v_coeff);

   // While the ste_coeff is not used for initialization in the Sedov case,
   // it is necessary for plotting the exact solution
   FunctionCoefficient ste_coeff(ste0_static);
   ste_coeff.SetTime(t_init);
   ste_gf.ProjectCoefficient(ste_coeff);


   // PLF to build mass vector
   FunctionCoefficient rho_coeff(rho0_static); 
   rho_coeff.SetTime(t_init);
   IntegrationRule ir = IntRules.Get(pmesh->GetElementBaseGeometry(0), 3 * H1FESpace.GetOrder(0) + L2FESpace.GetOrder(0) - 1);
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff, &ir));
   m->Assemble();

   /* Initialize rho0_gf */
   rho0_gf.ProjectCoefficient(rho_coeff);

   FunctionCoefficient p_coeff(p0_static);
   p_coeff.SetTime(t_init);

   /* Create Lagrangian Low Order Solver Object */
   LagrangianLOOperator hydro(dim, S.Size(), H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, rho_coeff, rho0_gf, m, ir, problem_class, offset, use_viscosity, elastic_eos, _mm, CFL);
   hydro.BuildCijMatrices();

   ParGridFunction L2_coords(&L2VFESpace);
   pmesh->GetNodes(L2_coords);
   cout << "L2 coordinates:\n";
   L2_coords.Print(cout);

   /* free memory */
   cout << "freeing memory..." << endl;
   delete problem_class;
   delete CRFEC;
   delete m;
   delete pmesh;
   cout << "done." << endl;
}

/*
The values of Cij should match the definition given in equation (4.17) in the paper:
Guermond, Popov, Tomas: 
"Invariant domain preserving discretization-independent schemes and convex limiting for hyperbolic systems."
*/
int ValidateCijComputationOrder1()
{
   // This function is a placeholder for validating the Cij computation order 1.
   // It can be implemented with specific tests or assertions as needed.
   cout << "Validating Cij computation order 1..." << endl;
   Mesh *mesh = new Mesh(Mesh::MakeCartesian2D(1, 1, Element::QUADRILATERAL, true));
   dim = mesh->Dimension();
   /* Optional serial refinement */
   int num_refinements = 1;
   for (int lev = 0; lev < num_refinements; lev++)
   {
      mesh->UniformRefinement();
   }
   
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;

   // Set up problem
   ProblemBase * problem_class = new SodProblem(dim);

   // Define the parallel finite element spaces. We use:
   // - H1 (Q2, continuous) for mesh movement.
   // - L2 (Q0, discontinuous) for state variables
   // - CR/RT for mesh reconstruction at nodes
   H1_FECollection H1FEC(1, dim);
   H1_FECollection H1FEC_L(1, dim);
   L2_FECollection L2FEC(1, dim, BasisType::Positive);
   FiniteElementCollection * CRFEC;
   if (dim == 1)
   {
      CRFEC = new CrouzeixRaviartFECollection();
   }
   else
   {
      CRFEC = new RT_FECollection(0, dim);
   }

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace H1FESpace_L(pmesh, &H1FEC_L, dim);
   /* Finite element space solely constructed for continuous representation of density field */
   ParFiniteElementSpace H1cFESpace(pmesh, &H1FEC, 1);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);
   ParFiniteElementSpace CRFESpace(pmesh, CRFEC, dim);
   const int NDofs_L2 = L2FESpace.GetNDofs();

   /* The monolithic BlockVector stores unknown fields as:
   *   - 0 -> position
   *   - 1 -> specific volume
   *   - 2 -> velocity (L2V)
   *   - 3 -> speific total energy
   */
   const int Vsize_l2 = L2FESpace.GetVSize();
   const int Vsize_l2v = L2VFESpace.GetVSize();
   const int Vsize_h1 = H1FESpace.GetVSize();
   Array<int> offset(5);
   offset[0] = 0;
   offset[1] = offset[0] + Vsize_h1;
   offset[2] = offset[1] + Vsize_l2;
   offset[3] = offset[2] + Vsize_l2v;
   offset[4] = offset[3] + Vsize_l2;
   BlockVector S(offset, Device::GetMemoryType());

   // Define GridFunction objects for the position, mesh velocity and specific
   // volume, velocity, and specific internal energy. At each step, each of 
   // these values will be updated.
   ParGridFunction x_gf, sv_gf, v_gf, ste_gf;
   ParGridFunction rho0_gf(&L2FESpace);
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   sv_gf.MakeRef(&L2FESpace, S, offset[1]);
   v_gf.MakeRef(&L2VFESpace, S, offset[2]);
   ste_gf.MakeRef(&L2FESpace, S, offset[3]);

   // Initialize x_gf using starting mesh positions
   pmesh->SetNodalGridFunction(&x_gf);
   // Sync the data location of x_gf with its base, S
   x_gf.SyncAliasMemory(S);

   // Change class variables into static std::functions since virtual static member functions are not an option
   // and Coefficient class requires std::function arguments
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> sv0_static = 
      std::bind(&ProblemBase::sv0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<void(const Vector &, const double, Vector &)> v0_static = 
      std::bind(&ProblemBase::v0, problem_class, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
   std::function<double(const Vector &,const double)> ste0_static = 
      std::bind(&ProblemBase::ste0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> rho0_static = 
      std::bind(&ProblemBase::rho0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> p0_static = 
      std::bind(&ProblemBase::p0, problem_class, std::placeholders::_1, std::placeholders::_2);

   // Initialize specific volume, velocity, and specific total energy
   FunctionCoefficient sv_coeff(sv0_static);
   sv_coeff.SetTime(t_init);
   sv_gf.ProjectCoefficient(sv_coeff);

   VectorFunctionCoefficient v_coeff(dim, v0_static);
   v_coeff.SetTime(t_init);
   v_gf.ProjectCoefficient(v_coeff);

   // While the ste_coeff is not used for initialization in the Sedov case,
   // it is necessary for plotting the exact solution
   FunctionCoefficient ste_coeff(ste0_static);
   ste_coeff.SetTime(t_init);
   ste_gf.ProjectCoefficient(ste_coeff);


   // PLF to build mass vector
   FunctionCoefficient rho_coeff(rho0_static); 
   rho_coeff.SetTime(t_init);
   IntegrationRule ir = IntRules.Get(pmesh->GetElementBaseGeometry(0), 3 * H1FESpace.GetOrder(0) + L2FESpace.GetOrder(0) - 1);
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff, &ir));
   m->Assemble();

   /* Initialize rho0_gf */
   rho0_gf.ProjectCoefficient(rho_coeff);

   FunctionCoefficient p_coeff(p0_static);
   p_coeff.SetTime(t_init);

   /* Create Lagrangian Low Order Solver Object */
   LagrangianLOOperator hydro(dim, S.Size(), H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, rho_coeff, rho0_gf, m, ir, problem_class, offset, use_viscosity, elastic_eos, _mm, CFL);
   hydro.BuildCijMatrices();
   Table L2Connectivity;
   hydro.GetL2ConnectivityTable(L2Connectivity);

   Vector cij(dim), cji(dim), _t(dim), sum(dim);;
   for (int i = 0; i < NDofs_L2; i++)
   {
      MFEM_WARNING("May want to check row sums for interior dofs.");
      /*
      I have manually checked that row sums for interior dofs is 0
      for smaller meshes.
      In order to do this at scaled I would need to rebase my fork of
      mfem to get access to the new function
      FiniteElementSpace::GetExteriorTrueDofs()
      */
      // sum = 0.;
      // for (int j = 0; j < NDofs_L2; j++)
      // {
      //    hydro.GetLocalCij(i,j,cij);
      //    cout << "Cij[" << i << "][" << j << "] = ";
      //    cij.Print(cout);
      //    sum.Add(1., cij);
      // }
      // cout << "Row " << i << " sum: " << sum.Norml2() << endl;

      /* Check that resulting matrices are skew symmetric */
      for (int j = i+1; j < NDofs_L2; j++)
      {
         hydro.GetLocalCij(i, j, cij);
         hydro.GetLocalCij(j,i,cji);
         add(cij, cji, _t);
         if (_t.Norml2() > 1.e-10)
         {
            cout << "Cij is not skew symmetric!\n";
            cout << "i: " << i << ", j: " << j << endl;
            cout << "cij: ";
            cij.Print(cout);
            cout << "cji: ";
            cji.Print(cout);
            return 1;
         }
      }
   }

   delete pmesh;
   return 0;
}

int TestDGNormalIntegrator()
{
   // This function is a placeholder for testing the DGNormalIntegrator.
   // It can be implemented with specific tests or assertions as needed.
   cout << "Testing DGNormalIntegrator..." << endl;

   Mesh *mesh = new Mesh(Mesh::MakeCartesian2D(1, 1, Element::QUADRILATERAL, true));
   dim = mesh->Dimension();
   /* Optional serial refinement */
   int num_refinements = 0;
   for (int lev = 0; lev < num_refinements; lev++)
   {
      mesh->UniformRefinement();
   }

   L2_FECollection L2FEC(1, dim, BasisType::Positive);
   FiniteElementSpace L2(mesh, &L2FEC);
   const int NDofs_L2 = L2.GetNDofs();
   BilinearForm xbf(&L2), ybf(&L2);
   ConstantCoefficient one(1.0);

   // xbf.AddDomainIntegrator(new DerivativeIntegrator(one, 0));
   xbf.AddInteriorFaceIntegrator(new DGNormalIntegrator(-1., 0));
   xbf.AddBdrFaceIntegrator(new DGNormalIntegrator(-1., 0));
   xbf.Assemble();
   xbf.Finalize();
   SparseMatrix xmat_sparse = xbf.SpMat();

   // ybf.AddDomainIntegrator(new DerivativeIntegrator(one, 1));
   ybf.AddInteriorFaceIntegrator(new DGNormalIntegrator(-1., 1));
   ybf.AddBdrFaceIntegrator(new DGNormalIntegrator(-1., 1));
   ybf.Assemble();
   ybf.Finalize();
   SparseMatrix ymat_sparse = ybf.SpMat();
   
   // cout << "xmat_sparse:\n";
   // xmat_sparse.Print(cout);
   // cout << "ymat_sparse:\n";
   // ymat_sparse.Print(cout);

   // Array<int> cols;
   // Vector vals, cij(dim);;
   // int col_index;

   /* Get cij */
   // for (int i = 0; i < NDofs_L2; i++)
   // {
   //    for (int j = 0; j < NDofs_L2; j++)
   //    {
   //       cij = 0.;
   //       xmat_sparse.GetRow(i, cols, vals);
   //       col_index = cols.Find(j);
   //       if (col_index != -1) {
   //          cij[0] = vals[col_index];
   //       } else {
   //          // cout << "col not found\n";
   //          cij[0] = 0.;
   //       }

   //       ymat_sparse.GetRow(i, cols, vals);
   //       col_index = cols.Find(j);
   //       if (col_index != -1) {
   //          cij[1] = vals[col_index];
   //       } else {
   //          // cout << "col not found for i: " << i << ", j: " << j << "\n";
   //          cij[1] = 0.;
   //       }
   //       cout << "Cij[" << i << "][" << j << "] = ";
   //       cij.Print(cout);
   //    }
   // }

   /* Free memory */
   delete mesh;

   /* Check symmetry */
   if (xmat_sparse.IsSymmetric() < 1.E-10 && ymat_sparse.IsSymmetric() < 1.E-10)
   {
      cout << "Both matrices are symmetric.\n";
      return 0;
   }
   else
   {
      cout << "At least one matrix is not symmetric.\n";
      cout << "symm val x: " << xmat_sparse.IsSymmetric() << endl;
      cout << "symm val y: " << ymat_sparse.IsSymmetric() << endl;
      cout << "xmat_sparse:\n";
      xmat_sparse.Print(cout);
      cout << "ymat_sparse:\n";
      ymat_sparse.Print(cout);
      return 1;
   }

   return 0;
}