#include "mfem.hpp"
#include "geometry.hpp"
#include "lagrange_multiplier.hpp"
#include "var-config.h"
#include <cassert>

using namespace std;
using namespace mfem;
using namespace hydrodynamics;


// Problem
const char *mesh_file_location = "/data/ref-square.mesh";
std::string result = std::string(LAGLOS_DIR) + std::string(mesh_file_location);
const char* mesh_file = result.c_str();
static int num_procs = 0;
static int myid = -1;

int test_lmsparse_lmcop_getgradient();


int main(int argc, char *argv[])
{
   Mpi::Init();
   num_procs = Mpi::WorldSize();
   myid = Mpi::WorldRank();

   int d = 0;

   /* Test LM-Sparse LMCOp::GetGradient */
   d += test_lmsparse_lmcop_getgradient();

   return d;
}

/***
 * Function: test_lmsparse_lmcop_getgradient
 * Purpose: A basic test where an instance of LocalMassConservationOperator from the 
 * file lagrange_multipllier.hpp is instantiated and the GetGradient function is validated.
 */
int test_lmsparse_lmcop_getgradient()
{
   cout << "=== Testing LocalMassConservationOperator::GetGradient ===\n";

   mfem::Mesh *mesh;
   mesh = new mfem::Mesh(mesh_file, true, true);
   mesh->UniformRefinement();
   assert(mesh->GetNE() == 4);

   // Construct ParMesh
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;

   H1_FECollection H1FEC(2 /*order_mv*/, 2 /*dim*/);
   H1_FECollection H1FEC_L(1, 2);
   L2_FECollection L2FEC(0/*order*/, 2/*dim*/, BasisType::Positive);

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, 2/*dim*/);
   ParFiniteElementSpace H1FESpace_L(pmesh, &H1FEC_L, 2);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, 2);

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

   Geometric<2/*dim*/> geom(offset,H1FESpace);

   ParGridFunction x_gf, mv_gf, mv_gf_l(&H1FESpace_L);
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   pmesh->SetNodalGridFunction(&x_gf);
   mv_gf.MakeRef(&H1FESpace, S, offset[1]);
   mv_gf = 0.;

   Vector onex(2); onex[0] = 1., onex[1] = 0;
   geom.UpdateNodeVelocity(mv_gf, 4, onex);
   geom.UpdateNodeVelocity(mv_gf, 6, onex);
   geom.UpdateNodeVelocity(mv_gf, 8, onex);
   mv_gf_l.ProjectGridFunction(mv_gf);

   /* Set up sparsity pattern */
   Array<int> I, J;
   SetHiopConstraintGradSparsityPattern(pmesh, pmesh->GetNE(), H1FESpace_L.GetNVDofs(), I, J);

   LocalMassConservationOperator<2/*dim*/> LMCoper(geom, x_gf, pmesh->GetNE(), 2*H1FESpace_L.GetNVDofs(), 1/*dt*/, I, J);
   const Operator &oper_C = LMCoper.GetGradient(mv_gf_l);
   const SparseMatrix *grad_C = dynamic_cast<const SparseMatrix *>(&oper_C);
   const DenseMatrix *grad_C_dense = grad_C->ToDenseMatrix();
   cout << "grad_C_dense\n";
   grad_C_dense->Print(cout);

   cout << "x_gf: ";
   x_gf.Print(cout);

   // Run ComputeGradient with a set velocity vector and check the values returned to see they match up with a SparseMatrix we define.

   return 0;
}