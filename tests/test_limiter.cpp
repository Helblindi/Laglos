#include "mfem.hpp"
#include "laglos_solver.hpp"
#include "limiter.h"
#include "test_problems_include.h"
#include "var-config.h"
#include <cassert>

using namespace std;
using namespace mfem;
using namespace hydroLO;

/* Parameters */
int dim = 2;
string result = string(LAGLOS_DIR) + "data/ref-square.mesh";
const char* mesh_file_location = result.c_str();

static int num_procs = 0;
static int myid = -1;
/* ---------------- End Parameters ---------------- */

/* Forward declarations */
int TestComputeLocalBoundsOT1OK2();
int TestComputeLocalBoundsOT2OK3();

int main(int argc, char *argv[])
{
   // Initialize MPI.
   Mpi::Init();
   num_procs = Mpi::WorldSize();
   myid = Mpi::WorldRank();

   int d = 0;
   d += TestComputeLocalBoundsOT1OK2();
   d += TestComputeLocalBoundsOT2OK3();

   return d;
}

/**
 * @brief Unit test for the IDPLimiter::ComputeLocalBounds method with order_t = 1 and order_k = 2.
 *
 * This test sets up a finite element mesh and associated spaces, constructs a low-order grid function,
 * and uses the IDPLimiter to compute the minimum and maximum values of the grid function over the mesh.
 * The computed min/max values are compared against exact reference values to validate correctness.
 *
 * The high order space:  Min:           Max: 
 *    *-----*-----*      *-----*-----*   *-----*-----*
 *    |14 15|10 11|      |12  8|8   8|   |15 15|15 11|
 *    |     |     |      |     |     |   |     |     |
 *    |12 13|8   9|      |2   2|3   6|   |15 15|15 11|
 *    *-----*-----*      *-----*-----*   *-----*-----*
 *    |2   3|6   7|      |0   0|1   4|   |13 13|13  9|
 *    |     |     |      |     |     |   |     |     |
 *    |0   1|4   5|      |0   0|1   4|   |3   6|7   7|
 *    *-----*-----*      *-----*-----*   *-----*-----*
 *
 * Steps performed:
 * - Initializes mesh and refines it.
 * - Sets up high-order and low-order finite element spaces.
 * - Constructs a constant mass vector.
 * - Initializes an IDPLimiter object.
 * - Fills a low-order grid function with test values.
 * - Computes min and max values using the limiter.
 * - Compares results to exact reference vectors.
 * - Reports errors and returns 1 on failure, 0 on success.
 *
 * @return int Returns 0 if the test passes, 1 if any min/max errors are detected.
 */
int TestComputeLocalBoundsOT1OK2()
{
   cout << "========================================================================\n"
        << "\tTestComputeLocalBoundsOT1OK2\n";

   /* Test specific parameters */
   int order_t = 1, order_k = 2;
   int rs_levels = 1;

   // Create mesh objects
   Mesh *mesh = new Mesh(mesh_file_location, true, true);
   for (int lev = 0; lev < rs_levels; lev++)
   {
      mesh->UniformRefinement();
   }
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   ParMesh *pmesh_lo = new ParMesh(ParMesh::MakeRefined(*pmesh, order_t+1, BasisType::ClosedUniform));

   /* Define FE spaces */
   L2_FECollection L2FEC_HO(order_t, dim, BasisType::Positive);
   H1_FECollection H1FEC_HO(order_k, dim);
   L2_FECollection L2FEC_LO(0, dim, BasisType::Positive);
   H1_FECollection H1FEC_LO(1, dim);

   ParFiniteElementSpace L2FESpace_HO(pmesh, &L2FEC_HO),
                         H1FESpace_HO(pmesh, &H1FEC_HO, dim),
                         L2FESpace_LO(pmesh_lo, &L2FEC_LO),
                         H1FESpace_LO(pmesh_lo, &H1FEC_LO, dim);
   ParFiniteElementSpace H1FESpace_proj_HO(pmesh, &H1FEC_HO, 1),
                         H1FESpace_proj_LO(pmesh_lo, &H1FEC_LO, 1);

   L2FESpace_HO.BuildDofToArrays();

   /* Construct mass vector */
   ConstantCoefficient one_coeff(1.0);
   ParLinearForm *mHO = new ParLinearForm(&L2FESpace_HO);
   mHO->AddDomainIntegrator(new DomainLFIntegrator(one_coeff));
   mHO->Assemble();
   HypreParVector *mHO_hpv = mHO->ParallelAssemble();

   /* Define IDPLimiter object */
   int order_q = 4; // Order of quadrature
   IDPLimiter *idpl = new IDPLimiter(L2FESpace_HO, L2FESpace_LO, H1FESpace_proj_LO, H1FESpace_proj_HO, *mHO_hpv, order_q);

   /* Construct LO grid function to toss into compute rho min/max */
   ParGridFunction rho_gf_LO(&L2FESpace_LO);
   for (int i = 0; i < L2FESpace_LO.GetNDofs(); i++)
   {
      rho_gf_LO(i) = (double)i; // Initialize with some value
   }
   idpl->ComputeLocalBounds(rho_gf_LO);
   
   /* Set exact min and max */
   // Take all neighbors
   double xmin_vec[16] = {0., 0., 0., 0., 1., 4., 1., 4., 3., 6., 8., 8., 2., 2., 12., 8.};
   double xmax_vec[16] = {3., 6., 13., 13., 7., 7., 13., 9., 15., 11., 15., 11., 15., 15., 15., 15.};
   const Vector xmin_ex(xmin_vec);
   const Vector xmax_ex(xmax_vec);

   /* Compute min and max*/
   Vector xmin, xmax;
   idpl->GetRhoMin(xmin);
   idpl->GetRhoMax(xmax);

   /* Validate */
   int num_min_error = 0, num_max_error = 0;
   for (int i = 0; i < L2FESpace_LO.GetNDofs(); i++)
   {
      if (xmin_ex[i] != xmin[i])
      {
         cout << "xmin_ex[" << i << "] = " << xmin_ex[i] << ", xmin[" << i << "] = " << xmin[i] << endl;
         num_min_error++;
      }
      if (xmax_ex[i] != xmax[i])
      {
         cout << "xmax_ex[" << i << "] = " << xmax_ex[i] << ", xmax[" << i << "] = " << xmax[i] << endl;
         num_max_error++;
      }
   }

   /* Clean up objects */
   delete pmesh;
   delete pmesh_lo;
   delete idpl;

   if (num_min_error > 0 || num_max_error > 0)
   {
      cout << "TestComputeLocalBoundsOT1OK2 failed with " 
           << num_min_error << " min errors and "
           << num_max_error << " max errors." << endl
           << "========================================================================\n\n";
      return 1;
   }
   else
   {
      cout << "TestComputeLocalBoundsOT1OK2 passed!" << endl
           << "========================================================================\n\n";
   }
        
   return 0;
}

int TestComputeLocalBoundsOT2OK3()
{
   cout << "========================================================================\n"
        << "\tTestComputeLocalBoundsOT2OK3\n";

   /* Test specific parameters */
   int order_t = 2, order_k = 3;
   int rs_levels = 1;

   // Create mesh objects
   Mesh *mesh = new Mesh(mesh_file_location, true, true);
   for (int lev = 0; lev < rs_levels; lev++)
   {
      mesh->UniformRefinement();
   }
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   ParMesh *pmesh_lo = new ParMesh(ParMesh::MakeRefined(*pmesh, order_t+1, BasisType::ClosedUniform));

   /* Define FE spaces */
   L2_FECollection L2FEC_HO(order_t, dim, BasisType::Positive);
   H1_FECollection H1FEC_HO(order_k, dim);
   L2_FECollection L2FEC_LO(0, dim, BasisType::Positive);
   H1_FECollection H1FEC_LO(1, dim);

   ParFiniteElementSpace L2FESpace_HO(pmesh, &L2FEC_HO),
                         H1FESpace_HO(pmesh, &H1FEC_HO, dim),
                         L2FESpace_LO(pmesh_lo, &L2FEC_LO),
                         H1FESpace_LO(pmesh_lo, &H1FEC_LO, dim);
   ParFiniteElementSpace H1FESpace_proj_HO(pmesh, &H1FEC_HO, 1),
                         H1FESpace_proj_LO(pmesh_lo, &H1FEC_LO, 1);

   L2FESpace_HO.BuildDofToArrays();

   /* Construct mass vector */
   ConstantCoefficient one_coeff(1.0);
   ParLinearForm *mHO = new ParLinearForm(&L2FESpace_HO);
   mHO->AddDomainIntegrator(new DomainLFIntegrator(one_coeff));
   mHO->Assemble();
   HypreParVector *mHO_hpv = mHO->ParallelAssemble();

   /* Define IDPLimiter object */
   int order_q = 4; // Order of quadrature
   IDPLimiter *idpl = new IDPLimiter(L2FESpace_HO, L2FESpace_LO, H1FESpace_proj_LO, H1FESpace_proj_HO, *mHO_hpv, order_q);

   /* Construct LO grid function to toss into compute rho min/max */
   ParGridFunction rho_gf_LO(&L2FESpace_LO);
   for (int i = 0; i < L2FESpace_LO.GetNDofs(); i++)
   {
      rho_gf_LO(i) = (double)i; // Initialize with some value
   }
   idpl->ComputeLocalBounds(rho_gf_LO);
   
   /* Set exact min and max */
   // Full adjacency, not using cell adjacency
   double xmin_vec[36] = {0., 0., 1., 0., 0., 1., 3., 3., 4., 2., 9., 10., 2., 9., 10., 5., 12., 13., 8., 15., 16., 18., 18., 19., 21., 21., 22., 6., 6., 7., 27., 27., 18., 30., 30., 21.};
   double xmax_vec[36] = {4., 5., 12., 7., 8., 15., 28., 29., 29., 13., 14., 14., 16., 17., 17., 29., 20., 20., 32., 23., 23., 35., 26., 26., 35., 26., 26., 31., 32., 32., 34., 35., 35., 34., 35., 35.};
   const Vector xmin_ex(xmin_vec);
   const Vector xmax_ex(xmax_vec);

   /* Compute min and max*/
   Vector xmin, xmax;
   idpl->GetRhoMin(xmin);
   idpl->GetRhoMax(xmax);

   /* Validate */
   int num_min_error = 0, num_max_error = 0;
   for (int i = 0; i < L2FESpace_LO.GetNDofs(); i++)
   {
      if (xmin_ex[i] != xmin[i])
      {
         cout << "xmin_ex[" << i << "] = " << xmin_ex[i] << ", xmin[" << i << "] = " << xmin[i] << endl;
         num_min_error++;
      }
      if (xmax_ex[i] != xmax[i])
      {
         cout << "xmax_ex[" << i << "] = " << xmax_ex[i] << ", xmax[" << i << "] = " << xmax[i] << endl;
         num_max_error++;
      }
   }

   /* Clean up objects */
   delete pmesh;
   delete pmesh_lo;
   delete idpl;

   if (num_min_error > 0 || num_max_error > 0)
   {
      cout << "TestComputeLocalBoundsOT1OK2 failed with " 
           << num_min_error << " min errors and "
           << num_max_error << " max errors." << endl
           << "========================================================================\n\n";
      return 1;
   }
   else
   {
      cout << "TestComputeLocalBoundsOT1OK2 passed!" << endl
           << "========================================================================\n\n";
   }
        
   return 0;
}