// Copyright (c) 2022, Texas A&M University. All Rights
// reserved. See files LICENSE and NOTICE for details.
//
//                     __                __
//                    / /   ____  ____  / /___  _____
//                   / /   / __ `/ __ `/ / __ \/ ___/
//                  / /___/ /_/ / /_/ / / /_/ (__  )
//                 /_____/\__,_/\__, /_/\____/____/
//                             /____/
//
//             Low-order Lagrangian Hydrodynamics Miniapp
//
// Laglos(LAGrangian Low-Order Solver) is a miniapp that solves the
// time-dependent Euler equation of compressible gas dynamics in a moving
// Lagrangian frame using unstructured low-order finite element spatial
// discretization and forward euler time-stepping.

#include "mfem.hpp"

using namespace mfem;
using namespace std;

// Choice for the problem setup.
static int dim;

// Forward declarations
double fRand(double fMin, double fMax);
void mesh_v(const Vector &x, const double t, Vector &res);
void move_mesh(ParFiniteElementSpace & pfes);

int main(int argc, char *argv[]) {
   // Initialize MPI.
   MPI_Session mpi(argc, argv);
   const int myid = mpi.WorldRank();

   // Parse command line options
   dim = 2;
   const char *mesh_file = "default";
   int rs_levels = 0;
   int rp_levels = 0;
   int mesh_order = 1;
   int order_mv = 2;  // Order of mesh movement approximation space
   int order_u = 0;
   int precision = 12;

   OptionsParser args(argc, argv);
   args.AddOption(&dim, "-dim", "--dimension", "Dimension of the problem.");
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&rs_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&rp_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");

   args.Parse();
   if (!args.Good())
   {
      if (mpi.Root()) { args.PrintUsage(cout); }
      return 1;
   }
   if (mpi.Root()) { args.PrintOptions(cout); }

   // On all processors, use the default builtin 1D/2D/3D mesh or read the
   // serial one given on the command line.
   Mesh *mesh;
   if (strncmp(mesh_file, "default", 7) != 0)
   {
      mesh = new Mesh(mesh_file, true, true);
   }
   else // Default mesh
   {
      if (dim == 1)
      {
         mesh = new Mesh(Mesh::MakeCartesian1D(2));
         // mesh->GetBdrElement(0)->SetAttribute(1);
         // mesh->GetBdrElement(1)->SetAttribute(1);
      }
      if (dim == 2)
      {
         mesh = new Mesh(Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL,
                                               true));
         // const int NBE = mesh->GetNBE();
         // for (int b = 0; b < NBE; b++)
         // {
         //    Element *bel = mesh->GetBdrElement(b);
         //    const int attr = (b < NBE/2) ? 2 : 1;
         //    bel->SetAttribute(attr);
         // }
      }
      if (dim == 3)
      {
         mesh = new Mesh(Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON,
                                               true));
         // const int NBE = mesh->GetNBE();
         // for (int b = 0; b < NBE; b++)
         // {
         //    Element *bel = mesh->GetBdrElement(b);
         //    const int attr = (b < NBE/3) ? 3 : (b < 2*NBE/3) ? 1 : 2;
         //    bel->SetAttribute(attr);
         // }
      }
   }
   dim = mesh->Dimension(); // Set the correct dimension if a mesh other than
                            // 'default' was provided.
   // mesh->SetCurvature(2);

   // Refine the mesh in serial to increase the resolution. In this example
   // we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
   // a command-line parameter. If the mesh is of NURBS type, we convert it
   // to a (piecewise-polynomial) high-order mesh.
   for (int lev = 0; lev < rs_levels; lev++)
   {
      mesh->UniformRefinement();
   }
   // if (mesh->NURBSext)
   // {
   //    mesh->SetCurvature(max(mesh_order, 1));
   // }

   // Define the parallel mesh by a partitioning of the serial mesh. Refine
   // this mesh further in parallel to increase the resolution. Once the
   // parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   for (int lev = 0; lev < rp_levels; lev++)
   {
      pmesh->UniformRefinement();
   }

   int NE = pmesh->GetNE(), ne_min, ne_max;
   MPI_Reduce(&NE, &ne_min, 1, MPI_INT, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&NE, &ne_max, 1, MPI_INT, MPI_MAX, 0, pmesh->GetComm());

   if (myid == 0)
   { cout << "Zones min/max: " << ne_min << " " << ne_max << endl; }

   // Define the parallel finite element spaces. We use:
   // - H1 (Q2, continuous) for mesh movement.
   // - L2 (Q0, discontinuous) for state variables
   H1_FECollection H1FEC(order_mv, dim);
   L2_FECollection L2FEC(order_u, dim, BasisType::Positive);

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, pmesh->Dimension());
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);

   cout << "N dofs: " << H1FESpace.GetNDofs() << endl;
   cout << "Num Scalar Vertex DoFs: " << H1FESpace.GetNVDofs() << endl;
   cout << "Num Scalar edge-interior DoFs: " << H1FESpace.GetNEDofs() << endl;
   cout << "Num Scalar face-interior DoFs: " << H1FESpace.GetNFDofs() << endl;
   cout << "Num faces: " << H1FESpace.GetNF() << endl;
   cout << "dim: " << dim << endl;

   // Print initialized mesh
   // Can be visualized with glvis -np # -m mesh-test-init
   {
      ostringstream mesh_name;
      mesh_name << "mesh-test-init." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   /*
   * Move Mesh
   */
   move_mesh(H1FESpace);

   // Print moved mesh
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      ostringstream mesh_name;
      mesh_name << "mesh-test-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }
   
   // FaceGeometricFactors FaceGFs = *pmesh->GetFaceGeometricFactors();

}

/* 
*
* Mesh Movement
*
*/

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void mesh_v(const Vector &x, const double t, Vector &res)
{
   int mesh_opt = 1; // Change this param to get different movement

   switch (dim)
   {
      case 1:
      {
         res[0] = .1;
         return;
      }
      case 3:
      case 2:
      {
         switch (mesh_opt)
         {
            case 0: // constant velocity
            {
               for (int i = 0; i < dim; i++)
               {
                  res[i] = 0.2;
               }
               return;
            }
            case 1: // Move non corners
            {
               // Set v[1]
               if (x[0] != 0 && x[0] != 1)
               {
                  res[1] = 0.1;
               }
               else
               {
                  res[1] = 0.;
               }

               // Set v[0]
               if (x[1] != 0 && x[1] != 1)
               {
                  res[0] = .2;
               }
               else
               { 
                  res[0] = 0.;
               }
               return;
            }
            case 2: // rotation about origin with angle theta
            {
               const double theta = M_PI/4;
               res[0] = x[0]*(cos(theta) - 1) - x[1]*sin(theta);
               res[1] = x[0]*sin(theta) + x[1]*(cos(theta) - 1);
            }
            case 3: // random
            {
               const double r_max = .125;
               res[0] = fRand(-r_max, r_max);
               res[1] = fRand(-r_max, r_max);
            }
         } // mesh_opt switch
      } // 2D velocity end
   }
}

void move_mesh(ParFiniteElementSpace & pfes)
{
   ParGridFunction *x_gf = new ParGridFunction(&pfes);
   VectorFunctionCoefficient mesh_velocity(dim, mesh_v);
   mesh_velocity.SetTime(0);
   ParGridFunction *v_gf = new ParGridFunction(&pfes);
   v_gf->ProjectCoefficient(mesh_velocity);

   // Initialize x_gf using the starting mesh coordinates.
   pfes.GetParMesh()->SetNodalGridFunction(x_gf);

   // Logic to change x_gf based on velocities in v_gf
   // V = v0 + 0.5 * dt * dv_dt;
   add(*x_gf, 1, *v_gf, *x_gf);

   // Update mesh gridfunction
   pfes.GetParMesh()->NewNodes(*x_gf, false);
}