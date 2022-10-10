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
void orthogonal(Vector &v);
void build_C(ParFiniteElementSpace &pfes, ParGridFunction & velocities, double dt);
void calc_outward_normal_int(int cell, int face, ParFiniteElementSpace &pfes, ParGridFunction & velocities, double dt, Vector & res);
Vector Get_Int_Der_Ref_Shap_Functions();
double fRand(double fMin, double fMax);
void mesh_v(const Vector &x, const double t, Vector &res);
void move_mesh(ParFiniteElementSpace & pfes, const ParGridFunction * velocities);

// Fortran subroutine
extern "C" {
   void __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      double *in_rhol, double *in_ul, double *in_el, double *in_pl,
      double *in_rhor, double *in_ur, double *in_er, double *in_pr,
      double *in_tol, bool *no_iter,double *lambda_maxl_out,
      double *lambda_maxr_out, double *pstar, int *k);
}


int main(int argc, char *argv[]) {
   // Initialize MPI.
   MPI_Session mpi(argc, argv);
   const int myid = mpi.WorldRank();

   // Parse command line options
   dim = 2;
   const char *mesh_file = "default";
   int rs_levels = 0;
   int rp_levels = 0;
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
   // Can be visualized with glvis -np # -m ./results/mesh-test-init
   {
      ostringstream mesh_name;
      mesh_name << "./results/mesh-test-init." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   /* Print Face Information */
   pmesh->ExchangeFaceNbrData();
   cout << "num interior faces: " << pmesh->GetNFbyType(FaceType::Interior) << endl;
   cout << "num boundary faces: " << pmesh->GetNFbyType(FaceType::Boundary) << endl;
   cout << "num total faces: " << pmesh->GetNumFaces() << endl;
   cout << "NDofs: " << H1FESpace.GetNDofs() << endl;

   /*
   * Compute Normals
   */
   VectorFunctionCoefficient mesh_velocity(dim, mesh_v);
   mesh_velocity.SetTime(0);
   ParGridFunction *velocities = new ParGridFunction(&H1FESpace);
   velocities->ProjectCoefficient(mesh_velocity);
   build_C(H1FESpace, *velocities, 0.1); // Don't pass dereferenced pointers.
   double in_rhol,in_ul,in_el,in_pl,in_rhor,in_ur,in_er,in_pr,in_tol,lambda_maxl_out,lambda_maxr_out,pstar;
   bool no_iter;
   int k;
   __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(&in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,&no_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&k);
   
   /*
   * Move Mesh
   */
   move_mesh(H1FESpace, velocities);

   // Print moved mesh
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      ostringstream mesh_name;
      mesh_name << "./results/mesh-test-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }
}

/*
*
* Build normal matrix
*
*/
void build_C(ParFiniteElementSpace &pfes, ParGridFunction & velocities, double dt)
{
   /*
   iterate over faces
   TODO: 
      - Need to do specific things for interface vs boundary faces.
      - Incorporate direction of normal somehow.
   */
   Vector res(dim);
   Array<int> fids, oris;

   for (int i = 0; i < pfes.GetNE(); i++) // Cell iterator
   {
      cout << "Element # " << i << endl;
      pfes.ExchangeFaceNbrData();
      pfes.GetParMesh()->GetElementEdges(i, fids, oris); // 3DTODO: GetElementFaces()
      for (int j=0; j < fids.Size(); j++) // Face iterator
      {
         cout << "Face # " << fids[j] << endl;
         calc_outward_normal_int(i, fids[j], pfes, velocities, dt, res);
      } // End Face iterator
      cout << endl;
   } // End cell iterator
}

void calc_outward_normal_int(int cell, int face, ParFiniteElementSpace &pfes, ParGridFunction & velocities, double dt, Vector & res)
{
   // Fill nodes with mesh nodes for retrieval
   ParGridFunction nodes(&pfes);
   pfes.GetParMesh()->GetNodes(nodes);
   const int nNodes = nodes.Size() / dim;

   Array<int> face_dofs;
   Vector a(dim);
   Vector face_node(dim); // Temp
   Vector v(dim);
   res = 0.;

   // Get int der shape functions
   Vector shapes = Get_Int_Der_Ref_Shap_Functions();
   assert (dim == 2); // "This code only works for dim==2."

   pfes.GetFaceDofs(face, face_dofs);

   // Compute C_{cc'} = 1/2 \Sum_{i\in {1:3}} [(a_i^n + dtV_i^n) \int_0^1 \theta_i(\xi) d\xi]
   for (int d = 0; d < face_dofs.Size(); d++)
   {
      a = 0.;
      v = 0.;
      // cout << "d: " << d << ", dof: " << face_dofs[d] << endl;
      for (int _dim = 0; _dim < dim; _dim++)
      {
         int index = _dim * nNodes + face_dofs[d];
         a[_dim] = nodes(index);
         v[_dim] = velocities(index);
         // cout << "dim: " << _dim << ", index: " << index << ", val: " << a[_dim] << endl;
      }
      if (d == face_dofs.Size() - 1)
      {
         face_node = a;
      }

      // Now that vectors are retrieved, sum
      add(a, dt, v, a);
      a *= shapes[d];
      add(res, a, res);
   }

   orthogonal(res);
   res *= 0.5;

   /* Ensure orientation of normal */
   const auto FI = pfes.GetParMesh()->GetFaceInformation(face);

   if (FI.IsInterior()) // Interior face
   {
      // Orientation of the normal vector depends on the indices of
      // the cells that share that face.  The normal points towards
      // the cell with the lower global index.  This normal must be
      // flipped if the cell we are working with is the one of the
      // lower index.
      if (cell == min(FI.element[0].index, FI.element[1].index))
      {
         res *= -1.;
      }
   }
   else // Boundary face
   {
      // By default, boundary normals point toward cell.  This always
      // needs to be reversed to compute outward normal.
      assert(FI.IsBoundary());
      res *= -1.;
   }

   // Output normal
   {
      cout << "Face node coords: ";
      face_node.Print(cout);
      cout << "Normal vector: ";
      res.Print(cout);
   }
}

void orthogonal(Vector &v)
{
   if (dim == 2)
   {
      double x = v(0), y = v(1);
      v(0) = -1 * y;
      v(1) = x;
   }
}

Vector Get_Int_Der_Ref_Shap_Functions() 
{
   /* 
   Test stuff 
   Can use this to evaluate shape functions on the reference face.
   */
   // const FiniteElement * test_fe_p = H1FEC.FiniteElementForGeometry(Geometry::SEGMENT); 
   // cout << "FE derive type: " << test_fe_p->GetDerivType() << endl;
   // const IntegrationRule ir = test_fe_p->GetNodes();
   // cout << "num integration points: " << ir.GetNPoints() << endl;

   // for (int p = 0; p < ir.GetNPoints(); p++)
   // {
   //    const IntegrationPoint &ip = ir.IntPoint(p);
   //    cout << "ip x: " << ip.x << ", ip y: " << ip.y << ", ip weight: " << ip.weight << endl;
   //    DenseMatrix dm(3,2);
   //    test_fe_p->CalcDShape(ip, dm);
   //    dm.Print(cout);
   // }

   assert (dim == 2); // "This code only works for dim==2."

   Vector r(3);
   r[0] = -1.;
   r[1] = 1.;
   r[2] = 0.;

   return r;
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
   int mesh_opt = 0; // Change this param to get different movement

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
               return;
            }
            case 3: // random
            {
               const double r_max = .02;
               res[0] = fRand(0, r_max);
               res[1] = fRand(0, r_max);
               return;
            }
         } // mesh_opt switch
      } // 2D velocity end
   }
}

void move_mesh(ParFiniteElementSpace & pfes, const ParGridFunction * velocities)
{
   ParGridFunction *x_gf = new ParGridFunction(&pfes);

   // Initialize x_gf using the starting mesh coordinates.
   pfes.GetParMesh()->SetNodalGridFunction(x_gf);

   // Logic to change x_gf based on velocities in v_gf
   // V = v0 + 0.5 * dt * dv_dt;
   add(*x_gf, 1, *velocities, *x_gf);

   // Update mesh gridfunction
   pfes.GetParMesh()->NewNodes(*x_gf, false);
}
