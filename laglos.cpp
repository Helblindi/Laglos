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
#include "laglos_solver.hpp"

using namespace mfem;
using namespace std;

// Choice for the problem setup.
static int problem, dim;

// Forward declarations
void orthogonal(Vector &v);
void build_C(ParFiniteElementSpace &pfes, ParGridFunction & velocities, double dt);
void calc_outward_normal_int(int cell, int face, ParFiniteElementSpace &pfes, ParGridFunction & velocities, double dt, Vector & res);
Vector Get_Int_Der_Ref_Shap_Functions();

// double compute_viscosity();
void move_mesh(ParFiniteElementSpace & pfes, const ParGridFunction * velocities);

double fRand(double fMin, double fMax);
void mesh_v(const Vector &x, const double t, Vector &res); // Initial mesh velocity
double rho0(const Vector &x); // from laghos
double sv0(const Vector &x); // inverse rho0
double ste0(const Vector &x);
void v0(const Vector &x, Vector &v);
static double rad(double x, double y);
double gamma_func(const Vector &x);


int main(int argc, char *argv[]) {
   // Initialize MPI.
   MPI_Session mpi(argc, argv);
   const int myid = mpi.WorldRank();

   // Parse command line options
   problem = 1;
   dim = 2;
   const char *mesh_file = "default";
   int rs_levels = 0;
   int rp_levels = 0;
   int order_mv = 2;  // Order of mesh movement approximation space
   int order_u = 0;
   int precision = 12;

   OptionsParser args(argc, argv);
   args.AddOption(&problem, "-p", "--problem", "Problem setup to use.");
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
      }
      if (dim == 2)
      {
         mesh = new Mesh(Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL,
                                               true));
      }
      if (dim == 3)
      {
         mesh = new Mesh(Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON,
                                               true));
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
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, pmesh->Dimension());

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

   // Print # DoFs
   const HYPRE_Int glob_size_l2 = L2FESpace.GlobalTrueVSize();
   const HYPRE_Int glob_size_h1 = H1FESpace.GlobalTrueVSize();
   if (mpi.Root())
   {
      cout << "Number of kinematic (position, mesh velocity) dofs: "
           << glob_size_h1 << endl
           << "Corresponding NDofs: "
           << H1FESpace.GetNDofs() << endl
           << "Corresponding VDim: " 
           << H1FESpace.GetVDim() << endl
           << "Corresponding Vsize: "
           << H1FESpace.GetVSize() << endl;

      cout << "Number of specific internal energy dofs: "
           << glob_size_l2 << endl
           << "Corresponding NDofs: "
           << L2FESpace.GetNDofs() << endl
           << "Corresponding VDim: " 
           << L2FESpace.GetVDim() << endl
           << "Corresponding Vsize: "
           << L2FESpace.GetVSize() << endl;
      
      cout << "Number of velocity dofs: "
           << L2VFESpace.GlobalTrueVSize() << endl
           << "Corresponding NDofs: "
           << L2VFESpace.GetNDofs() << endl
           << "Corresponding VDim: " 
           << L2VFESpace.GetVDim() << endl
           << "Corresponding Vsize: "
           << L2VFESpace.GetVSize() << endl;
   }

   /* Print Face Information */
   pmesh->ExchangeFaceNbrData();
   cout << "num interior faces: " << pmesh->GetNFbyType(FaceType::Interior) << endl;
   cout << "num boundary faces: " << pmesh->GetNFbyType(FaceType::Boundary) << endl;
   cout << "num total faces: " << pmesh->GetNumFaces() << endl;
   

   /* The monolithic BlockVector stores unknown fields as:
   *   - 0 -> position
   *   - 1 -> mesh velocity
   *   - 2 -> specific volume
   *   - 3 -> velocity (L2V)
   *   - 4 -> speific total energy
   */
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

   // Define GridFunction objects for the position, mesh velocity and specific
   // volume, velocity, and specific internal energy. At each step, each of 
   // these values will be updated.
   ParGridFunction x_gf, mv_gf, sv_gf, v_gf, ste_gf;
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   mv_gf.MakeRef(&H1FESpace, S, offset[1]);
   sv_gf.MakeRef(&L2FESpace, S, offset[2]);
   v_gf.MakeRef(&L2VFESpace, S, offset[3]);
   ste_gf.MakeRef(&L2FESpace, S, offset[4]);

   // Initialize x_gf using starting mesh positions
   pmesh->SetNodalGridFunction(&x_gf);
   // Sync the data location of x_gf with its base, S
   x_gf.SyncAliasMemory(S);

   // Initialize mesh velocity
   VectorFunctionCoefficient mesh_velocity(dim, mesh_v);
   mesh_velocity.SetTime(0);
   mv_gf.ProjectCoefficient(mesh_velocity);
   // Sync the data location of mv_gf with its base, S
   mv_gf.SyncAliasMemory(S);

   // Initialize specific volume, velocity, and specific total energy
   FunctionCoefficient sv_coeff(sv0);
   sv_gf.ProjectCoefficient(sv_coeff);
   sv_gf.SyncAliasMemory(S);

   VectorFunctionCoefficient v_coeff(pmesh->Dimension(), v0);
   v_gf.ProjectCoefficient(v_coeff);
   // for (int i = 0; i < ess_vdofs.Size(); i++)
   // {
   //    v_gf(ess_vdofs[i]) = 0.0;
   // }
   // Sync the data location of v_gf with its base, S
   v_gf.SyncAliasMemory(S);

   FunctionCoefficient ste_coeff(ste0);
   ste_gf.ProjectCoefficient(ste_coeff);
   ste_gf.SyncAliasMemory(S);

   // /*
   // * Compute Normals
   // */
   build_C(H1FESpace, mv_gf, 0.1); // Don't pass dereferenced pointers.
   // double _visc = compute_viscosity();
   LagrangianLOOperator hydro(S.Size(), H1FESpace, L2FESpace, L2VFESpace);
   Vector U;
   // Testing functions
   // hydro.GetCellStateVector(S, 0, U);
   // U.Print(cout);
   // U = 1;
   // hydro.SetCellStateVector(S, 0, U);
   // U.Print(cout);


   /*
   * Move Mesh
   */
   move_mesh(H1FESpace, &mv_gf);

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
* Viscosity computation
*
*/
// double compute_viscosity()
// {
//    double in_rhol,in_ul,in_el,in_pl,in_rhor,in_ur,in_er,in_pr,in_tol,lambda_maxl_out,lambda_maxr_out,pstar;
//    bool no_iter;
//    int k;
//    __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
//       &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,
//       &in_tol,&no_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&k);

//    return 1.0;
// }

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

/*
*
* Initial Conditions
*
*/
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

double sv0(const Vector &x)
{
   double val = rho0(x);
   assert(val != 0.);
   return 1./val;
}

double rho0(const Vector &x)
{
   /* Copied from Laghos */
   switch (problem)
   {
      case 0: return 1.0;
      case 1: return 1.0;
      case 2: return (x(0) < 0.5) ? 1.0 : 0.1;
      case 3: return (dim == 2) ? (x(0) > 1.0 && x(1) > 1.5) ? 0.125 : 1.0
                        : x(0) > 1.0 && ((x(1) < 1.5 && x(2) < 1.5) ||
                                         (x(1) > 1.5 && x(2) > 1.5)) ? 0.125 : 1.0;
      case 4: return 1.0;
      case 5:
      {
         if (x(0) >= 0.5 && x(1) >= 0.5) { return 0.5313; }
         if (x(0) <  0.5 && x(1) <  0.5) { return 0.8; }
         return 1.0;
      }
      case 6:
      {
         if (x(0) <  0.5 && x(1) >= 0.5) { return 2.0; }
         if (x(0) >= 0.5 && x(1) <  0.5) { return 3.0; }
         return 1.0;
      }
      case 7: return x(1) >= 0.0 ? 2.0 : 1.0;
      default: MFEM_ABORT("Bad number given for problem id!"); return 0.0;
   }
}

double gamma_func(const Vector &x)
{
   switch (problem)
   {
      case 0: return 5.0 / 3.0;
      case 1: return 1.4;
      case 2: return 1.4;
      case 3: return (x(0) > 1.0 && x(1) <= 1.5) ? 1.4 : 1.5;
      case 4: return 5.0 / 3.0;
      case 5: return 1.4;
      case 6: return 1.4;
      case 7: return 5.0 / 3.0;
      default: MFEM_ABORT("Bad number given for problem id!"); return 0.0;
   }
}

static double rad(double x, double y) { return sqrt(x*x + y*y); }

double ste0(const Vector &x)
{
   switch (problem)
   {
      case 0:
      {
         const double denom = 2.0 / 3.0;  // (5/3 - 1) * density.
         double val;
         if (x.Size() == 2)
         {
            val = 1.0 + (cos(2*M_PI*x(0)) + cos(2*M_PI*x(1))) / 4.0;
         }
         else
         {
            val = 100.0 + ((cos(2*M_PI*x(2)) + 2) *
                           (cos(2*M_PI*x(0)) + cos(2*M_PI*x(1))) - 2) / 16.0;
         }
         return val/denom;
      }
      case 1: return 0.0; // This case in initialized in main().
      case 2: return (x(0) < 0.5) ? 1.0 / rho0(x) / (gamma_func(x) - 1.0)
                        : 0.1 / rho0(x) / (gamma_func(x) - 1.0);
      case 3: return (x(0) > 1.0) ? 0.1 / rho0(x) / (gamma_func(x) - 1.0)
                        : 1.0 / rho0(x) / (gamma_func(x) - 1.0);
      case 4:
      {
         const double r = rad(x(0), x(1)), rsq = x(0) * x(0) + x(1) * x(1);
         const double gamma = 5.0 / 3.0;
         if (r < 0.2)
         {
            return (5.0 + 25.0 / 2.0 * rsq) / (gamma - 1.0);
         }
         else if (r < 0.4)
         {
            const double t1 = 9.0 - 4.0 * log(0.2) + 25.0 / 2.0 * rsq;
            const double t2 = 20.0 * r - 4.0 * log(r);
            return (t1 - t2) / (gamma - 1.0);
         }
         else { return (3.0 + 4.0 * log(2.0)) / (gamma - 1.0); }
      }
      case 5:
      {
         const double irg = 1.0 / rho0(x) / (gamma_func(x) - 1.0);
         if (x(0) >= 0.5 && x(1) >= 0.5) { return 0.4 * irg; }
         if (x(0) <  0.5 && x(1) >= 0.5) { return 1.0 * irg; }
         if (x(0) <  0.5 && x(1) <  0.5) { return 1.0 * irg; }
         if (x(0) >= 0.5 && x(1) <  0.5) { return 1.0 * irg; }
         MFEM_ABORT("Error in problem 5!");
         return 0.0;
      }
      case 6:
      {
         const double irg = 1.0 / rho0(x) / (gamma_func(x) - 1.0);
         if (x(0) >= 0.5 && x(1) >= 0.5) { return 1.0 * irg; }
         if (x(0) <  0.5 && x(1) >= 0.5) { return 1.0 * irg; }
         if (x(0) <  0.5 && x(1) <  0.5) { return 1.0 * irg; }
         if (x(0) >= 0.5 && x(1) <  0.5) { return 1.0 * irg; }
         MFEM_ABORT("Error in problem 5!");
         return 0.0;
      }
      case 7:
      {
         const double rho = rho0(x), gamma = gamma_func(x);
         return (6.0 - rho * x(1)) / (gamma - 1.0) / rho;
      }
      default: MFEM_ABORT("Bad number given for problem id!"); return 0.0;
   }
}

void v0(const Vector &x, Vector &v)
{
   const double atn = pow((x(0)*(1.0-x(0))*4*x(1)*(1.0-x(1))*4.0),0.4);
   switch (problem)
   {
      case 0:
         v(0) =  sin(M_PI*x(0)) * cos(M_PI*x(1));
         v(1) = -cos(M_PI*x(0)) * sin(M_PI*x(1));
         if (x.Size() == 3)
         {
            v(0) *= cos(M_PI*x(2));
            v(1) *= cos(M_PI*x(2));
            v(2) = 0.0;
         }
         break;
      case 1: v = 0.0; break;
      case 2: v = 0.0; break;
      case 3: v = 0.0; break;
      case 4:
      {
         v = 0.0;
         const double r = rad(x(0), x(1));
         if (r < 0.2)
         {
            v(0) =  5.0 * x(1);
            v(1) = -5.0 * x(0);
         }
         else if (r < 0.4)
         {
            v(0) =  2.0 * x(1) / r - 5.0 * x(1);
            v(1) = -2.0 * x(0) / r + 5.0 * x(0);
         }
         else { }
         break;
      }
      case 5:
      {
         v = 0.0;
         if (x(0) >= 0.5 && x(1) >= 0.5) { v(0)=0.0*atn, v(1)=0.0*atn; return;}
         if (x(0) <  0.5 && x(1) >= 0.5) { v(0)=0.7276*atn, v(1)=0.0*atn; return;}
         if (x(0) <  0.5 && x(1) <  0.5) { v(0)=0.0*atn, v(1)=0.0*atn; return;}
         if (x(0) >= 0.5 && x(1) <  0.5) { v(0)=0.0*atn, v(1)=0.7276*atn; return; }
         MFEM_ABORT("Error in problem 5!");
         return;
      }
      case 6:
      {
         v = 0.0;
         if (x(0) >= 0.5 && x(1) >= 0.5) { v(0)=+0.75*atn, v(1)=-0.5*atn; return;}
         if (x(0) <  0.5 && x(1) >= 0.5) { v(0)=+0.75*atn, v(1)=+0.5*atn; return;}
         if (x(0) <  0.5 && x(1) <  0.5) { v(0)=-0.75*atn, v(1)=+0.5*atn; return;}
         if (x(0) >= 0.5 && x(1) <  0.5) { v(0)=-0.75*atn, v(1)=-0.5*atn; return;}
         MFEM_ABORT("Error in problem 6!");
         return;
      }
      case 7:
      {
         v = 0.0;
         v(1) = 0.02 * exp(-2*M_PI*x(1)*x(1)) * cos(2*M_PI*x(0));
         break;
      }
      default: MFEM_ABORT("Bad number given for problem id!");
   }
}