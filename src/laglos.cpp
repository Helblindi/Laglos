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
#include "initial_vals.hpp"
#include "compile_time_vals.h"

using namespace mfem;
using namespace hydrodynamics;
using namespace std;

// Forward declarations
void orthogonal(Vector &v);
void build_C(ParFiniteElementSpace &pfes, ParGridFunction & velocities, double dt);
void calc_outward_normal_int(int cell, int face, ParFiniteElementSpace &pfes, ParGridFunction & velocities, double dt, Vector & res);
Vector Get_Int_Der_Ref_Shap_Functions();

// double compute_viscosity();
void move_mesh(ParFiniteElementSpace & pfes, const ParGridFunction * velocities);

double fRand(double fMin, double fMax);
void mesh_v(const Vector &x, const double t, Vector &res); // Initial mesh velocity


int main(int argc, char *argv[]) {
   // Initialize MPI.
   MPI_Session mpi(argc, argv);
   const int myid = mpi.WorldRank();

   const int dim = CompileTimeVals::dim;
   const int problem = CompileTimeVals::problem;
   // const int dim = 2;
   // const int problem = 1;

   // Parse command line options
   const char *mesh_file = "default";
   int rs_levels = 0;
   int rp_levels = 0;
   int order_mv = 2;  // Order of mesh movement approximation space
   int order_u = 0;
   double t_final = 0.6;
   int max_tsteps = -1;
   bool visualization = false;
   int vis_steps = 5;
   int precision = 12;

   OptionsParser args(argc, argv);
   // args.AddOption(&problem, "-p", "--problem", "Problem setup to use.");
   // args.AddOption(&dim, "-dim", "--dimension", "Dimension of the problem.");
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&rs_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&rp_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&max_tsteps, "-ms", "--max-steps",
                  "Maximum number of steps (negative means no restriction).");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");

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
   if (dim != mesh->Dimension())
   {
      cout << "Mesh dimension does not match compile time vals dimension\n"
           << "Aborting program.\n";
      return -1;
   }
   cout << "Meshing done\n";
   // dim = mesh->Dimension(); // Set the correct dimension if a mesh other than
   //                          // 'default' was provided.

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

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);

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
   FunctionCoefficient sv_coeff(InitialValues<problem, dim>::sv0);
   sv_gf.ProjectCoefficient(sv_coeff);
   sv_gf.SyncAliasMemory(S);

   VectorFunctionCoefficient v_coeff(dim, InitialValues<problem, dim>::v0);
   v_gf.ProjectCoefficient(v_coeff);
   // for (int i = 0; i < ess_vdofs.Size(); i++)
   // {
   //    v_gf(ess_vdofs[i]) = 0.0;
   // }
   // Sync the data location of v_gf with its base, S
   v_gf.SyncAliasMemory(S);

   FunctionCoefficient ste_coeff(InitialValues<problem, dim>::sie0);
   ste_gf.ProjectCoefficient(ste_coeff);
   ste_gf.SyncAliasMemory(S);

   // PLF to build mass vector
   FunctionCoefficient rho0_coeff(InitialValues<problem, dim>::rho0); // To build mass vector
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(rho0_coeff));
   m->Assemble();

   /* Create Lagrangian Low Order Solver Object */
   LagrangianLOOperator<dim> hydro(H1FESpace, L2FESpace, L2VFESpace, m);

   socketstream vis_sv, vis_v, vis_ste;
   char vishost[] = "localhost";
   int  visport   = 19916;

   if (visualization)
   {
      // Make sure all MPI ranks have sent their 'v' solution before initiating
      // another set of GLVis connections (one from each rank):
      MPI_Barrier(pmesh->GetComm());
      vis_sv.precision(8);
      vis_v.precision(8);
      vis_ste.precision(8);
      int Wx = 0, Wy = 0; // window position
      const int Ww = 350, Wh = 350; // window size
      int offx = Ww+10; // window offsets
      if (problem != 0 && problem != 4) // TODO: Remove??
      {
         VisualizeField(vis_sv, vishost, visport, sv_gf,
                                       "Specific Volume", Wx, Wy, Ww, Wh);
      }
      Wx += offx;
      VisualizeField(vis_v, vishost, visport, v_gf,
                                    "Velocity", Wx, Wy, Ww, Wh);
      Wx += offx;
      VisualizeField(vis_ste, vishost, visport, ste_gf,
                                    "Specific Total Energy", Wx, Wy, Ww, Wh);
   }

   // Perform the time-integration by looping over time iterations
   // ti with a time step dt.  The main function call here is the
   // LagrangianLOOperator.MakeTimeStep() funcall.
   double t = 0.0, dt = 0.0001, t_old;
   bool last_step = false;
   int steps = 0;
   BlockVector S_old(S);

   for (int ti = 1; !last_step; ti++)
   {
      hydro.CalculateTimestep(S);
      dt = hydro.GetTimestep();

      if (t + dt >= t_final)
      {
         dt = t_final - t;
         last_step = true;
      }
      if (steps == max_tsteps) { last_step = true; }

      S_old = S;
      t_old = t;

      hydro.MakeTimeStep(S, t, dt); // testing
      steps++;
      if (S_old == S) { cout << "\t ########## We have a problem!\n"; }
      // TODO: Adaptive time step control?

      // Ensure the sub-vectors x_gf, v_gf, and ste_gf know the location of the
      // data in S. This operation simply updates the Memory validity flags of
      // the sub-vectors to match those of S.
      x_gf.SyncAliasMemory(S);
      mv_gf.SyncAliasMemory(S);
      sv_gf.SyncAliasMemory(S);
      v_gf.SyncAliasMemory(S);
      ste_gf.SyncAliasMemory(S);

      // Make sure that the mesh corresponds to the new solution state. This is
      // needed, because some time integrators use different S-type vectors
      // and the oper object might have redirected the mesh positions to those.
      pmesh->NewNodes(x_gf, false);

      if (last_step || (ti % vis_steps) == 0)
      {
         double lnorm = ste_gf * ste_gf, norm;
         MPI_Allreduce(&lnorm, &norm, 1, MPI_DOUBLE, MPI_SUM, pmesh->GetComm());
         // if (mem_usage)
         // {
         //    mem = GetMaxRssMB();
         //    MPI_Reduce(&mem, &mmax, 1, MPI_LONG, MPI_MAX, 0, pmesh->GetComm());
         //    MPI_Reduce(&mem, &msum, 1, MPI_LONG, MPI_SUM, 0, pmesh->GetComm());
         // }
         // const double internal_energy = hydro.InternalEnergy(e_gf);
         // const double kinetic_energy = hydro.KineticEnergy(v_gf);
         if (mpi.Root())
         {
            const double sqrt_norm = sqrt(norm);

            cout << std::fixed;
            cout << "step " << std::setw(5) << ti
                 << ",\tt = " << std::setw(5) << std::setprecision(4) << t
                 << ",\tdt = " << std::setw(5) << std::setprecision(6) << dt
                 << ",\t|e| = " << std::setprecision(10) << std::scientific
                 << sqrt_norm;
            //  << ",\t|IE| = " << std::setprecision(10) << std::scientific
            //  << internal_energy
            //   << ",\t|KE| = " << std::setprecision(10) << std::scientific
            //  << kinetic_energy
            //   << ",\t|E| = " << std::setprecision(10) << std::scientific
            //  << kinetic_energy+internal_energy;
            // cout << std::fixed;
            // if (mem_usage)
            // {
            //    cout << ", mem: " << mmax << "/" << msum << " MB";
            // }
            cout << endl;
         }

         if (visualization)
         {
            int Wx = 0, Wy = 0; // window position
            int Ww = 350, Wh = 350; // window size
            int offx = Ww+10; // window offsets
            if (problem != 0 && problem != 4)
            {
               VisualizeField(vis_sv, vishost, visport, sv_gf,
                              "Specific Volume", Wx, Wy, Ww, Wh);
            }
            Wx += offx;
            VisualizeField(vis_v, vishost, visport,
                           v_gf, "Velocity", Wx, Wy, Ww, Wh);
            Wx += offx;
            VisualizeField(vis_ste, vishost, visport, ste_gf,
                           "Specific Total Energy",
                           Wx, Wy, Ww,Wh);
            Wx += offx;
         }
      }
   } // End time step iteration

   /*
   * Move Mesh
   */
   // move_mesh(H1FESpace, &mv_gf);

   // Print moved mesh
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      ostringstream mesh_name;
      mesh_name << "./results/mesh-test-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   delete pmesh;

   return 0;
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
   res.SetSize(x.Size());
   res = 0.;
}
// void mesh_v(const Vector &x, const double t, Vector &res)
// {
//    int mesh_opt = 0; // Change this param to get different movement

//    switch (dim)
//    {
//       case 1:
//       {
//          res[0] = .1;
//          return;
//       }
//       case 3:
//       case 2:
//       {
//          switch (mesh_opt)
//          {
//             case 0: // constant velocity
//             {
//                for (int i = 0; i < dim; i++)
//                {
//                   res[i] = 0.2;
//                }
//                return;
//             }
//             case 1: // Move non corners
//             {
//                // Set v[1]
//                if (x[0] != 0 && x[0] != 1)
//                {
//                   res[1] = 0.1;
//                }
//                else
//                {
//                   res[1] = 0.;
//                }

//                // Set v[0]
//                if (x[1] != 0 && x[1] != 1)
//                {
//                   res[0] = .2;
//                }
//                else
//                { 
//                   res[0] = 0.;
//                }
//                return;
//             }
//             case 2: // rotation about origin with angle theta
//             {
//                const double theta = M_PI/4;
//                res[0] = x[0]*(cos(theta) - 1) - x[1]*sin(theta);
//                res[1] = x[0]*sin(theta) + x[1]*(cos(theta) - 1);
//                return;
//             }
//             case 3: // random
//             {
//                const double r_max = .02;
//                res[0] = fRand(0, r_max);
//                res[1] = fRand(0, r_max);
//                return;
//             }
//          } // mesh_opt switch
//       } // 2D velocity end
//    }
// }