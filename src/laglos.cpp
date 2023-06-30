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

/*
* Example run time parameters:  [Remember to change the problem in compile-time-vals.h]
*
* ./Laglos -m data/ref-square-c0.mesh -tf 2 -cfl 0.5 -ot -visc -mm -vis -rs 3 [problem = 2, dim = 2] Isentropic vortex, stationary center
* ./Laglos -m ../data/ref-square-c0.mesh -tf 2 -cfl 0.5 -ot -visc -mm -vis -rs 5 [problem =3, dim=2]
* ./Laglos -m data/ref-square.mesh -tf 1. -cfl 0.5 -ot -visc -mm -vis -rs 3 [problem = 5, dim = 2]
* ./Laglos -m ../data/ref-tube.mesh -tf 0.225 -cfl 0.5 -ot -visc -mm -vis -rs 5 [problem = 6, dim = 2, shocktube = 1] // Sod
* ./Laglos -m ../data/ref-square-tube.mesh -tf 0.67 -cfl 0.2 -ot -visc -mm -vis -rs 5 [problem = 6, dim = 2, shocktube = 3]
* ./Laglos -m ../data/rectangle_saltzman.mesh -rs 3 -visc -mm -vis -tf 0.6 -ot -cfl 0.01 [problem = 7, dim = 2] // Saltzman problem
* ./Laglos -m ../data/ref-square.mesh -rs 3 -visc -mm -tf 0.5 -ot -vis [problem = 8, dim = 2]
* 
*/


#include "compile_time_vals.h"
#include "mfem.hpp"
#include "laglos_solver.hpp"
#include <iostream>
#include <fstream>
#include <chrono>

#define _USE_MATH_DEFINES
#include <cmath>
#include "matplotlibcpp.h"

#ifdef DMALLOC  // Memory allocation checks
#include "dmalloc.h"
#endif

using namespace mfem;
using namespace hydrodynamics;
using namespace std;
namespace plt = matplotlibcpp;

// // Forward declarations
// double free_steam_conditions(const Vector &x);
// void free_steam_conditions_v(const Vector &x, Vector &v);
// void orthogonal(Vector &v);
// void build_C(ParFiniteElementSpace &pfes, ParGridFunction & velocities, double dt);
// void calc_outward_normal_int(int cell, int face, ParFiniteElementSpace &pfes, ParGridFunction & velocities, double dt, Vector & res);
// Vector Get_Int_Der_Ref_Shap_Functions();

// // double compute_viscosity();
// void move_mesh(ParFiniteElementSpace & pfes, const ParGridFunction * velocities);

// double fRand(double fMin, double fMax);
// void mesh_v(const Vector &x, const double t, Vector &res); // Initial mesh velocity


int main(int argc, char *argv[]) {
   // Initialize MPI.
   Mpi::Init();
   const int num_procs = Mpi::WorldSize();
   const int myid = Mpi::WorldRank();
   // Hypre::Init();

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
   double dt = 0.001;
   bool visualization = false;
   int vis_steps = 5;
   int precision = 12;
   bool use_viscosity = true;
   bool mm = false;
   bool optimize_timestep = false;
   bool convergence_testing = false;
   double CFL = 0.5;

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
   args.AddOption(&use_viscosity, "-visc", "--use-viscosity", "-no-visc",
                  "--no-viscosity",
                  "Enable or disable the use of artificial viscosity.");
   args.AddOption(&mm, "-mm", "--move-mesh", "-no-mm", "--no-move-mesh",
                  "Enable or disable mesh movement.");
   args.AddOption(&optimize_timestep, "-ot", "--optimize-timestep", "-no-ot",
                  "--no-optimize-timestep",
                  "Enable or disable timestep optimization using CFL.");
   args.AddOption(&CFL, "-cfl", "--time-step-restriction",
                  "CFL value to use.");
   args.AddOption(&convergence_testing, "-ct", "--set-dt-to-h", "-no-ct", 
                  "--no-set-dt-to-h", 
                  "Enable or disable convergence testing timestep.");
   args.AddOption(&dt, "-dt", "--timestep", "Timestep to use.");

   args.Parse();
   if (!args.Good())
   {
      if (Mpi::Root()) { args.PrintUsage(cout); }
      return 1;
   }
   if (Mpi::Root()) { args.PrintOptions(cout); }

   // Check that convergence testing and optimizing the timestep are not both set to true
   if (optimize_timestep && convergence_testing)
   {
      cout << "Cannot both optimize the timestep and set the timestep for convergence testing.\n";
      return -1;
   }

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

   cout << "done with serial refinement\n";
   // Define the parallel mesh by a partitioning of the serial mesh. Refine
   // this mesh further in parallel to increase the resolution. Once the
   // parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);            
   delete mesh;
   for (int lev = 0; lev < rp_levels; lev++)
   {
      pmesh->UniformRefinement();
   }
   cout << "done with parallel refinement\n";

   cout << "dim: " << pmesh->Dimension() << ", spacedim: " << pmesh->SpaceDimension() << endl;

   int NE = pmesh->GetNE(), ne_min, ne_max;
   MPI_Reduce(&NE, &ne_min, 1, MPI_INT, MPI_MIN, 0, pmesh->GetComm());
   MPI_Reduce(&NE, &ne_max, 1, MPI_INT, MPI_MAX, 0, pmesh->GetComm());
   double hmin, hmax, kmin, kmax;
   pmesh->GetCharacteristics(hmin, hmax, kmin, kmax);
   if (myid == 0)
   { cout << "Zones min/max: " << ne_min << " " << ne_max << endl; }

   // Define the parallel finite element spaces. We use:
   // - H1 (Q2, continuous) for mesh movement.
   // - L2 (Q0, discontinuous) for state variables
   H1_FECollection H1FEC(order_mv, dim);
   // H1Ser_FECollection H1FEC(order_mv, dim);
   L2_FECollection L2FEC(order_u, dim, BasisType::Positive);
   CrouzeixRaviartFECollection CRFEC;

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);
   ParFiniteElementSpace CRFESpace(pmesh, &CRFEC, dim);

   HYPRE_BigInt global_vSize = L2FESpace.GlobalTrueVSize();
   cout << "global_vSize: " << global_vSize << endl;
   cout << "N dofs: " << H1FESpace.GetNDofs() << endl;
   cout << "Num Scalar Vertex DoFs: " << H1FESpace.GetNVDofs() << endl;
   cout << "Num Scalar edge-interior DoFs: " << H1FESpace.GetNEDofs() << endl;
   cout << "Num Scalar face-interior DoFs: " << H1FESpace.GetNFDofs() << endl;
   cout << "Num faces: " << H1FESpace.GetNF() << endl;
   cout << "dim: " << dim << endl;

   // Print # DoFs
   const HYPRE_Int glob_size_l2 = L2FESpace.GlobalTrueVSize();
   const HYPRE_Int glob_size_h1 = H1FESpace.GlobalTrueVSize();
   if (Mpi::Root())
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

   cout << "Block Vector constructed\n";

   // Define GridFunction objects for the position, mesh velocity and specific
   // volume, velocity, and specific internal energy. At each step, each of 
   // these values will be updated.
   ParGridFunction x_gf, mv_gf, sv_gf, v_gf, ste_gf;
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   mv_gf.MakeRef(&H1FESpace, S, offset[1]);
   sv_gf.MakeRef(&L2FESpace, S, offset[2]);
   v_gf.MakeRef(&L2VFESpace, S, offset[3]);
   ste_gf.MakeRef(&L2FESpace, S, offset[4]);

   cout << "grid functions associated\n";

   // Initialize x_gf using starting mesh positions
   pmesh->SetNodalGridFunction(&x_gf);
   // Sync the data location of x_gf with its base, S
   x_gf.SyncAliasMemory(S);

   // cout << "Printing initial x_gf positions\n";
   // x_gf.Print(cout);

   /* Distort Mesh for Saltzman Problem */
   if (CompileTimeVals::distort_mesh)
   {
      cout << "distorting mesh\n";
      switch (problem)
      {
         case 7: // Saltzman problem
         {
            Array<double> coords(dim);
            for (int vertex = 0; vertex < H1FESpace.GetNDofs(); vertex++)
            {
               for (int i = 0; i < dim; i++)
               {
                  int index = vertex + i * H1FESpace.GetNDofs();
                  coords[i] = x_gf[index];
               }
               // Only the x-coord is distorted, according to Guermond, Popov, Saavedra
               double x_new = coords[0] + (0.1 - coords[1]) * sin(M_PI * coords[0]); // 

               x_gf[vertex] = x_new;
            }
         }
         // case 8: // Linear node movement
         // {
         //    Array<double> coords(dim), coords_new(dim);
         //    cout << "pmesh call: " << pmesh->GetNV() << endl;
         //    cout << "H1 call: " << H1FESpace.GetNDofs() << endl;
         //    cout << "pmesh num faces: " << pmesh->GetNumFaces() << endl;
         
         //    // Distort all corner nodes
         //    // pmesh->GetNV will give # corner nodes 
         //    for (int vertex = 0; vertex < pmesh->GetNV(); vertex++)
         //    {
         //       cout << "vertex: " << vertex << endl;
         //       for (int i = 0; i < dim; i++)
         //       {
         //          int index = vertex + i * H1FESpace.GetNDofs();
         //          cout << "index: " << index << endl;

         //          coords_new[i] = coords[i]; //*pow(coords[i], 2);
         //          x_gf[index] = coords_new[i];
         //       }
         //    }

         //    // Adjust all face nodes to be averages of their adjacent corners
         //    mfem::Mesh::FaceInformation FI;
         //    for (int face = 0; face < pmesh->GetNumFaces(); face++)
         //    {
         //       cout << "Iterating on face: " << face << endl;
         //       FI = pmesh->GetFaceInformation(face);
         //       Array<int> face_dof_row;
         //       H1FESpace.GetFaceDofs(face, face_dof_row);

         //       cout << "Face dof row for face: " << face << ":\n";
         //       face_dof_row.Print(cout);

         //       int face_dof = face_dof_row[2];
         //       int index0 = face_dof_row[0];
         //       int index1 = face_dof_row[1];
         //       for (int i = 0; i < dim; i++)
         //       {
         //          int face_dof_index = face_dof + i * H1FESpace.GetNDofs();
         //          int node0_index = index0 + i * H1FESpace.GetNDofs();
         //          int node1_index = index1 + i * H1FESpace.GetNDofs();
         //          x_gf[face_dof_index] = 0.5 * (x_gf[node0_index] + x_gf[node1_index]);
         //       }
               
         //    }
         // }
         default: // do nothing for all other problems
         {
         }
      }
      cout << "Mesh distorted\n";
   }
   
   // cout << "Printing initial x_gf positions\n";
   // x_gf.Print(cout);

   // Print initialized mesh
   // Can be visualized with glvis -np # -m ./results/mesh-test-init
   {
      ostringstream mesh_name;
      mesh_name << "./results/mesh-test-init." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   // Initialize mesh velocity
   ConstantCoefficient zero(0.0);
   mv_gf.ProjectCoefficient(zero);
   // Sync the data location of mv_gf with its base, S
   mv_gf.SyncAliasMemory(S);

   // Initialize specific volume, velocity, and specific total energy
   FunctionCoefficient sv_coeff(InitialValues<problem, dim>::sv0);
   sv_coeff.SetTime(0.);
   sv_gf.ProjectCoefficient(sv_coeff);
   sv_gf.SyncAliasMemory(S);

   VectorFunctionCoefficient v_coeff(dim, InitialValues<problem, dim>::v0);
   v_coeff.SetTime(0.);
   v_gf.ProjectCoefficient(v_coeff);
   // Sync the data location of v_gf with its base, S
   v_gf.SyncAliasMemory(S);

   FunctionCoefficient ste_coeff(InitialValues<problem, dim>::ste0);
   ste_coeff.SetTime(0.);
   ste_gf.ProjectCoefficient(ste_coeff);
   ste_gf.SyncAliasMemory(S);

   // PLF to build mass vector
   FunctionCoefficient rho_coeff(InitialValues<problem, dim>::rho0); 
   rho_coeff.SetTime(0.);
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff));
   m->Assemble();

   /* Create Lagrangian Low Order Solver Object */
   LagrangianLOOperator<dim, problem> hydro(H1FESpace, L2FESpace, L2VFESpace, CRFESpace, m, use_viscosity, mm, CFL);
   cout << "Solver created.\n";

   /* Set up visualiztion object */
   socketstream vis_rho, vis_v, vis_ste, vis_rho_ex, vis_v_ex, vis_ste_ex, vis_rho_err, vis_v_err, vis_ste_err;
   char vishost[] = "localhost";
   int  visport   = 19916;

   ParGridFunction rho_gf(&L2FESpace);

   if (visualization)
   {
      // Compute Density
      for (int i = 0; i < sv_gf.Size(); i++)
      {
         rho_gf[i] = 1./sv_gf[i];
      } 

      // Make sure all MPI ranks have sent their 'v' solution before initiating
      // another set of GLVis connections (one from each rank):
      MPI_Barrier(pmesh->GetComm());

      vis_rho.precision(8);
      vis_v.precision(8);
      vis_ste.precision(8);

      vis_rho_ex.precision(8);
      vis_v_ex.precision(8);
      vis_ste_ex.precision(8);

      vis_rho_err.precision(8);
      vis_v_err.precision(8);
      vis_ste_err.precision(8);

      int Wx = 0, Wy = 0; // window position
      const int Ww = 350, Wh = 350; // window size
      int offx = Ww+10, offy = Wh+45;; // window offsets

      VisualizeField(vis_rho, vishost, visport, rho_gf,
                     "Density", Wx, Wy, Ww, Wh);
      Wx += offx;
      VisualizeField(vis_v, vishost, visport, v_gf,
                     "Velocity", Wx, Wy, Ww, Wh);
      Wx += offx;
      VisualizeField(vis_ste, vishost, visport, ste_gf,
                     "Specific Total Energy", Wx, Wy, Ww, Wh);
      Wx = 0;
      Wy += offy;

      switch(problem)
      {
         case 8:
         case 7:
         case 6:
         case 5:
         case 4: // Noh Problem
         case 3:
         case 2:
         case 1:
         case 0:
         {
            // Compute errors
            ParGridFunction *rho_ex = new ParGridFunction(rho_gf.ParFESpace());
            ParGridFunction *vel_ex = new ParGridFunction(v_gf.ParFESpace());
            ParGridFunction *ste_ex = new ParGridFunction(ste_gf.ParFESpace());

            rho_ex->ProjectCoefficient(rho_coeff);
            vel_ex->ProjectCoefficient(v_coeff);
            ste_ex->ProjectCoefficient(ste_coeff);

            ParGridFunction rho_err(rho_gf), vel_err(v_gf), ste_err(ste_gf);
            rho_err -= *rho_ex;
            vel_err -= *vel_ex;
            ste_err -= *ste_ex;

            // Visualize difference between exact and approx
            VisualizeField(vis_rho_ex, vishost, visport, *rho_ex,
                           "Exact: Density", Wx, Wy, Ww, Wh);
            
            Wx += offx;
            VisualizeField(vis_v_ex, vishost, visport, *vel_ex,
                           "Exact: Velocity", Wx, Wy, Ww, Wh);
            
            Wx += offx;
            VisualizeField(vis_ste_ex, vishost, visport, *ste_ex,
                           "Exact: Specific Total Energy", Wx, Wy, Ww, Wh);
            Wx = 0;
            Wy += offy;

            // Visualize difference between exact and approx
            VisualizeField(vis_rho_err, vishost, visport, rho_err,
                           "Error: Density", Wx, Wy, Ww, Wh);
            
            Wx += offx;
            VisualizeField(vis_v_err, vishost, visport, vel_err,
                           "Error: Velocity", Wx, Wy, Ww, Wh);
            
            Wx += offx;
            VisualizeField(vis_ste_err, vishost, visport, ste_err,
                           "Error: Specific Total Energy", Wx, Wy, Ww, Wh);

            delete rho_ex;
            delete vel_ex;
            delete ste_ex;
         }
      }
   }

   // Perform the time-integration by looping over time iterations
   // ti with a time step dt.  The main function call here is the
   // LagrangianLOOperator.MakeTimeStep() funcall.
   double t = 0.0, t_old;
   if (convergence_testing)
   {
      dt = hmin; // Set timestep to smalled h val for convergence testing
   }

   bool last_step = false;
   int steps = 0;
   BlockVector S_old(S);

   cout << "Entering time loop\n";

   for (int ti = 1; !last_step; ti++)
   {
      /* Check if we need to change CFL */
      if (CompileTimeVals::change_CFL && problem == 7 && t > CompileTimeVals::CFL_time_change && hydro.GetCFL() != CompileTimeVals::CFL_second)
      {
         cout << "Changing CFL for Saltzman at time = " << t << endl;
         double CFL_new = CompileTimeVals::CFL_second;
         hydro.SetCFL(CFL_new);
      }

      if (optimize_timestep)
      {
         hydro.CalculateTimestep(S);
         dt = hydro.GetTimestep();
      }

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
      
      if (S_old.GetData() == S.GetData()) { cout << "\t State has not changed with step.\n"; return -1; }

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
         cout << "Current time: " << t << endl;
         double lnorm = ste_gf * ste_gf, norm;
         MPI_Allreduce(&lnorm, &norm, 1, MPI_DOUBLE, MPI_SUM, pmesh->GetComm());
 
         if (Mpi::Root())
         {
            const double sqrt_norm = sqrt(norm);

            cout << std::fixed;
            cout << "step " << std::setw(5) << ti
                 << ",\tt = " << std::setw(5) << std::setprecision(4) << t
                 << ",\tdt = " << std::setw(5) << std::setprecision(6) << dt
                 << ",\t|e| = " << std::setprecision(10) << std::scientific
                 << sqrt_norm;
            cout << endl;
         }

         if (visualization)
         {
            // Compute Density
            for (int i = 0; i < sv_gf.Size(); i++)
            {
               rho_gf[i] = 1./sv_gf[i];
            } 
            
            int Wx = 0, Wy = 0; // window position
            int Ww = 350, Wh = 350; // window size
            int offx = Ww+10, offy = Wh + 45; // window offsets

            VisualizeField(vis_rho, vishost, visport, rho_gf,
                           "Density", Wx, Wy, Ww,Wh); 
            Wx += offx;

            VisualizeField(vis_v, vishost, visport,
                           v_gf, "Velocity", Wx, Wy, Ww, Wh);
            Wx += offx;

            VisualizeField(vis_ste, vishost, visport, ste_gf,
                           "Specific Total Energy",
                           Wx, Wy, Ww,Wh);
            
            Wx = 0;
            Wy += offy;
            
            switch(problem)
            {
               case 8:
               case 7:
               case 6:
               case 5:
               case 4: // Noh Problem
               case 3:
               case 2:
               case 1:
               case 0:
               {
                  // Compute errors
                  ParGridFunction *rho_ex = new ParGridFunction(rho_gf.ParFESpace());
                  ParGridFunction *vel_ex = new ParGridFunction(v_gf.ParFESpace());
                  ParGridFunction *ste_ex = new ParGridFunction(ste_gf.ParFESpace());

                  rho_coeff.SetTime(t);
                  v_coeff.SetTime(t);
                  ste_coeff.SetTime(t);

                  rho_ex->ProjectCoefficient(rho_coeff);
                  vel_ex->ProjectCoefficient(v_coeff);
                  ste_ex->ProjectCoefficient(ste_coeff);

                  ParGridFunction rho_err(rho_gf), vel_err(v_gf), ste_err(ste_gf);
                  rho_err -= *rho_ex;
                  vel_err -= *vel_ex;
                  ste_err -= *ste_ex;

                  // Visualize difference between exact and approx
                  VisualizeField(vis_rho_ex, vishost, visport, *rho_ex,
                                 "Exact: Density", Wx, Wy, Ww, Wh);
                  
                  Wx += offx;
                  VisualizeField(vis_v_ex, vishost, visport, *vel_ex,
                                 "Exact: Velocity", Wx, Wy, Ww, Wh);
                  
                  Wx += offx;
                  VisualizeField(vis_ste_ex, vishost, visport, *ste_ex,
                                 "Exact: Specific Total Energy", Wx, Wy, Ww, Wh);
                  Wx = 0;
                  Wy += offy;


                  // Visualize difference between exact and approx
                  VisualizeField(vis_rho_err, vishost, visport, rho_err,
                                 "Error: Density", Wx, Wy, Ww, Wh);
                  
                  Wx += offx;
                  VisualizeField(vis_v_err, vishost, visport, vel_err,
                                 "Error: Velocity", Wx, Wy, Ww, Wh);
                  
                  Wx += offx;

                  VisualizeField(vis_ste_err, vishost, visport, ste_err,
                                 "Error: Specific Total Energy", Wx, Wy, Ww, Wh);

                  delete rho_ex;
                  delete vel_ex;
                  delete ste_ex;
               }
            }
         }
      }
   } // End time step iteration

   /* Plots end y velocity in the case of the Saltzman problem */
   if (visualization)
   {
      switch (problem)
      {
         case 7:
         {
            /* Prepare data */
            // Velocity plot
            int nv_py = pmesh->GetNV();
            std::vector<double> xgf_py(nv_py), mv_y_py(nv_py);
            // for(int i = 0; i < n_py; i++) {
            //    x.at(i) = i*i;
            //    y.at(i) = sin(2*M_PI*i/360.0);
            //    z.at(i) = log(i);
            // }
            for(int i = 0; i < nv_py; i++) {
               xgf_py[i] = x_gf[i];
               // mv_y_py[i] = mv_gf[i];
               mv_y_py[i] = mv_gf[i + H1FESpace.GetNDofs()];
            }

            // Density plot 
            int nc_py = L2FESpace.GetNE();
            cout << "nc_py: " << nc_py << endl;
            Vector center(dim);
            std::vector<double> rho_x_py(nc_py), rho_py(nc_py);
            double _min = 0.6, _max = 1.;
            for (int ci = 0; ci < nc_py; ci++) // cell iterator
            {
               pmesh->GetElementCenter(ci, center);
               rho_x_py[ci] = center[0];
               rho_py[ci] = rho_gf[ci];
            }

            std::vector<double> rho_x_exact_py(200), rho_exact_py(200);
            double _h = (_max - _min) / 200.;
            for (int i = 0; i < 200; i++)
            {
               rho_x_exact_py[i] = _min + i*_h;
               Vector _x(dim);
               _x = 0.;
               _x[0] = rho_x_exact_py[i];
               rho_exact_py[i] = InitialValues<problem, dim>::rho0(_x, t);
            }

            
            // Set the size of output image = 1200x780 pixels
            plt::figure_size(1200, 780);

            // Plot line from given x and y data. Color is selected automatically.
            // plt::scatter(xgf_py, mv_y_py);
            // plt::scatter(rho_x_py, rho_py);

            const long nrows=1, ncols=2;
            long row = 0, col = 0;

            plt::subplot2grid(nrows, ncols, row, col);
            plt::scatter(rho_x_py, rho_py);
            plt::named_plot("Exact", rho_x_exact_py, rho_exact_py, "black");
            plt::legend();
            plt::title("Density");

            col = 1;
            plt::subplot2grid(nrows, ncols, row, col);
            plt::scatter(xgf_py, mv_y_py);
            plt::title("Y Velocity");

            // Plot a red dashed line from given x and y data.
            // plt::plot(x, w,"r--");

            // Plot a line whose name will show up as "log(x)" in the legend.
            // plt::named_plot("log(x)", x, z);

            // Enable legend.
            plt::show();

            // save figure
            // const char* filename = "./basic.png";
            // std::cout << "Saving result to " << filename << std::endl;;
            // plt::save(filename);
            break;
         }
         default: // do nothing for all other problems
         {
         }
      }
   }

   // Print moved mesh
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      ostringstream mesh_name;
      mesh_name << "./results/mesh-test-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   // Last time recomputing density
   for (int i = 0; i < sv_gf.Size(); i++)
   {
      rho_gf[i] = 1./sv_gf[i];
   } 

   /* When the exact solution is known, print out an error file */
   switch(problem)
   {
      // case 8: // 8TODO
      case 7:
      case 6:
      case 4: // Noh problem
      case 3:
      case 2:
      {
         // Compute errors
         ParGridFunction *rho_ex = new ParGridFunction(rho_gf.ParFESpace());
         ParGridFunction *vel_ex = new ParGridFunction(v_gf.ParFESpace());
         ParGridFunction *ste_ex = new ParGridFunction(ste_gf.ParFESpace());

         rho_coeff.SetTime(t);
         v_coeff.SetTime(t);
         ste_coeff.SetTime(t);

         rho_ex->ProjectCoefficient(rho_coeff);
         vel_ex->ProjectCoefficient(v_coeff);
         ste_ex->ProjectCoefficient(ste_coeff);

         /* Compute denoms */
         double rho_ex_L1_norm = rho_ex->ComputeL1Error(zero);
         double rho_ex_L2_norm = rho_ex->ComputeL2Error(zero);
         double rho_ex_max_norm = rho_ex->ComputeMaxError(zero);

         double vel_ex_L1_norm = vel_ex->ComputeL1Error(zero);
         double vel_ex_L2_norm = vel_ex->ComputeL2Error(zero);
         double vel_ex_max_norm = vel_ex->ComputeMaxError(zero);

         double ste_ex_L1_norm = ste_ex->ComputeL1Error(zero);
         double ste_ex_L2_norm = ste_ex->ComputeL2Error(zero);
         double ste_ex_max_norm = ste_ex->ComputeMaxError(zero);
         
         /* Compute relative errors */
         double rho_L1_error_n = rho_gf.ComputeL1Error(rho_coeff) / rho_ex_L1_norm;
         double vel_L1_error_n = v_gf.ComputeL1Error(v_coeff) / vel_ex_L1_norm;
         double ste_L1_error_n = ste_gf.ComputeL1Error(ste_coeff) / ste_ex_L1_norm;
         const double L1_error = rho_L1_error_n + vel_L1_error_n + ste_L1_error_n;
 
         double rho_L2_error_n = rho_gf.ComputeL2Error(rho_coeff) / rho_ex_L2_norm;
         double vel_L2_error_n = v_gf.ComputeL2Error(v_coeff) / vel_ex_L2_norm;
         double ste_L2_error_n = ste_gf.ComputeL2Error(ste_coeff) / ste_ex_L2_norm;
         const double L2_error = rho_L2_error_n + vel_L2_error_n + ste_L2_error_n;

         double rho_Max_error_n = rho_gf.ComputeMaxError(rho_coeff) / rho_ex_max_norm;
         double vel_Max_error_n = v_gf.ComputeMaxError(v_coeff) / vel_ex_max_norm;
         double ste_Max_error_n = ste_gf.ComputeMaxError(ste_coeff) / ste_ex_max_norm;
         const double Max_error = rho_Max_error_n + vel_Max_error_n + ste_Max_error_n;

         if (Mpi::Root())
         {
            ostringstream convergence_filename;
            convergence_filename << "/Users/madisonsheridan/Workspace/Laglos/saved/convergence/temp_output/np" << num_procs;
            if (rs_levels != 0) {
               convergence_filename << "_s" << setfill('0') << setw(2) << rs_levels;
            }
            if (rp_levels != 0) {
               convergence_filename << "_p" << setfill('0') << setw(2) << rp_levels;
            }
            convergence_filename << "_refinement_"
                                 << setfill('0') << setw(2)
                                 << to_string(rp_levels + rs_levels)
                                 << ".out";
            ofstream convergence_file(convergence_filename.str().c_str());
            convergence_file.precision(precision);
            convergence_file << "Processor_Runtime " << "1." << "\n"
                             << "n_processes " << num_procs << "\n"
                             << "n_refinements "
                             << to_string(rp_levels + rs_levels) << "\n"
                             << "n_Dofs " << global_vSize << "\n"
                             << "h " << hmin << "\n"
                             << "L1_Error " << L1_error << "\n"
                             << "L2_Error " << L2_error << "\n"
                             << "Linf_Error " << Max_error << "\n"
                             << "mass_loss " << hydro.CalcMassLoss(S) << "\n"
                             << "dt " << dt << "\n"
                             << "Endtime " << t << "\n";
                        
            convergence_file.close();

            delete rho_ex;
            delete vel_ex;
            delete ste_ex;
         }
      }
      default: // do nothing for all other problems
      {
      }
   }
   
   delete pmesh;
   delete m;

   return 0;
}
