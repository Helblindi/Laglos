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

/***
* --- To generate images for movies
* -of "output-dir" = directory where files should be saved.
* -print           = flag that ensures the files are saved to the designated directory
* -vs #            = number of steps between each visualization or saved image
*
* --- To generate movies
* python create_glvis_animation_script.py "/Users/madisonsheridan/Workspace/Laglos/build/vortex/gfs_r03/"
* glvis -run 
*
* --- Current work in progress, new mesh velocity runs ---
* Current issue is singular matrix on boundary nodes
* Sod, smooth, 
* ./Laglos -m ../data/ref-square.mesh -p 1 -tf 0.225 -cfl 0.5 -rs 4 -mv 3
*
* CAVEAT Weight LS method for corner vertex calculation 
* Sod
* ./Laglos -m ../data/ref-square.mesh -p 1 -tf 0.225 -cfl 0.5 -rs 4 -mv 4 
*
* Test problems:
*     p = 0  --> taylor green (smooth problem)
*     p = 1  --> sedov
*     p = 2  --> sod
*     p = 3  --> tp
*     p = 4  --> gresho vortex (smooth problem)
*     p = 5  --> 2D Riemann problem, config. 12 of doi.org/10.1002/num.10025
*     p = 6  --> 2D Riemann problem, config.  6 of doi.org/10.1002/num.10025
*     p = 7  --> 2D Rayleigh-Taylor instability problem.
*     p = 8  --> sod radial
*     p = 9  --> isentropic vortex
*     p = 10 --> noh
*     p = 11 --> saltzman
*     p = 12 --> vdw1
*     p = 13 --> vdw2
*     p = 14 --> vdw3
*     p = 15 --> vdw4
*     p = 16 --> Kidder shell
*     p = 17 --> kidder ball
*     p = 18 --> ICF
*     p = 21 --> Sedov (Laura)
*     p = 40 --> smooth wave
*     p = 41 --> lax
*     p = 42 --> leblanc
*     p = 43 --> riemann
*     p = 50 --> elastic shocktube
*     p = 51 --> elastic impact
*     p = 52 --> elastic shear
*     p = 53 --> elastic human leg impact (Bonet-Burton 1997)
*     p = 54 --> elastic projectile plate
*     p = 55 --> elastic shear in y direction
*     p = 56 --> elastic impact with shear
*     p = 57 --> elastic rotation (Favrie 2014 Sec. 5.4)
*     p = 58 --> elast Noh (Experimental)
*     p = 59 --> elastic projectile impact (Vilar-Maire-Shu 2D)
*     p = 60 --> elastic test
*     p = 100 --> TestBCs 
***/
#include "mfem.hpp"
#include "var-config.h"
#include "compile_time_vals.h"
#include "laglos_solver.hpp"
#include "test_problems_include.h"
#include "laglos_tools.hpp"
#include "riemann1D.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <boost/filesystem.hpp>
#include <sstream>
#include <string>

#define _USE_MATH_DEFINES
#include <cmath>

#ifdef DMALLOC  // Memory allocation checks
#include "dmalloc.h"
#endif

using namespace mfem;
using namespace hydroLO;
using namespace std;

/* forward declarations */
static void display_banner(std::ostream&);


int main(int argc, char *argv[]) {
   StopWatch chrono;

   // Initialize MPI.
   Mpi::Init();
   const int num_procs = Mpi::WorldSize();
   const int myid = Mpi::WorldRank();
   // Hypre::Init();

   // Print the banner.
   if (Mpi::Root()) { display_banner(cout); }

   int dim = 2;
   const string results_dir = CompileTimeVals::results_dir;

   // Parse command line options
   const char *mesh_file_location = "default";
   const char *output_flag = "default";
   int problem = 0;
   int rs_levels = 0;
   int rp_levels = 0;
   int order_k = 1;  // Order of mesh movement approximation space
   int order_t = 0;
   int ode_solver_type = 1;
   double t_init = 0.0;
   double t_final = 0.6;
   int max_tsteps = -1;
   double dt = 0.001;
   bool visualization = false;
   int vis_steps = 5;
   bool pview = false;
   bool gfprint = false;
   int precision = 12;
   /* This parameter describes how often to compute mass corrective face velocities, 0 indicates to always take the average. */
   bool mc = false;
   int visc = 1; // Default to GMS-GV
   bool greedy = false;
   int gmv_to_greedy_steps = 0;
   bool mm = true;
   bool check_mesh = true;
   bool post_process_density = true;
   int elastic_eos = 0;
   int mv_option = 2;
   double mv_target_visc_coeff = 0.;
   bool do_mv_linearization = false;
   int fv_option = 0;
   int mv_it_option = 0;
   int mv_n_iterations = 0;
   double mm_visc_face = 0., mm_cell = 0.;
   bool convergence_testing = false;
   bool suppress_output = false;
   double CFL = 0.5;
   double dm_val = 0.; // Parameter to distort the mesh, overrides file val
   string sv_output_prefix;
   string ts_output_prefix;
   string output_path;

   // Debug options
   bool debug_mesh_velocity = false;

   OptionsParser args(argc, argv);

   args.AddOption(&mesh_file_location, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&output_flag, "-of", "--output-flag", 
                  "Directory that output files should be placed");
   args.AddOption(&order_t, "-ot", "--order-thermo",
                  "Order (degree) of the thermodynamic finite element space.");
   args.AddOption(&order_k, "-ok", "--order-kinetic",
                  "Order (degree) of the mesh velocity finite element space.");
   args.AddOption(&rs_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&rp_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&problem, "-p", "--problem", "Test problem to run.");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  "ODE solver: 1 - Forward Euler,\n\t"
                  "            2 - RK2 SSP, 3 - RK3 SSP, 4 - RK4, 6 - RK6,\n\t"
                  "            7 - RK2Avg.");
   args.AddOption(&t_init, "-ti", "--t-init",
                  "Start time; default is 0.");               
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&max_tsteps, "-ms", "--max-steps",
                  "Maximum number of steps (negative means no restriction).");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.AddOption(&pview, "-pview", "--paraview", "-no-pview", "--no-paraview",
                  "Enable or disable ParaView visualization.");
   args.AddOption(&gfprint, "-print", "--print", "-no-print", "--no-print",
                  "Enable or disable result output (files in mfem format).");
   args.AddOption(&mc, "-mc", "--mass-correct", "-no-mc", "--no-mass-correct",
                  "Enable or disable mass correction.");
   args.AddOption(&visc, "-visc", "--viscosity",
                  "Type of viscosity to use:"
                  "\n\t 0 - No viscosity,"
                  "\n\t 1 - GMS-GV (Guaranteed Maximum Speed Graph Viscosity),"
                  "\n\t 2 - Greedy Viscosity,"
                  "\n\t 3 - Artificial graph viscosity for HO (Binder),"
                  "\n\t 4 - Artificial graph viscosity for HO (Non-binder),");
   args.AddOption(&greedy, "-greedy", "--use-greedy-viscosity", "-no-greedy", "--no-greedy-viscosity",
                  "Enable or disable greedy viscosity computation.");
   args.AddOption(&gmv_to_greedy_steps, "-ggn", "--gmv-to-greedy-num-steps",
                  "Number of steps of GMV to take before switching to greedy viscosity.");
   args.AddOption(&mm, "-mm", "--move-mesh", "-no-mm", "--no-move-mesh",
                  "Enable or disable mesh movement.");
   args.AddOption(&check_mesh, "-cm", "--check-mesh", "-no-cm", "--no-check-mesh",
                  "Enable or disable checking if the mesh has twisted.");
   args.AddOption(&post_process_density, "-ppd", "--post-process-density", "-no-ppd", "--no-post-process-density",
                  "Enable or disable density post processing to guarantee conservation of mass.");
   args.AddOption(&elastic_eos, "-ue", "--use-elasticity",
                  "Choose equation of state for shear energy:"
                  "\n\t 00 - No shear energy,"
                  "\n\t 01 - Neo Hookean EOS,"
                  "\n\t 02 - Mooney Rivlin EOS,"
                  "\n\t 03 - Aortic EOS,"
                  "\n\t 04 - Transversely Isotropic EOS"
                  "\n\t 05 - Multiple shear EOS");
   args.AddOption(&mv_option, "-mv", "--mesh-velocity-option",
                  "Choose how to compute mesh velocities:"
                  "\n\t 00 - Arithmetic avg of adj cells,"
                  "\n\t 01 - Arithmetic avg of adj cells with distributed viscosity,"
                  "\n\t 02 - Cell Face Normal with viscosity,"
                  "\n\t 1* - Sparse HiOp LM with * mv as target"
                  );
                  // "\n\t 4 - CAVEAT Weighted LS,"
                  // "\n\t 5 - CAVEAT/Cell Face combo,"
                  // "\n\t 6 - 5 with weights");
   args.AddOption(&mv_target_visc_coeff, "-tvc", "--mv-target-viscosity-coeff",
                  "Set the coefficient for the viscosity term in the target mesh velocity approach.");
   args.AddOption(&do_mv_linearization, "-mv-lin", "--mesh-velocity-linearization",
                  "-no-mv-lin", "--no-mesh-velocity-linearization",
                  "Enable mesh velocity linearization using alpha parameter.");
   args.AddOption(&fv_option, "-fv", "--face-velocity-option",
                  "Choose how to compute face velocities:"
                  "\n\t 0 - Do nothing,"
                  "\n\t 1 - Mass conservative bubble, Q2"
                  "\n\t 2 - Average, Q1 type,"
                  "\n\t 3 - Butterfly, Q2");
   args.AddOption(&mv_it_option, "-mv-it-op", "--mv-it-op",
                  "Set the desired type of mesh velocity iteration to use");
   args.AddOption(&mv_n_iterations, "-mv-iter-n", "--mv-num-iterations",
                  "Set the number of times to iterate on the corner node mesh velocities.");
   args.AddOption(&mm_visc_face, "-mmvf", "--mesh-motion-viscosity-coefficient-on-face",
                  "Set the amount of viscosity to add to the mesh motion iteration on the faces.");
   args.AddOption(&mm_cell, "-mmcc", "--mesh-motion-consistency-coefficient-on-cells",
                  "Set the coefficient for the consistency term defined on the adjacent cells.");
   args.AddOption(&CFL, "-cfl", "--CFL",
                  "CFL value to use.");
   args.AddOption(&convergence_testing, "-ct", "--set-dt-to-h", "-no-ct", 
                  "--no-set-dt-to-h", 
                  "Enable or disable convergence testing timestep.");
   args.AddOption(&dt, "-dt", "--timestep", "Timestep to use.");
   args.AddOption(&suppress_output, "-so", "--suppress-output", "-no-so",
                  "--no-suppress-output",
                  "Enable or disable output during runtime.");
   args.AddOption(&dm_val, "-dm", "--distort-mesh",
                  "Mesh distortion parameter.");
   args.AddOption(&debug_mesh_velocity, "-dmv", "--debug-mesh-velocity", "-no-dmv",
                  "--no-debug-mesh-velocity",
                  "Enable or disable debug mesh velocity visualization in glvis.");

   args.Parse();
   if (!args.Good())
   {
      if (Mpi::Root()) { args.PrintUsage(cout); }
      return 1;
   }
   if (Mpi::Root()) { args.PrintOptions(cout); }

   /***** Arg checks *****/
   if (elastic_eos > 0 && order_t > 0)
   {
      MFEM_ABORT("Elasticity cannot be used with order_t > 0. "
                 "Set order_t = 0 to use elasticity.");
   }
   // if (visc == 2)
   // {
   //    if (!greedy)
   //    {
   //       MFEM_ABORT("Greedy viscosity must be enabled if visc = 2.");
   //    }
   // }
   // if (greedy && visc != 2)
   // {
   //    MFEM_ABORT("Greedy viscosity can only be used if visc = 2.");
   // }
   // if (order_t > 0 && order_k != order_t)
   // {
   //    MFEM_ABORT("In the case of higher order approximations, the mesh velocity order must be equal to that of the thermodynamic space.");
   // }
   if (fv_option > 0 && order_k != order_t + 2)
   {
      MFEM_ABORT("If a different face velocity is to be used, the mesh velocity order must be equal to the thermodynamic order + 2.");
   }
   /**** End arg checks *****/

   // Set output_flag string
   if (strncmp(output_flag, "default", 7) != 0)
   {
      output_path = std::string(output_flag) + "/";
   }
   else
   {
      output_path = results_dir + "temp/";
   }

   // We must manually create all corresponding output directories
   ostringstream _refinement_path_ss;
   _refinement_path_ss << output_path
                       << "r"
                       << setfill('0')
                       << setw(2)
                       << to_string(rp_levels + rs_levels)
                       << "/";
   string _refinement_path = _refinement_path_ss.str();

   string _convergence = output_path + "convergence";
   string _temp_output = _convergence + "/temp_output";
   string _state_vectors = output_path + "state_vectors";
   string gfprint_path = _refinement_path + "gfs/";

   const char* path = output_path.c_str();
   boost::filesystem::path output_path_dir(path);
   boost::filesystem::create_directory(output_path_dir);

   path = _convergence.c_str();
   boost::filesystem::path convergence_path_dir(path);
   boost::filesystem::create_directory(convergence_path_dir);

   path = _temp_output.c_str();
   boost::filesystem::path temp_output_dir(path);
   boost::filesystem::create_directory(temp_output_dir);

   sv_output_prefix = output_path + "state_vectors/";
   path = _state_vectors.c_str();
   boost::filesystem::path state_vectors_dir(path);
   boost::filesystem::create_directory(state_vectors_dir);

   ts_output_prefix = output_path + "timeseries/";
   path = ts_output_prefix.c_str();
   boost::filesystem::path ts_dir(path);
   boost::filesystem::create_directory(ts_dir);

   path = _refinement_path.c_str();
   boost::filesystem::path refinement_dir(path);
   boost::filesystem::create_directory(refinement_dir);

   path = gfprint_path.c_str();
   boost::filesystem::path gfprint_output_dir(path);
   boost::filesystem::create_directory(gfprint_output_dir);

   // On all processors, use the default builtin 1D/2D/3D mesh or read the
   // serial one given on the command line.
   Mesh *mesh;
   if (strncmp(mesh_file_location, "default", 7) != 0)
   {
      mesh = new Mesh(mesh_file_location, true, true);
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

   dim = mesh->Dimension();
   cout << "Meshing done\n";

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
   ParMesh *pmesh0 = new ParMesh(MPI_COMM_WORLD, *mesh);
   ParMesh *smesh = new ParMesh(ParMesh::MakeRefined(*pmesh, order_t+1, BasisType::ClosedUniform));;
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

   // Set up problem
   ProblemBase * problem_class = NULL;
   switch (problem)
   {
      case 0: // Taylor-Green
         problem_class = new TaylorGreenProblem(dim);
         break;
      case 1: // Sedov
         problem_class = new SedovLLNLProblem(dim);
         break;
      case 2: // Sod
         problem_class = new SodProblem(dim);
         break;
      case 3: // Triple Point
         problem_class = new TriplePoint(dim);
         break;
      case 4: // gresh vortex
         break;
      case 5: // 2D Riemann problem
         break;
      case 6: // 2D Riemann problem
         break;
      case 7: // 2D Rayleigh-Taylor instability
         break;
      case 8: // Radial Sod
         problem_class = new SodRadial(dim);
         break;
      case 9: // Isentropic Vortex, stationary center
         problem_class = new IsentropicVortex(dim);
         break;
      case 10: // Noh
         problem_class = new NohProblem(dim);
         break;
      case 11: // Saltzmann
         problem_class = new SaltzmannProblem(dim);
         break;
      /* VDW */
      case 12:
         problem_class = new VdwTest1(dim);
         break;
      case 13:
         problem_class = new VdwTest2(dim);
         break;
      case 14:
         problem_class = new VdwTest3(dim);
         break;
      case 15:
         problem_class = new VdwTest4(dim);
         break;
      case 16: // Kidder shell
         problem_class = new KidderProblem(dim);
         break;
      case 17: // Kidder ball
         problem_class = new KidderBallProblem(dim);
         break;
      case 18: // ICF
         problem_class = new ICFProblem(dim);
         break;
      case 19: // VDW Isentropic Vortex
         problem_class = new VDWIsentropicVortex(dim);
         break;
      case 21: // Sedov
      {
         MFEM_ABORT("Not implemented\n");
         assert(hmin == hmax);
         Vector params(2);
         params[0] = hmax, params[1] = pmesh->GetElementVolume(0);

         problem_class = new SedovProblem(dim);
         problem_class->update(params, t_init);
         // TODO: Will need to modify initialization of internal energy
         //       if distorted meshes are used.
         break;
      }
      case 40: // Smooth
         problem_class = new SmoothWave(dim);
         break;
      case 41: // Lax
         problem_class = new LaxProblem(dim);
         break;
      case 42: // Leblanc
         problem_class = new LeblancProblem(dim);
         break;
      case 43: // Riemann Problem
         problem_class = new RiemannProblem(dim);
         break;
      case 50: // Elastic shocktube
         problem_class = new ElasticShocktube(dim);
         break;
      case 51: // Elastic impact
         problem_class = new ElasticImpact(dim);
         break;
      case 52: // Elastic shear
         problem_class = new ElasticShear(dim);
         break;
      case 53: // Elastic human leg impact
         problem_class = new ElasticHumanLegImpact(dim);
         break;
      case 54: // Elastic projectile plate
         problem_class = new ElasticProjectilePlate(dim);
         break;
      case 55: // Elastic shear rotate in y direction
         problem_class = new ElasticShearY(dim);
         break;
      case 56: // Elastic impact + shear
         problem_class = new ElasticImpactShear(dim);
         break;
      case 57: // Elastic 2D, Favrie 2014 section 5.4
         problem_class = new ElasticRotation(dim);
         break;
      case 58: // Elastic noh
         problem_class = new ElasticNoh(dim);
         break;
      case 59: // Elastic projectile impact (vilar-mair-shu 2d)
         problem_class = new ElasticProjectileImpact(dim);
         break;
      case 60: // Elastic test
         problem_class = new ElasticTest(dim);
         break;
      case 100:
         problem_class = new TestBCs(dim);
         break;
      default:
         MFEM_ABORT("Failed to initiate a problem.\n");
   }

   // Define the parallel finite element spaces. We use:
   // - H1 (continuous) for mesh movement. (Q1 / order_k)
   // - L2 (discontinuous) for state variables (Q0 / order_t)
   // - CR/RT for mesh velocity reconstruction at nodes
   H1_FECollection H1FEC(order_k, dim);
   H1_FECollection H1FEC_L(1, dim);
   L2_FECollection L2FEC(order_t, dim, BasisType::GaussLobatto);
   L2_FECollection L20FEC(0, dim);
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
   ParFiniteElementSpace L20FESpace(pmesh, &L20FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);
   ParFiniteElementSpace CRFESpace(smesh, CRFEC, dim);

   // Define the explicit ODE solver used for time integration.
   ODESolver *ode_solver = NULL;
   switch (ode_solver_type)
   {
      case 1: ode_solver = new ForwardEulerSolver; break;
      case 2: ode_solver = new RK2AvgSolver; break;
      case 3: ode_solver = new HydroRK3SSPSolver; break;
      // case 2: ode_solver = new RK2Solver(0.5); break;
      // case 3: ode_solver = new RK3SSPSolver; break;
      // case 4: ode_solver = new RK4Solver; break;
      // case 6: ode_solver = new RK6Solver; break;
      // case 7: ode_solver = new RK2AvgSolver; break;
      case 11: ode_solver = new ForwardEulerIDPSolver; break;
      case 12: ode_solver = new RK2IDPSolver; break;
      case 13: ode_solver = new RK3IDPSolver; break;
      case 14: ode_solver = new RK4IDPSolver; break;
      case 16: ode_solver = new RK6IDPSolver; break;
      default:
         MFEM_ABORT("ode_solver_type not implemented\n");
         // if (myid == 0)
         // {
         //    cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         // }
         // delete pmesh;
         // MPI_Finalize();
         // return 3;
   }

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

   /* Print element information */
   cout << "num boundary elements: " << pmesh->GetNBE() << endl;
   cout << "num elements: " << pmesh->GetNE() << endl;

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
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   sv_gf.MakeRef(&L2FESpace, S, offset[1]);
   v_gf.MakeRef(&L2VFESpace, S, offset[2]);
   ste_gf.MakeRef(&L2FESpace, S, offset[3]);

   // Initialize x_gf using starting mesh positions
   pmesh->SetNodalGridFunction(&x_gf);

   /* Distort Mesh for Saltzman Problem */
   if (problem_class->get_distort_mesh() || dm_val != 0.)
   {
      cout << "distort mesh.\n";
      switch (problem)
      {
      case 11: // Saltzmann
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
         break;
      }
      default: // All others
      {
         double mesh_distortion_factor = 0.1;
         // Random distortion based on min mesh size
         Array<double> coords(dim);
         for (int vertex = 0; vertex < H1FESpace.GetNDofs(); vertex++)
         {
            for (int i = 0; i < dim; i++)
            {
               int index = vertex + i * H1FESpace.GetNDofs();
               coords[i] = x_gf[index];

               // Generate a random number between -1 and 1
               double random_number = mesh_distortion_factor * (-0.5 + (1.0 * std::rand()) / RAND_MAX);

               // Calculate new coord and store it
               double coord_new = coords[i] + random_number * hmin;
               x_gf[index] = coord_new;
            }
         }
         break;
      }
      } // switch

      // Update pmesh reference grid function
      pmesh->NewNodes(x_gf, false);
      cout << "Mesh distorted\n";
   } // End distort mesh

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

   /* Specific volume, velocity, and energy initialization */
   // Similar to Laghos, we interpolate in a non-positive basis to get
   // the correct values at the dofs. Then we do an L2 projection to 
   // the positive basis in which we actually compute.
   // L2_FECollection l2_fec(order_t, pmesh->Dimension(), BasisType::GaussLegendre);
   // ParFiniteElementSpace l2_fes(pmesh, &l2_fec);
   // ParFiniteElementSpace l2_vfes(pmesh, &l2_fec, dim);
   // ParGridFunction l2_ste(&l2_fes), l2_sv(&l2_fes), l2_v(&l2_vfes);

   /* sv */
   FunctionCoefficient sv_coeff(sv0_static);
   sv_coeff.SetTime(t_init);
   // l2_sv.ProjectCoefficient(sv_coeff);
   // sv_gf.ProjectGridFunction(l2_sv);
   sv_gf.ProjectCoefficient(sv_coeff);
   sv_gf.SyncAliasMemory(S);

   /* v */
   VectorFunctionCoefficient v_coeff(dim, v0_static);
   v_coeff.SetTime(t_init);
   // l2_v.ProjectCoefficient(v_coeff);
   // v_gf.ProjectGridFunction(l2_v);
   v_gf.ProjectCoefficient(v_coeff);
   v_gf.SyncAliasMemory(S);

   /* ste */
   // While the ste_coeff is not used for initialization in the Sedov case,
   // it is necessary for plotting the exact solution
   FunctionCoefficient ste_coeff(ste0_static);
   ste_coeff.SetTime(t_init);
   if (problem == 1)
   {
      double blast_energy = 0.25;
      double blast_position[] = {0.0, 0.0, 0.0};
      DeltaCoefficient e_coeff(blast_position[0], blast_position[1],
                               blast_position[2], blast_energy);
      // l2_ste.ProjectCoefficient(e_coeff);
      ste_gf.ProjectCoefficient(e_coeff);
   }
   else
   {
      // l2_ste.ProjectCoefficient(ste_coeff);
      ste_gf.ProjectCoefficient(ste_coeff);
   }
   // ste_gf.ProjectGridFunction(l2_ste);
   ste_gf.SyncAliasMemory(S);

   // PLF to build mass vector
   FunctionCoefficient rho0_coeff(rho0_static), rho_coeff(rho0_static);
   rho0_coeff.SetTime(t_init);
   int ir_order = 3 * H1FESpace.GetOrder(0) + L2FESpace.GetOrder(0) - 1;
   /* 
   Use default IntRules object, which uses Gauss-Legendre basis 
   For some reason, using Gauss-Lobatto for this basis prohibits
   the Sedov problem from running to completion.
   */
   IntegrationRule ir = IntRules.Get(pmesh->GetElementBaseGeometry(0), ir_order);
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(rho0_coeff,&ir));
   m->Assemble();

   /* Initialize rho0_gf */
   // ParGridFunction l2_rho0(&l2_fes), l2_rho(&l2_fes);
   ParGridFunction rho0_gf(&L2FESpace);
   // l2_rho0.ProjectCoefficient(rho0_coeff);
   // rho0_gf.ProjectGridFunction(l2_rho0);
   rho0_gf.ProjectCoefficient(rho0_coeff);

   FunctionCoefficient p_coeff(p0_static);
   p_coeff.SetTime(t_init);

   /* Create Lagrangian Low Order Solver Object */
   LagrangianLOOperator hydro(dim, S.Size(), H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, rho0_coeff, rho0_gf, m, ir, problem_class, offset, visc, elastic_eos, mm, CFL);
   hydro.SetInitialMassesAndVolumes(S);

   /* Set parameters of the LagrangianLOOperator */
   hydro.SetMVOption(mv_option);
   hydro.SetMVLinOption(do_mv_linearization);
   hydro.SetFVOption(fv_option);
   if (elastic_eos) //NF//MS
   {
      if (problem_class->get_shear_modulus() < 1.e-12)
      {
         MFEM_WARNING("Elasticity has been chosen, but the shear modulus is 0, meaning this is the fluid case.");
      }
   }

   /* 
   If opting to use greedy viscosity, verify that you have not opted to use a few steps of GMV first. 
   This options proves to be useful in the case of the Sedov problem, where the sie is initialized 
   with a delta function.
   */
   if (gmv_to_greedy_steps == 0 && greedy)
   {
      hydro.SetViscOption(greedy);
   }

   hydro.SetProblem(problem);
   hydro.SetDensityPP(post_process_density);
   hydro.SetMVTargetViscCoeff(mv_target_visc_coeff);
   if (mv_n_iterations != 0)
   {
      // If an iterative method is desired for computation of the mesh velocity
      // at the lagrange nodes, the type of iteration method must be specified.
      if (mv_it_option == 0)
      {
         MFEM_ABORT("Must set the mesh velocity iteration option.\n");
      }
      hydro.SetMVIterationOption(mv_it_option);
      hydro.SetMVIteration(mv_n_iterations);
      hydro.SetMMViscFace(mm_visc_face);
      hydro.SetMMCell(mm_cell);
   }

   /* Set up visualiztion object */
   socketstream vis_rho, vis_v, vis_ste, vis_press, vis_gamma, vis_mc;
   socketstream vis_sig, vis_f, vis_frho, vis_esheer;
   socketstream vis_rho_ex, vis_v_ex, vis_ste_ex, vis_p_ex;
   socketstream vis_rho_err, vis_v_err, vis_ste_err, vis_p_err;
   /* debug */
   socketstream vis_mv;
   char vishost[] = "localhost";
   int  visport   = 19916;

   ParGridFunction rho_gf(&L2FESpace), rho_cont_gf(&H1cFESpace);
   ParGridFunction mc_gf(&L20FESpace); // Gridfunction to show mass conservation
   mc_gf = 0.;  // if a cells value is 0, mass is conserved
   ParGridFunction mass_gf(&L20FESpace); // Gridfunction to show mass per cell
   ParGridFunction vol_gf(&L20FESpace); // Gridfunction to show volume per cell

   /* Gridfunctions used in Triple Point */
   ParGridFunction press_gf(&L2FESpace), gamma_gf(&L2FESpace);

   /* Elasticity gridfunctions */
   ParGridFunction sigma_gf(&L2FESpace), f_gf(&L2FESpace), frho_gf(&L2FESpace), e_sheer_gf(&L2FESpace);

   // Compute values for visualization
   if (visualization || gfprint || pview)
   {
      /* Density */
      for (int i = 0; i < sv_gf.Size(); i++)
      {
         rho_gf[i] = 1./sv_gf[i];
      } 
      hydro.MassesAndVolumesAtPosition(rho_gf, x_gf, mass_gf, vol_gf);

      /* Continuous density projection */
      GridFunctionCoefficient rho_gf_coeff(&rho_gf);
      rho_cont_gf.ProjectDiscCoefficient(rho_gf_coeff, mfem::ParGridFunction::AvgType::ARITHMETIC);

      /* Pressure */
      hydro.ComputePressGF(S, press_gf);

      /* Elastic attributes*/
      if (elastic_eos)
      {
         // Compute Sigma and F
         hydro.ComputeSigmaGF(S, sigma_gf);
         hydro.ComputeFGF(f_gf);
         hydro.ComputeESheerGF(e_sheer_gf);

         for (int i = 0; i < L2FESpace.GetNDofs(); i++)
         {
            frho_gf[i] = f_gf[i] * rho_gf[i];
         }
      }
   }

   /* Gamma */
   if (problem == 1 || problem == 3 || elastic_eos)
   {
      for (int i = 0; i < press_gf.Size(); i++)
      {
         int el_i = L2FESpace.GetElementForDof(i);
         gamma_gf[i] = problem_class->get_gamma(pmesh->GetAttribute(el_i));
      }
   }

   /* Initial visualization */
   if (visualization)
   {
      // Make sure all MPI ranks have sent their 'v' solution before initiating
      // another set of GLVis connections (one from each rank):
      MPI_Barrier(pmesh->GetComm());

      vis_rho.precision(8);
      vis_v.precision(8);
      vis_ste.precision(8);

      vis_rho_ex.precision(8);
      vis_v_ex.precision(8);
      vis_ste_ex.precision(8);
      vis_p_ex.precision(8);

      vis_rho_err.precision(8);
      vis_v_err.precision(8);
      vis_ste_err.precision(8);
      vis_p_err.precision(8);

      vis_press.precision(8);
      vis_gamma.precision(8);
      vis_mc.precision(8);

      vis_sig.precision(8);
      vis_f.precision(8);
      vis_frho.precision(8);
      vis_esheer.precision(8);

      /* Debug */
      vis_mv.precision(8);

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

      if (problem == 1 || problem == 3 || problem == 16 || elastic_eos)
      {
         Wx += offx;
         VisualizeField(vis_press, vishost, visport, press_gf,
                        "Pressure", Wx, Wy, Ww, Wh);
      }
      if (problem == 3)
      {
         Wx += offx; 
         VisualizeField(vis_gamma, vishost, visport, gamma_gf,
                        "Gamma", Wx, Wy, Ww, Wh);
      }
      if (debug_mesh_velocity)
      {
         ParGridFunction mv_gf(&H1FESpace);
         hydro.GetMV(mv_gf);
         Wx += offx;
         VisualizeField(vis_mv, vishost, visport, mv_gf,
                        "Mesh Velocity", Wx, Wy, Ww, Wh);
      }

      //NF//MS
      if (elastic_eos)
      {
         // Visualize
         Wx = 0;
         Wy += offy;

         VisualizeField(vis_sig, vishost, visport, sigma_gf,
                        "Sigma", Wx, Wy, Ww, Wh);
         Wx += offx;
         VisualizeField(vis_f, vishost, visport, f_gf,
                        "F", Wx, Wy, Ww, Wh);
         Wx += offx;
         VisualizeField(vis_frho, vishost, visport, frho_gf,
                        "F/rho", Wx, Wy, Ww, Wh);
         Wx += offx;
         VisualizeField(vis_esheer, vishost, visport, e_sheer_gf,
                        "e shear", Wx, Wy, Ww, Wh);
      }
      
      Wx = 0;
      Wy += offy;

      if (problem_class->has_exact_solution())
      {
         ParGridFunction rho_err(rho_gf), vel_err(v_gf), ste_err(ste_gf), p_err(press_gf);

         // Set current time for exact solution coefficient variables
         rho_coeff.SetTime(t_init);
         v_coeff.SetTime(t_init);
         ste_coeff.SetTime(t_init);
         p_coeff.SetTime(t_init);

         if (problem_class->get_indicator() == "Kidder")
         {
            ParFiniteElementSpace L2FESpace_kidder(pmesh0, &L2FEC);
            ParFiniteElementSpace L2VFESpace_kidder(pmesh0, &L2FEC, dim);

            ParGridFunction rho_k_ex_gf(&L2FESpace_kidder), vel_k_ex_gf(&L2VFESpace_kidder), ste_k_ex_gf(&L2FESpace_kidder), p_k_ex_gf(&L2FESpace_kidder);
            rho_k_ex_gf.ProjectCoefficient(rho_coeff);
            vel_k_ex_gf.ProjectCoefficient(v_coeff);
            ste_k_ex_gf.ProjectCoefficient(ste_coeff);
            p_k_ex_gf.ProjectCoefficient(p_coeff);

            rho_err -= rho_k_ex_gf;
            vel_err -= vel_k_ex_gf;
            ste_err -= ste_k_ex_gf;
            p_err -= p_k_ex_gf;
         }
         else
         {
            if (problem_class->get_indicator() == "Vdw1")
            {
               problem_class->update(x_gf);
            }

            // Compute errors
            ParGridFunction rho_ex_gf(&L2FESpace), vel_ex_gf(&L2VFESpace), ste_ex_gf(&L2FESpace), p_ex_gf(&L2FESpace);
            /* rho */
            // l2_rho.ProjectCoefficient(rho_coeff);
            // rho_ex_gf.ProjectGridFunction(l2_rho);
            rho_ex_gf.ProjectCoefficient(rho_coeff);
            /* v */
            // l2_v.ProjectCoefficient(v_coeff);
            // vel_ex_gf.ProjectGridFunction(l2_v);
            vel_ex_gf.ProjectCoefficient(v_coeff);
            /* ste */
            // l2_ste.ProjectCoefficient(ste_coeff);
            // ste_ex_gf.ProjectGridFunction(l2_ste);
            ste_ex_gf.ProjectCoefficient(ste_coeff);
            /* pressure */
            // ParGridFunction l2_p(&l2_fes);
            // l2_p.ProjectCoefficient(p_coeff);
            // p_ex_gf.ProjectGridFunction(l2_p);
            p_ex_gf.ProjectCoefficient(p_coeff);

            rho_err -= rho_ex_gf;
            vel_err -= vel_ex_gf;
            ste_err -= ste_ex_gf;
            p_err -= p_ex_gf;

            // Visualize difference between exact and approx
            VisualizeField(vis_rho_ex, vishost, visport, rho_ex_gf,
                           "Exact: Density", Wx, Wy, Ww, Wh);
            
            Wx += offx;
            VisualizeField(vis_v_ex, vishost, visport, vel_ex_gf,
                           "Exact: Velocity", Wx, Wy, Ww, Wh);
            
            Wx += offx;
            VisualizeField(vis_ste_ex, vishost, visport, ste_ex_gf,
                           "Exact: Specific Total Energy", Wx, Wy, Ww, Wh);

            if (problem == 1 || problem == 3)
            {
               Wx += offx;
               VisualizeField(vis_p_ex, vishost, visport, p_ex_gf,
                              "Exact: Pressure", Wx, Wy, Ww, Wh);
            }
            
            Wx = 0;
            Wy += offy;
         }

         // Visualize difference between exact and approx
         VisualizeField(vis_rho_err, vishost, visport, rho_err,
                        "Error: Density", Wx, Wy, Ww, Wh);
         
         Wx += offx;
         VisualizeField(vis_v_err, vishost, visport, vel_err,
                        "Error: Velocity", Wx, Wy, Ww, Wh);
         
         Wx += offx;
         VisualizeField(vis_ste_err, vishost, visport, ste_err,
                        "Error: Specific Total Energy", Wx, Wy, Ww, Wh);
         
         if (problem == 1 || problem == 3 || problem == 16)
         {
            Wx += offx;
            VisualizeField(vis_p_err, vishost, visport, p_err,
                           "Error: Pressure", Wx, Wy, Ww, Wh);
         }
      }
   }

   // Print initialized mesh and gridfunctions
   // Can be visualized with glvis -np # -m *.mesh
   //                        glvis -m *.mesh -g *.gf
   std::ostringstream smesh_name, mesh_name, rho_name, v_name, ste_name, mv_name;
   mesh_name << gfprint_path 
               << setfill('0') 
               << setw(6)
               << 0
               << ".mesh";
   smesh_name << gfprint_path 
               << setfill('0') 
               << setw(6)
               << 0
               << ".smesh";
   rho_name  << gfprint_path 
               << setfill('0') 
               << setw(6)
               << 0
               << "_rho.gf";
   v_name << gfprint_path 
            << setfill('0') 
            << setw(6)
            << 0
            << "_v.gf";
   ste_name << gfprint_path 
            << setfill('0') 
            << setw(6)
            << 0
            << "_ste.gf";
   mv_name << gfprint_path 
            << setfill('0') 
            << setw(6)
            << 0
            << "_mv.gf";

   std::ofstream mesh_ofs(mesh_name.str().c_str());
   mesh_ofs.precision(8);
   pmesh->PrintAsOne(mesh_ofs);
   mesh_ofs.close();

   std::ofstream smesh_ofs(smesh_name.str().c_str());
   smesh_ofs.precision(8);
   smesh->PrintAsOne(smesh_ofs);
   smesh_ofs.close();

   std::ofstream rho_ofs(rho_name.str().c_str());
   rho_ofs.precision(8);
   rho_gf.SaveAsOne(rho_ofs);
   rho_ofs.close();

   std::ofstream v_ofs(v_name.str().c_str());
   v_ofs.precision(8);
   v_gf.SaveAsOne(v_ofs);
   v_ofs.close();

   std::ofstream ste_ofs(ste_name.str().c_str());
   ste_ofs.precision(8);
   ste_gf.SaveAsOne(ste_ofs);
   ste_ofs.close();

   /* Print gamma/pressure grid function for Triple Point problem */
   if (problem == 1 || problem == 3 || problem == 16 || elastic_eos)
   {
      std::ostringstream _press_name;
      _press_name << gfprint_path 
                  << "press.gf";

      std::ofstream press_ofs(_press_name.str().c_str());
      press_ofs.precision(8);
      press_gf.SaveAsOne(press_ofs);
      press_ofs.close();
   }
   if (problem == 3)
   {
      std::ostringstream gamma_name;
      gamma_name << gfprint_path 
                  << "gamma.gf";

      std::ofstream gamma_ofs(gamma_name.str().c_str());
      gamma_ofs.precision(8);
      gamma_gf.SaveAsOne(gamma_ofs);
      gamma_ofs.close();
   }

   ParaViewDataCollection * paraview_dc;
   if (pview)
   {
      paraview_dc = new ParaViewDataCollection("ParaView", pmesh);
      paraview_dc->SetLevelsOfDetail(order_k);
      if (order_t > 0)
      {
         paraview_dc->SetHighOrderOutput(true);
      }
      paraview_dc->SetDataFormat(VTKFormat::BINARY);
      paraview_dc->SetPrefixPath(_refinement_path);
      paraview_dc->RegisterField("Specific Volume", &sv_gf);
      paraview_dc->RegisterField("Density", &rho_gf);
      paraview_dc->RegisterField("Density c", &rho_cont_gf);
      paraview_dc->RegisterField("Velocity", &v_gf);
      paraview_dc->RegisterField("Specific Total Energy", &ste_gf);
      paraview_dc->RegisterField("Mass Loss", &mc_gf);
      paraview_dc->RegisterField("Mass per Cell", &mass_gf);
      paraview_dc->RegisterField("Volume per Cell", &vol_gf);
      paraview_dc->RegisterField("Pressure", &press_gf);
      paraview_dc->RegisterField("Gamma", &gamma_gf);
      paraview_dc->RegisterField("Deviatoric Stress Frobenius Norm", &sigma_gf);
      paraview_dc->RegisterField("e shear", &e_sheer_gf);
      paraview_dc->SetCycle(0);
      paraview_dc->SetTime(0.0);
      paraview_dc->Save();
   }

   // Perform time-integration (looping over the time iterations, ti, with a
   // time-step dt). The object oper is of type LagrangianLOOperator that
   // defines the Mult() method that used by the time integrators.
   ode_solver->Init(hydro);
   double t = t_init;
   bool last_step = false;
   int steps = 0;
   BlockVector S_old(S);

   if (convergence_testing)
   {
      dt = hmin; // Set timestep to smalled h val for convergence testing
   }

   cout << "Entering time loop\n";
   chrono.Clear();
   chrono.Start();

   // Optionally suppress std::ios output to console
   ofstream file("/dev/null");

   //save cout stream buffer
   streambuf* strm_buffer = cout.rdbuf();
   if (suppress_output)
   {
      // redirect cout to /dev/null
      cout.rdbuf(file.rdbuf());
   }

   int ti = 1;
   // If at any point the mesh twists, break out of the loop and output results
   bool isCollapsed = false;
   for (; !last_step; ti++)
   {
      /* Optionally enable greedy viscosity */
      if (gmv_to_greedy_steps == ti && greedy)
      {
         hydro.SetViscOption(greedy);
      }
      hydro.chrono_dij.Start();
      hydro.BuildCijMatrices();
      hydro.BuildDijMatrix(S);
      hydro.chrono_dij.Stop();
      /* Check if we need to change CFL */
      if (problem_class->get_cfl_change() && t > problem_class->get_cfl_time_change() && hydro.GetCFL() != problem_class->get_cfl_second())
      {
         cout << "Changing CFL for Saltzman at time = " << t << endl;
         double CFL_new = problem_class->get_cfl_second();
         hydro.SetCFL(CFL_new);
      }

      if (ti == 1 && greedy)
      {
         dt = 1.E-4;
      }
      else
      {
         hydro.CalculateTimestep(S);
         dt = hydro.GetTimestep();
      }
      
      if (t + dt >= t_final)
      {
         dt = t_final - t;
         last_step = true;
      }

      if (ti == max_tsteps) { last_step = true; }

      S_old = S;
      hydro.UpdateMeshVelocityBCs(t,dt);
      ode_solver->Step(S, t, dt);
      hydro.EnforceL2BC(S, t, dt);
      steps++;

      // Make sure that the mesh corresponds to the new solution state. This is
      // needed, because some time integrators use different S-type vectors
      // and the oper object might have redirected the mesh positions to those.
      x_gf.SyncAliasMemory(S);
      pmesh->NewNodes(x_gf, false);
      
      // Optionally, post process the density
      if (post_process_density)
      {
         hydro.SetMassConservativeDensity(S);
      }

      isCollapsed = hydro.ComputeTimeSeriesData(S,t,dt);

      // MFEM_WARNING("Want adaptive time step control");

      /* End iteration if the mesh has collapsed */
      if (check_mesh && isCollapsed)
      {
         cout << "\n\n!!!!!!!!!!!!!!! Mesh has collapsed. !!!!!!!!!!!!!!!\n\n";
         last_step = true;
      }
      
      if (S_old.GetData() == S.GetData()) { cout << "\t State has not changed with step.\n"; return -1; }

      // Ensure the sub-vectors x_gf, v_gf, and e_gf know the location of the
      // data in S. This operation simply updates the Memory validity flags of
      // the sub-vectors to match those of S.
      sv_gf.SyncAliasMemory(S);
      v_gf.SyncAliasMemory(S);
      ste_gf.SyncAliasMemory(S);

      if (last_step || (ti % vis_steps) == 0)
      {
         // Turn off ios supppression to print status message
         if (suppress_output)
         {
            // restore cout stream buffer
            cout.rdbuf(strm_buffer);
         }

         double lnorm = ste_gf * ste_gf, norm;
         MPI_Allreduce(&lnorm, &norm, 1, MPI_DOUBLE, MPI_SUM, pmesh->GetComm());
 
         if (Mpi::Root())
         {
            const double sqrt_norm = sqrt(norm);

            cout << std::fixed;
            cout << "step " << std::setw(5) << ti
                 << ",\tt = " << std::setw(8) << std::setprecision(6) << t
                 << ",\tdt = " << std::setw(8) << std::setprecision(8) << dt
                 << ",\t|e| = " << std::setprecision(10) << std::scientific
                 << sqrt_norm
                 << endl;
         }
         // Fill grid function with mass information
         double mass_loss;
         hydro.ValidateMassConservation(S, mc_gf, mass_loss);

         // Turn back on suppression
         if (suppress_output)
         {
            // redirect cout to /dev/null
            cout.rdbuf(file.rdbuf());
         }

         // Compute values for visualization
         if (visualization || gfprint || pview)
         {
            /* Density */
            for (int i = 0; i < sv_gf.Size(); i++)
            {
               rho_gf[i] = 1./sv_gf[i];
            } 
            hydro.MassesAndVolumesAtPosition(rho_gf, x_gf, mass_gf, vol_gf);

            /* Continuous density projection */
            GridFunctionCoefficient rho_gf_coeff(&rho_gf);
            rho_cont_gf.ProjectDiscCoefficient(rho_gf_coeff, mfem::ParGridFunction::AvgType::ARITHMETIC);

            /* Pressure */
            hydro.ComputePressGF(S, press_gf);

            /* Elastic attributes*/
            if (elastic_eos)
            {
               // Compute Sigma and F
               hydro.ComputeSigmaGF(S, sigma_gf);
               hydro.ComputeFGF(f_gf);
               hydro.ComputeESheerGF(e_sheer_gf);

               for (int i = 0; i < L2FESpace.GetNDofs(); i++)
               {
                  frho_gf[i] = f_gf[i] * rho_gf[i];
               }
            }
         }

         /* Visualize at time t */
         if (visualization)
         {            
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

            if (problem == 1 || problem == 3 || problem == 16 || elastic_eos) // Visualize pressure
            {
               Wx += offx;

               // Pressure
               VisualizeField(vis_press, vishost, visport, press_gf,
                              "Pressure", Wx, Wy, Ww, Wh);
            }
            if (problem == 3) // Visualize gamma
            {
               vis_gamma.precision(8);
               ParGridFunction gamma_gf(&L2FESpace);
               for (int i = 0; i < gamma_gf.Size(); i++)
               {
                  int el_i = L2FESpace.GetElementForDof(i);
                  gamma_gf[i] = problem_class->get_gamma(pmesh->GetAttribute(el_i));
               }
               // gamma
               Wx += offx;
               VisualizeField(vis_gamma, vishost, visport, gamma_gf,
                              "Gamma", Wx, Wy, Ww, Wh);
            }

            if (debug_mesh_velocity)
            {
               ParGridFunction mv_gf(&H1FESpace);
               hydro.GetMV(mv_gf);
               Wx += offx;
               VisualizeField(vis_mv, vishost, visport, mv_gf,
                              "Mesh Velocity", Wx, Wy, Ww, Wh);
            }

            Wx += offx;
            VisualizeField(vis_mc, vishost, visport, mc_gf,
                           "Mass Conservation",
                           Wx, Wy, Ww, Wh);
            
            Wx = 0;
            Wy += offy;

            // MFEM_ABORT("Need to implement get_gamma with a cell_attr variable.");

            //NF//MS
            if (elastic_eos)
            {
               // Visualize
               Wx = 0;
               Wy += offy;

               VisualizeField(vis_sig, vishost, visport, sigma_gf,
                              "Sigma", Wx, Wy, Ww, Wh);
               Wx += offx;
               VisualizeField(vis_f, vishost, visport, f_gf,
                              "F", Wx, Wy, Ww, Wh);
               Wx += offx;
               VisualizeField(vis_frho, vishost, visport, frho_gf,
                              "F/rho", Wx, Wy, Ww, Wh);
               Wx += offx;
               VisualizeField(vis_esheer, vishost, visport, e_sheer_gf,
                              "e shear", Wx, Wy, Ww, Wh);
            }

            if (problem_class->has_exact_solution())
            {
               ParGridFunction rho_err(rho_gf), vel_err(v_gf), ste_err(ste_gf), p_err(press_gf);

               // Set current time for exact solution coefficient variables
               rho_coeff.SetTime(t);
               v_coeff.SetTime(t);
               ste_coeff.SetTime(t);
               p_coeff.SetTime(t);

               /* Kidder exact solution is defined in terms of initial partical locations, hence must use Finite Element Spaces defined on the undisturbed mesh */
               if (problem_class->get_indicator() == "Kidder")
               {
                  ParFiniteElementSpace L2FESpace_kidder(pmesh0, &L2FEC);
                  ParFiniteElementSpace L2VFESpace_kidder(pmesh0, &L2FEC, dim);

                  ParGridFunction rho_k_ex_gf(&L2FESpace_kidder), vel_k_ex_gf(&L2VFESpace_kidder), ste_k_ex_gf(&L2FESpace_kidder), p_k_ex_gf(&L2FESpace_kidder);
                  rho_k_ex_gf.ProjectCoefficient(rho_coeff);
                  rho_k_ex_gf[0] = rho_gf[0];
                  vel_k_ex_gf.ProjectCoefficient(v_coeff);
                  ste_k_ex_gf.ProjectCoefficient(ste_coeff);
                  p_k_ex_gf.ProjectCoefficient(p_coeff);

                  rho_err -= rho_k_ex_gf;
                  vel_err -= vel_k_ex_gf;
                  ste_err -= ste_k_ex_gf;
                  p_err -= p_k_ex_gf;

                  /* Visualizing the exact solution is not possible since the grid functions
                   representing the exact solution are defined based on the initial mesh, so
                   the solution that would be plotted would not correspond to the mesh. */
               }
               else
               {
                  if (problem_class->get_indicator() == "Vdw1")
                  {
                     problem_class->update(x_gf, t);
                  }

                  ParGridFunction rho_ex_gf(&L2FESpace), vel_ex_gf(&L2VFESpace), ste_ex_gf(&L2FESpace), p_ex_gf(&L2FESpace);
                  /* rho */
                  // l2_rho.ProjectCoefficient(rho_coeff);
                  // rho_ex_gf.ProjectGridFunction(l2_rho);
                  rho_ex_gf.ProjectCoefficient(rho_coeff);
                  // rho_ex_gf[0] = rho_gf[0];
                  /* v */
                  // l2_v.ProjectCoefficient(v_coeff);
                  // vel_ex_gf.ProjectGridFunction(l2_v);
                  vel_ex_gf.ProjectCoefficient(v_coeff);
                  /* ste */
                  // l2_ste.ProjectCoefficient(ste_coeff);
                  // ste_ex_gf.ProjectGridFunction(l2_ste);
                  ste_ex_gf.ProjectCoefficient(ste_coeff);
                  /* pressure */
                  // ParGridFunction l2_p(&l2_fes);
                  // l2_p.ProjectCoefficient(p_coeff);
                  // p_ex_gf.ProjectGridFunction(l2_p);
                  p_ex_gf.ProjectCoefficient(p_coeff);

                  rho_err -= rho_ex_gf;
                  vel_err -= vel_ex_gf;
                  ste_err -= ste_ex_gf;
                  p_err -= p_ex_gf;

                  // Visualize exact solution
                  VisualizeField(vis_rho_ex, vishost, visport, rho_ex_gf,
                                 "Exact: Density", Wx, Wy, Ww, Wh);
                  
                  Wx += offx;
                  VisualizeField(vis_v_ex, vishost, visport, vel_ex_gf,
                                 "Exact: Velocity", Wx, Wy, Ww, Wh);
                  Wx += offx;
                  VisualizeField(vis_ste_ex, vishost, visport, ste_ex_gf,
                                 "Exact: Specific Total Energy", Wx, Wy, Ww, Wh);
                  if (problem == 1 || problem == 3)
                  {
                     Wx += offx;
                     VisualizeField(vis_p_ex, vishost, visport, p_ex_gf,
                                    "Exact: Pressure", Wx, Wy, Ww, Wh);
                  }
                  
                  Wx = 0;
                  Wy += offy;
               }

               // Visualize difference between exact and approx
               VisualizeField(vis_rho_err, vishost, visport, rho_err,
                              "Error: Density", Wx, Wy, Ww, Wh);
               Wx += offx;
               VisualizeField(vis_v_err, vishost, visport, vel_err,
                              "Error: Velocity", Wx, Wy, Ww, Wh);
               Wx += offx;
               VisualizeField(vis_ste_err, vishost, visport, ste_err,
                              "Error: Specific Total Energy", Wx, Wy, Ww, Wh);
               
               if (problem == 1 || problem == 3 || problem == 16)
               {
                  Wx += offx;
                  VisualizeField(vis_p_err, vishost, visport, p_err,
                                 "Error: Pressure", Wx, Wy, Ww, Wh);
               }
               
            }
         } /* END VIS */

         if (gfprint)
         {
            // Save mesh and gfs to files
            std::ostringstream mesh_name, rho_name, v_name, ste_name, press_name, mass_name, mv_name;
            mesh_name << gfprint_path 
                      << setfill('0') 
                      << setw(6)
                      << ti 
                      << ".mesh";
            rho_name  << gfprint_path 
                      << setfill('0') 
                      << setw(6)
                      << ti 
                      << "_rho.gf";
            v_name << gfprint_path 
                   << setfill('0') 
                   << setw(6)
                   << ti 
                   << "_v.gf";
            ste_name << gfprint_path 
                     << setfill('0') 
                     << setw(6)
                     << ti  
                     << "_ste.gf";
            press_name << gfprint_path 
                       << setfill('0') 
                       << setw(6)
                       << ti  
                       << "_press.gf";
            mass_name << gfprint_path 
                      << setfill('0') 
                      << setw(6)
                      << ti  
                      << "_mass_loss.gf";
            mv_name << gfprint_path 
                    << setfill('0') 
                    << setw(6)
                    << ti 
                    << "_mv.gf";

            std::ofstream mesh_ofs(mesh_name.str().c_str());
            mesh_ofs.precision(8);
            pmesh->PrintAsOne(mesh_ofs);
            mesh_ofs.close();

            std::ofstream rho_ofs(rho_name.str().c_str());
            rho_ofs.precision(8);
            rho_gf.SaveAsOne(rho_ofs);
            rho_ofs.close();

            std::ofstream v_ofs(v_name.str().c_str());
            v_ofs.precision(8);
            v_gf.SaveAsOne(v_ofs);
            v_ofs.close();

            std::ofstream ste_ofs(ste_name.str().c_str());
            ste_ofs.precision(8);
            ste_gf.SaveAsOne(ste_ofs);
            ste_ofs.close();

            std::ofstream press_ofs(press_name.str().c_str());
            press_ofs.precision(8);
            press_gf.SaveAsOne(press_ofs);
            press_ofs.close();

            // Print mass loss
            std::ofstream mass_ofs(mass_name.str().c_str());
            mass_ofs.precision(8);
            mc_gf.SaveAsOne(mass_ofs);
            mass_ofs.close();

            // Print continuous interpolation of density
            std::ostringstream rho_cont_name;
            rho_cont_name  << gfprint_path 
                           << setfill('0') 
                           << setw(6)
                           << ti 
                           << "_rho_c.gf";
            std::ofstream rho_cont_ofs(rho_cont_name.str().c_str());
            rho_cont_ofs.precision(8);
            rho_cont_gf.SaveAsOne(rho_cont_ofs);
            rho_cont_ofs.close();       
         } // gfprint
         if (pview)
         {
            paraview_dc->SetCycle(ti);
            paraview_dc->SetTime(t);
            paraview_dc->Save();
         } // pview  
      }
   } // End time step iteration
   chrono.Stop();

   switch (ode_solver_type)
   {
      case 12:
      case 2: steps *= 2; break;
      case 13:
      case 3: steps *= 3; break;
      case 14:
      case 4: steps *= 4; break;
      case 16:
      case 6: steps *= 6; break;
      case 7: steps *= 2;
   }

   // Last time computing density
   for (int i = 0; i < sv_gf.Size(); i++)
   {
      rho_gf[i] = 1./sv_gf[i];
   } 

   // Print grid functions to files
   ostringstream sv_filename_suffix;
   sv_filename_suffix << "sv_"
                      << setfill('0')
                      << setw(2)
                      << to_string(rp_levels + rs_levels)
                      << ".csv";

   ostringstream ts_filename_suffix;
   ts_filename_suffix << "tsData_"
                      << setfill('0')
                      << setw(2)
                      << to_string(rp_levels + rs_levels)
                      << ".csv";

   hydro.SaveStateVecsToFile(S, sv_output_prefix, sv_filename_suffix.str());
   hydro.SaveTimeSeriesArraysToFile(ts_output_prefix, ts_filename_suffix.str());

   // In all cases, print the final grid functions
   {
      // Save mesh and gfs to files
      std::ostringstream mesh_name, smesh_name, rho_name, v_name, ste_name, massC_name, mv_name;
      mesh_name << gfprint_path 
                << "final.mesh";
      smesh_name << gfprint_path 
                 << "final.smesh";
      rho_name  << gfprint_path 
                << "rho_final.gf";
      v_name << gfprint_path 
             << "v_final.gf";
      ste_name << gfprint_path 
               << "ste_final.gf";
      massC_name << gfprint_path
                 << "mass_loss.gf";
      mv_name << gfprint_path 
             << "mv_final.gf";

      std::ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      pmesh->PrintAsOne(mesh_ofs);
      mesh_ofs.close();

      std::ofstream smesh_ofs(smesh_name.str().c_str());
      smesh_ofs.precision(8);
      smesh->PrintAsOne(smesh_ofs);
      smesh_ofs.close();

      std::ofstream rho_ofs(rho_name.str().c_str());
      rho_ofs.precision(8);
      rho_gf.SaveAsOne(rho_ofs);
      rho_ofs.close();

      std::ofstream v_ofs(v_name.str().c_str());
      v_ofs.precision(8);
      v_gf.SaveAsOne(v_ofs);
      v_ofs.close();

      std::ofstream ste_ofs(ste_name.str().c_str());
      ste_ofs.precision(8);
      ste_gf.SaveAsOne(ste_ofs);
      ste_ofs.close();

      std::ofstream mc_ofs(massC_name.str().c_str());
      mc_ofs.precision(12); // Need mass loss to machine error precision
      mc_gf.SaveAsOne(mc_ofs);
      mc_ofs.close();

      // Print continuous interpolation of density
      GridFunctionCoefficient rho_gf_coeff(&rho_gf);
      rho_cont_gf.ProjectDiscCoefficient(rho_gf_coeff, mfem::ParGridFunction::AvgType::ARITHMETIC);
      std::ostringstream rho_cont_name;
      rho_cont_name  << gfprint_path 
                     << "rho_c_final.gf";
      std::ofstream rho_cont_ofs(rho_cont_name.str().c_str());
      rho_cont_ofs.precision(8);
      rho_cont_gf.SaveAsOne(rho_cont_ofs);
      rho_cont_ofs.close();
   } // End print final grid functions

   if (problem_class->has_exact_solution())
   {
      // Save exact to file as well
      BlockVector S_exact(offset, Device::GetMemoryType());
      ParGridFunction x_gf_exact, sv_ex_gf, vel_ex_gf, ste_ex_gf;

      x_gf_exact.MakeRef(&H1FESpace, S_exact, offset[0]);
      sv_ex_gf.MakeRef(&L2FESpace, S_exact, offset[1]);
      vel_ex_gf.MakeRef(&L2VFESpace, S_exact, offset[2]);
      ste_ex_gf.MakeRef(&L2FESpace, S_exact, offset[3]);

      // Project exact solution
      if (problem_class->get_indicator() == "Vdw1")
      {
         problem_class->update(x_gf, t);
      }

      sv_coeff.SetTime(t);
      v_coeff.SetTime(t);
      ste_coeff.SetTime(t);

      /* sv */
      // l2_sv.ProjectCoefficient(sv_coeff);
      // sv_ex_gf.ProjectGridFunction(l2_sv);
      sv_ex_gf.ProjectCoefficient(sv_coeff);
      /* v */
      // l2_v.ProjectCoefficient(v_coeff);
      // vel_ex_gf.ProjectGridFunction(l2_v);
      vel_ex_gf.ProjectCoefficient(v_coeff);
      /* ste */
      // l2_ste.ProjectCoefficient(ste_coeff);
      // ste_ex_gf.ProjectGridFunction(l2_ste);
      ste_ex_gf.ProjectCoefficient(ste_coeff);

      // Print grid functions to files
      ostringstream sv_ex_filename_suffix;
      sv_ex_filename_suffix << setfill('0') << setw(2)
                        << "sv_exact"
                        << ".out";
      
      hydro.SaveStateVecsToFile(S_exact, sv_output_prefix, sv_ex_filename_suffix.str());
   }

   /* 
   * Whether or not the exact solution is known, we are still interested
   * in the mass defect convergence.  In the case an exact solution is known,
   * print the cooresponding L1, L2, Linf errors as well.
   */
   ostringstream convergence_filename;
   convergence_filename << output_path << "convergence/temp_output/np" << num_procs;

   /* Coefficient to assist in computation of errors */
   ConstantCoefficient zero(0.0);
   
   /* Values to store numerators, to be computed on case by case basis since exact solutions vary */
   double sv_L1_error_n = 0., rho_L1_error_n = 0., vel_L1_error_n = 0., ste_L1_error_n = 0.,
          sv_L2_error_n = 0., rho_L2_error_n = 0., vel_L2_error_n = 0., ste_L2_error_n = 0.,
          sv_Max_error_n = 0., rho_Max_error_n = 0., vel_Max_error_n = 0., ste_Max_error_n = 0.;

   if (problem_class->has_exact_solution())
   {
      /* Set coefficients to final time */
      rho_coeff.SetTime(t);
      v_coeff.SetTime(t);
      ste_coeff.SetTime(t);
      sv_coeff.SetTime(t);

      if (problem_class->get_indicator() == "Kidder")
      {
         ParFiniteElementSpace L2FESpace_kidder(pmesh0, &L2FEC);
         ParFiniteElementSpace L2VFESpace_kidder(pmesh0, &L2FEC, dim);

         ParGridFunction sv_k_ex_gf(&L2FESpace_kidder), vel_k_ex_gf(&L2VFESpace_kidder), ste_k_ex_gf(&L2FESpace_kidder);
         sv_k_ex_gf.ProjectCoefficient(sv_coeff);
         vel_k_ex_gf.ProjectCoefficient(v_coeff);
         ste_k_ex_gf.ProjectCoefficient(ste_coeff);

         /* Compute relative errors */
         GridFunctionCoefficient vel_ex_coeff(&vel_k_ex_gf), ste_ex_coeff(&ste_k_ex_gf), sv_ex_coeff(&sv_k_ex_gf);

         // Velocity errors
         vel_L1_error_n = v_gf.ComputeL1Error(vel_ex_coeff) / vel_k_ex_gf.ComputeL1Error(zero);
         vel_L2_error_n = v_gf.ComputeL2Error(vel_ex_coeff) / vel_k_ex_gf.ComputeL2Error(zero);
         vel_Max_error_n = v_gf.ComputeMaxError(vel_ex_coeff) / vel_k_ex_gf.ComputeMaxError(zero);

         // Specific volume and specific total energy
         sv_L1_error_n = sv_gf.ComputeL1Error(sv_ex_coeff) / sv_k_ex_gf.ComputeL1Error(zero);
         sv_L2_error_n = sv_gf.ComputeL2Error(sv_ex_coeff) / sv_k_ex_gf.ComputeL2Error(zero);
         sv_Max_error_n = sv_gf.ComputeMaxError(sv_ex_coeff) / sv_k_ex_gf.ComputeMaxError(zero);

         ste_L1_error_n = ste_gf.ComputeL1Error(ste_ex_coeff) / ste_k_ex_gf.ComputeL1Error(zero);
         ste_L2_error_n = ste_gf.ComputeL2Error(ste_ex_coeff) / ste_k_ex_gf.ComputeL2Error(zero);
         ste_Max_error_n = ste_gf.ComputeMaxError(ste_ex_coeff) / ste_k_ex_gf.ComputeMaxError(zero);
      }
      else
      {
         if (problem_class->get_indicator() == "Vdw1")
         {
            problem_class->update(x_gf, t);
         }
         // Compute errors
         ParGridFunction rho_ex_gf(&L2FESpace), vel_ex_gf(&L2VFESpace), ste_ex_gf(&L2FESpace), sv_ex_gf(&L2FESpace);
         /* sv */
         // l2_sv.ProjectCoefficient(sv_coeff);
         // sv_ex_gf.ProjectGridFunction(l2_sv);
         sv_ex_gf.ProjectCoefficient(sv_coeff);
         /* rho */
         // l2_rho.ProjectCoefficient(rho_coeff);
         // rho_ex_gf.ProjectGridFunction(l2_rho);
         rho_ex_gf.ProjectCoefficient(rho_coeff);
         /* v */
         // l2_v.ProjectCoefficient(v_coeff);
         // vel_ex_gf.ProjectGridFunction(l2_v);
         vel_ex_gf.ProjectCoefficient(v_coeff);
         /* ste */
         // l2_ste.ProjectCoefficient(ste_coeff);
         // ste_ex_gf.ProjectGridFunction(l2_ste);
         ste_ex_gf.ProjectCoefficient(ste_coeff);

         // In the case of the Noh Problem, project 0 on the boundary of approx and exact
         if ((problem_class->get_indicator() == "ElasticNoh" || problem_class->get_indicator() == "Noh"))
         {
            cout << "[Noh] Projecting zero on the boundary cells.\n";
            ParGridFunction cell_bdr_flag_gf;
            hydro.GetCellBdrFlagGF(cell_bdr_flag_gf);

            for (int i = 0; i < pmesh->GetNE(); i++)
            {
               if (cell_bdr_flag_gf[i] != -1)
               {
                  // We have a boundary cell
                  rho_gf[i] = 0.;
                  ste_gf[i] = 0.;
                  rho_ex_gf[i] = 0.;
                  ste_ex_gf[i] = 0.;
                  sv_ex_gf[i] = 0.;
                  for (int j = 0; j < dim; j++)
                  {
                     int index = i + j*pmesh->GetNE();
                     v_gf[index] = 0.;
                     vel_ex_gf[index] = 0.;
                  }
               }
            }
         }

         /* Project 0 on all extrapolated cells, marked with attr = 99 */
         if (pmesh->attributes.Find(99) != -1)
         {
            cout << "projecting zero on cells with attr 99\n";
            Vector _vec_zero(dim);
            _vec_zero = 0.;
            VectorConstantCoefficient _zero_vcc(_vec_zero);
            // onto approx
            sv_gf.ProjectCoefficient(_zero_vcc, 99);
            rho_gf.ProjectCoefficient(_zero_vcc, 99);
            v_gf.ProjectCoefficient(_zero_vcc, 99);
            ste_gf.ProjectCoefficient(_zero_vcc, 99);

            // onto exact
            sv_ex_gf.ProjectCoefficient(_zero_vcc, 99);
            rho_ex_gf.ProjectCoefficient(_zero_vcc, 99);
            vel_ex_gf.ProjectCoefficient(_zero_vcc, 99);
            ste_ex_gf.ProjectCoefficient(_zero_vcc, 99);
         }

         /* Compute relative errors */
         GridFunctionCoefficient rho_ex_coeff(&rho_ex_gf), vel_ex_coeff(&vel_ex_gf), ste_ex_coeff(&ste_ex_gf), sv_ex_coeff(&sv_ex_gf);

         // Velocity errors
         vel_L1_error_n = v_gf.ComputeL1Error(vel_ex_coeff) / vel_ex_gf.ComputeL1Error(zero);
         vel_L2_error_n = v_gf.ComputeL2Error(vel_ex_coeff) / vel_ex_gf.ComputeL2Error(zero);
         vel_Max_error_n = v_gf.ComputeMaxError(vel_ex_coeff) / vel_ex_gf.ComputeMaxError(zero);

         // Specific volume and specific total energy
         /* In the Sedov case, we do not get convergence of the specific total energy, but the total energy */
         if (problem == 1)
         {
            ParGridFunction te_ex_gf(&L2FESpace), te_gf(&L2FESpace);
            for (int i = 0; i < ste_ex_gf.Size(); i++)
            {
               te_gf[i] = rho_gf[i] * ste_gf[i];
               te_ex_gf[i] = rho_ex_gf[i] * ste_ex_gf[i];
            }
            sv_L1_error_n = rho_gf.ComputeL1Error(rho_ex_coeff) / rho_ex_gf.ComputeL1Error(zero);
            sv_L2_error_n = rho_gf.ComputeL2Error(rho_ex_coeff) / rho_ex_gf.ComputeL2Error(zero);
            sv_Max_error_n = rho_gf.ComputeMaxError(rho_ex_coeff) / rho_ex_gf.ComputeMaxError(zero);

            GridFunctionCoefficient te_ex_coeff(&te_ex_gf);
            ste_L1_error_n = te_gf.ComputeL1Error(te_ex_coeff) / te_ex_gf.ComputeL1Error(zero);
            ste_L2_error_n = te_gf.ComputeL2Error(te_ex_coeff) / te_ex_gf.ComputeL2Error(zero);
            ste_Max_error_n = te_gf.ComputeMaxError(te_ex_coeff) / te_ex_gf.ComputeMaxError(zero);
         }
         else 
         {
            sv_L1_error_n = sv_gf.ComputeL1Error(sv_ex_coeff) / sv_ex_gf.ComputeL1Error(zero);
            sv_L2_error_n = sv_gf.ComputeL2Error(sv_ex_coeff) / sv_ex_gf.ComputeL2Error(zero);
            sv_Max_error_n = sv_gf.ComputeMaxError(sv_ex_coeff) / sv_ex_gf.ComputeMaxError(zero);

            ste_L1_error_n = ste_gf.ComputeL1Error(ste_ex_coeff) / ste_ex_gf.ComputeL1Error(zero);
            ste_L2_error_n = ste_gf.ComputeL2Error(ste_ex_coeff) / ste_ex_gf.ComputeL2Error(zero);
            ste_Max_error_n = ste_gf.ComputeMaxError(ste_ex_coeff) / ste_ex_gf.ComputeMaxError(zero);
         }

      } // End all other test problems
      
   } // End if problem has exact solution

   /* Get composite errors values, will return 0 if exact solution is not known */
   const double L1_error = (sv_L1_error_n + vel_L1_error_n + ste_L1_error_n) / 3.;
   const double L2_error = (sv_L2_error_n + vel_L2_error_n + ste_L2_error_n) / 3.;
   const double Max_error = (sv_Max_error_n + vel_Max_error_n + ste_Max_error_n) / 3.;

   /* Calculate mass loss */
   double mass_loss;
   hydro.ValidateMassConservation(S, mc_gf, mass_loss);

   /* In either case, write convergence file. */
   if (Mpi::Root())
   {
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
                        // sv
                        << "sv_L1_Error " << sv_L1_error_n << "\n"
                        << "sv_L2_Error " << sv_L2_error_n << "\n"
                        << "sv_Linf_Error " << sv_Max_error_n << "\n"
                        // vel
                        << "vel_L1_Error " << vel_L1_error_n << "\n"
                        << "vel_L2_Error " << vel_L2_error_n << "\n"
                        << "vel_Linf_Error " << vel_Max_error_n << "\n"
                        // ste
                        << "ste_L1_Error " << ste_L1_error_n << "\n"
                        << "ste_L2_Error " << ste_L2_error_n << "\n"
                        << "ste_Linf_Error " << ste_Max_error_n << "\n"
                        // total
                        << "L1_Error " << L1_error << "\n"
                        << "L2_Error " << L2_error << "\n"
                        << "Linf_Error " << Max_error << "\n"
                        << "mass_loss " << mass_loss << "\n"
                        << "dt " << dt << "\n"
                        << "Endtime " << t << "\n";
                  
      convergence_file.close();
   } // Writing convergence file
   
   if (suppress_output)
   {
      // restore cout stream buffer
      cout.rdbuf(strm_buffer);
   }

   /* Close visualization object */
   if (visualization)
   {
      // Make sure all MPI ranks have sent their 'v' solution before initiating
      // another set of GLVis connections (one from each rank):
      MPI_Barrier(pmesh->GetComm());

      vis_rho.close();
      vis_v.close();
      vis_ste.close();

      vis_rho_ex.close();
      vis_v_ex.close();
      vis_ste_ex.close();
      vis_p_ex.close();

      vis_rho_err.close();
      vis_v_err.close();
      vis_ste_err.close();
      vis_p_err.close();

      vis_press.close();
      vis_gamma.close();
      vis_mc.close();

      vis_sig.close();
      vis_f.close();
      vis_frho.close();
      vis_esheer.close();
   }

   cout << "Program took " << chrono.RealTime() << "s.\n";
   double temp_denom = ti * L2FESpace.GetNE();
   double time_per_gp_ts = chrono.RealTime() / temp_denom;
   cout << "This amounts to " << time_per_gp_ts << " s per timestep per gridpoint.\n";

   cout << "==================== Runtime breakdown ====================\n"
        << setw(59) <<  " total time | t / (dof*ts) \n"
        << setprecision(3)
        << "  Total time: " << setw(28) << chrono.RealTime() << " | " << chrono.RealTime() / temp_denom << endl
        << "  Dij calculation: " << setw(23) << hydro.chrono_dij.RealTime() << " | " << hydro.chrono_dij.RealTime() / temp_denom << endl
        << "  Mesh Motion calculation: " << setw(15) << hydro.chrono_mm.RealTime() << " | " << hydro.chrono_mm.RealTime() / temp_denom << endl
        << "    Hiop calculation (if used): " << setw(10) << hydro.chrono_hiop.RealTime() << " | " << hydro.chrono_hiop.RealTime() / temp_denom << endl
        << "  Mesh Motion Linearization: " << setw(13) << hydro.chrono_mm_lin.RealTime() << " | " << hydro.chrono_mm_lin.RealTime() / temp_denom << endl
        << "  State Update calculation: " << setw(14) << hydro.chrono_state.RealTime() << " | " << hydro.chrono_state.RealTime() / temp_denom << endl
        << "===========================================================\n";
   
   delete pmesh;
   delete smesh;
   delete pmesh0;
   delete CRFEC;
   delete problem_class;
   delete m;
   delete paraview_dc;

   return 0;
}

static void display_banner(std::ostream &os)
{
   os << endl
      << "       __                __                 " << endl
      << "      / /   ____  ____  / /____  _____   " << endl
      << "     / /   / __ `/ __ `/ / __ \\/ ___/ " << endl
      << "    / /___/ /_/ / /_/ / / /_/ (__  )    " << endl
      << "   /_____/\\__,_/\\__,_/_/\\____/____/  " << endl
      << "               /____/                       " << endl << endl;
}
