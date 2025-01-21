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
* Example run time parameters:
*
* ----- 1D -----
* ./Laglos -m data/ref-segment.mesh -p 0 -tf 0.6 -cfl 0.5 -rs 8            ## Smooth
* ./Laglos -m data/ref-segment.mesh -p 1 -tf 0.225 -cfl 2. -rs 8           ## Sod
* ./Laglos -m data/ref-segment.mesh -p 2 -tf 0.15 -cfl 0.5 -rs 8           ## Lax
# ./Laglos -m data/ref-segment.mesh -p 3 -tf 0.667 -cfl 0.2 -rs 8          ## Leblanc
*
* ----- vdw -----
* ./Laglos -m data/ref-segment-c0.mesh -p 8 -cfl 0.5 -tf 0.5 -rs 8 -vis    ## Vdw1
* ./Laglos -m data/segment-nhalf-1.mesh -p 9 -cfl 0.5 -tf 1.25 -rs 8 -vis  ## Vdw2 
* ./Laglos -m data/segment-nhalf-1.mesh -p 10 -cfl 0.5 -tf 0.4 -rs 8 -vis  ## Vdw3 
* ./Laglos -m data/segment-n1p7-1.mesh -p 11 -cfl 1.3 -tf 0.005 -rs 8 -vis ## Vdw4
*
* --------------------------------------------------------------
* ----- 2D -----
* --------------------------------------------------------------
* ./Laglos -m data/ref-square.mesh -p 0 -tf 0.6 -cfl 0.5 -rs 4             ## Smooth wave in 2D
* ./Laglos -m data/ref-square.mesh -p 1 -tf 0.225 -cfl 0.5 -rs 4           ## Sod in 2D
* ./Laglos -m data/distorted-square.mesh -p 1 -tf 0.225 -cfl 0.5 -rs 4     ## Sod Distorted
* ./Laglos -m data/ref-square.mesh -p 13 -tf 0.2 -cfl 0.25 -rs 4           ## Sod Radial
* ./Laglos -m data/square-vortex.mesh -p 5 -tf 10 -cfl 0.5 -rs 3           ## Isentropic Vortex
* ./Laglos -m data/noh-nonuniform.mesh -p 4 -tf 0.6 -cfl 1 -rs 0           ## Noh nonuniform (likely need to change BCs)
* ./Laglos -m data/ref-square-c0.mesh -p 4 -tf 0.6 -cfl 0.25 -rs 6         ## Noh
* ./Laglos -m data/ref-square.mesh -p 6 -tf .9 -cfl 0.1 -rs 5              ## Sedov
* ./Laglos -m data/rectangle_saltzmann.mesh -p 7 -tf 0.6 -cfl 0.01 -rs 3   ## Saltzman problem
* ----- Fails -----
* ./Laglos -m data/triple-point.mesh -p 12 -tf 5. -cfl 0.5 -rs 2           ## Triple Point [Fails]
* ----- Untested -----
* ./Laglos -m data/ref-square.mesh -p 14 -tf 0.6 -cfl 0.5 -rs 4            ## ICF [Untested]
* ./Laglos -m data/full_ring_r0.mesh -p 15 -tf 0.1887 -cfl 0.5 -rs 0       ## Kidder [Untested]
*
* ----- vdw -----
* ./Laglos -m data/tube-np5-1.mesh -p 9 -cfl 0.5 -tf 1.25 -rs 2 -vis        ## Vdw2 
* ./Laglos -m data/tube-np5-1.mesh -p 10 -cfl 0.5 -tf 0.4 -rs 2 -vis        ## Vdw3 
* ./Laglos -m data/tube-n1p7-1.mesh -p 11 -cfl 1.3 -tf 0.005 -rs 8 -vis     ## Vdw4
* 
* --- General Riemann Problem, change riemann_problem.h ---
* ./Laglos -m data/ref-segment-c0.mesh -p 20 -cfl 0.5 -tf 1 -rs 8 -vis     ## General Riemann Problem
*
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
* ./Laglos -m data/ref-square.mesh -p 1 -tf 0.225 -cfl 0.5 -rs 4 -mv 3
*
* CAVEAT Weight LS method for corner vertex calculation 
* Sod
* ./Laglos -m data/ref-square.mesh -p 1 -tf 0.225 -cfl 0.5 -rs 4 -mv 4 
*
***/
#include "mfem.hpp"
#include "var-config.h"
#include "compile_time_vals.h"
#include "laglos_solver.hpp"
#include "test_problems_include.h"
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
using namespace hydrodynamics;
using namespace std;


int main(int argc, char *argv[]) {
   StopWatch chrono;

   // Initialize MPI.
   Mpi::Init();
   const int num_procs = Mpi::WorldSize();
   const int myid = Mpi::WorldRank();
   // Hypre::Init();

   const int dim = CompileTimeVals::dim;
   const string results_dir = CompileTimeVals::results_dir;

   // Parse command line options
   const char *mesh_file_location = "default";
   const char *output_flag = "default";
   int problem = 0;
   int rs_levels = 0;
   int rp_levels = 0;
   int order_mv = 2;  // Order of mesh movement approximation space
   int order_u = 0;
   double t_init = 0.0;
   double t_final = 0.6;
   int max_tsteps = -1;
   double dt = 0.001;
   bool visualization = false;
   int vis_steps = 5;
   bool pview = false;
   bool visit = false;
   bool gfprint = false;
   int precision = 12;
   /* This parameter describes how often to compute mass corrective face velocities, 0 indicates to always take the average. */
   bool mc = false;
   bool use_viscosity = true;
   bool greedy = false;
   bool mm = true;
   bool check_mesh = true;
   bool post_process_density = false;
   int mv_option = 2;
   double mv_target_visc_coeff = 0.;
   bool do_mv_linearization = false;
   int fv_option = 2;
   int mv_it_option = 0;
   int mv_n_iterations = 0;
   double mm_visc_face = 0., mm_cell = 0.;
   bool optimize_timestep = true;
   bool convergence_testing = false;
   bool suppress_output = false;
   double CFL = 0.5;
   double dm_val = 0.; // Parameter to distort the mesh, overrides file val
   string sv_output_prefix;
   string ts_output_prefix;
   string output_path;

   OptionsParser args(argc, argv);

   args.AddOption(&problem, "-p", "--problem", "Test problem to run.");
   args.AddOption(&mesh_file_location, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&output_flag, "-of", "--output-flag", 
                  "Directory that output files should be placed");
   args.AddOption(&rs_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&rp_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
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
   args.AddOption(&visit, "-visit", "--visit", "-no-visit", "--no-visit",
                  "Enable or disable VisIt visualization.");
   args.AddOption(&gfprint, "-print", "--print", "-no-print", "--no-print",
                  "Enable or disable result output (files in mfem format).");
   args.AddOption(&mc, "-mc", "--mass-correct", "-no-mc", "--no-mass-correct",
                  "Enable or disable mass correction.");
   args.AddOption(&use_viscosity, "-visc", "--use-viscosity", "-no-visc",
                  "--no-viscosity",
                  "Enable or disable the use of artificial viscosity.");
   args.AddOption(&greedy, "-greedy", "--use-greedy-viscosity", "-no-greedy", "--no-greedy-viscosity",
                  "Enable or disable greedy viscosity computation.");
   args.AddOption(&mm, "-mm", "--move-mesh", "-no-mm", "--no-move-mesh",
                  "Enable or disable mesh movement.");
   args.AddOption(&check_mesh, "-cm", "--check-mesh", "-no-cm", "--no-check-mesh",
                  "Enable or disable checking if the mesh has twisted.");
   args.AddOption(&post_process_density, "-ppd", "--post-process-density", "-no-ppd", "--no-post-process-density",
                  "Enable or disable density post processing to guarantee conservation of mass.");
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
   args.AddOption(&optimize_timestep, "-ot", "--optimize-timestep", "-no-ot",
                  "--no-optimize-timestep",
                  "Enable or disable timestep optimization using CFL.");
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

   // Set output_flag string
   if (strncmp(output_flag, "default", 7) != 0)
   {
      output_path = results_dir + std::string(output_flag) + "/";
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
   string _visit = _refinement_path + "VisIt/";
   string _pview = _refinement_path + "ParaView/";
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

   const char* _visit_basename = _visit.c_str();
   if (visit)
   {
      boost::filesystem::path visit_output_dir(_visit_basename);
      boost::filesystem::create_directory(visit_output_dir);
   }

   const char* _pview_basename = _pview.c_str();
   if (pview)
   {
      boost::filesystem::path pview_output_dir(_pview_basename);
      boost::filesystem::create_directory(pview_output_dir);
   }

   // On all processors, use the default builtin 1D/2D/3D mesh or read the
   // serial one given on the command line.
   Mesh *mesh;
   if (strncmp(mesh_file_location, "default", 7) != 0)
   {
      std::string result = std::string(LAGLOS_DIR) + std::string(mesh_file_location);
      const char* mesh_file = result.c_str();
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
   ProblemBase<dim> * problem_class = NULL;
   switch (problem)
   {
      case 0: // Smooth
      {
         problem_class = new SmoothWave<dim>();
         break;
      }
      case 1: // Sod
      {
         problem_class = new SodProblem<dim>();
         break;
      }
      case 2: // Lax
      {
         problem_class = new LaxProblem<dim>();
         break;
      }
      case 3: // Leblanc
      {
         problem_class = new LeblancProblem<dim>();
         break;
      }
      case 4: // Noh
      {
         problem_class = new NohProblem<dim>();
         break;
      }
      case 5: // Isentropic Vortex, stationary center
      {
         problem_class = new IsentropicVortex<dim>();
         break;
      }
      case 6: // Sedov
      {
         problem_class = new SedovLLNLProblem<dim>();
         break;
      }
      case 7: // Saltzmann
      {
         problem_class = new SaltzmannProblem<dim>();
         break;
      }
      /* VDW */
      case 8:
      {
         problem_class = new VdwTest1<dim>();
         break;
      }
      case 9:
      {
         problem_class = new VdwTest2<dim>();
         break;
      }
      case 10:
      {
         problem_class = new VdwTest3<dim>();
         break;
      }
      case 11:
      {
         problem_class = new VdwTest4<dim>();
         break;
      }
      case 12: // Triple Point
      {
         problem_class = new TriplePoint<dim>();
         break;
      }
      case 13: // Radial Sod
      {
         problem_class = new SodRadial<dim>();
         break;
      }
      case 14: // ICF
      {
         problem_class = new ICFProblem<dim>();
         break;
      }
      case 15: // Kidder
      {
         problem_class = new KidderProblem<dim>();
         break;
      }
      case 20: // Riemann Problem
      {
         problem_class = new RiemannProblem<dim>();
         break;
      }
      case 21: // Sedov
      {
         MFEM_ABORT("Not implemented\n");
         assert(hmin == hmax);
         Vector params(2);
         params[0] = hmax, params[1] = pmesh->GetElementVolume(0);

         problem_class = new SedovProblem<dim>();
         problem_class->update(params, t_init);
         // TODO: Will need to modify initialization of internal energy
         //       if distorted meshes are used.
         break;
      }
      case 100:
      {
         problem_class = new TestBCs<dim>();
         break;
      }
      default:
      {
         MFEM_ABORT("Failed to initiate a problem.\n");
      }
   }

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

   /* Distort Mesh for Saltzman Problem */
   if (problem_class->get_distort_mesh() || dm_val != 0.)
   {
      cout << "distort mesh.\n";
      switch (problem)
      {
      case 1: // Sod
      {
         // Random distortion based on min mesh size
         Array<double> coords(dim);
         for (int vertex = 0; vertex < H1FESpace.GetNDofs(); vertex++)
         {
            int index = vertex + H1FESpace.GetNDofs();
            coords[0] = x_gf[vertex];
            coords[1] = x_gf[index];

            double dy = 0.;
            if (coords[1] >= 0.75 ) {
               dy = -4 * (coords[1] - 1);
            } else if (coords[1] >= 0.5) {
               dy = 4. * (coords[1] - 0.5);
               dy *= -1;
            } else if (coords[1] >= 0.25) {
               dy = -4. * (coords[1] - 0.5);
            } else if (coords[1] >= 0.) {
               dy = 4. * coords[1];
               dy *= -1.;
            }
            dy *= hmin;

            double y_new = coords[1] + 2 * (1.-coords[0]) * coords[0] * dm_val * dy;

            double x_new = coords[0] + dm_val * hmin * (1. - coords[1]) * sin(M_PI * coords[0]);
            x_gf[index] = y_new;
            x_gf[vertex] = x_new;
         }
         break;
      }
      case 7: // Saltzmann
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

      // Adjust all face nodes to be averages of their adjacent corners
      mfem::Mesh::FaceInformation FI;
      for (int face = 0; face < pmesh->GetNumFaces(); face++)
      {
         FI = pmesh->GetFaceInformation(face);
         Array<int> face_dof_row;
         H1FESpace.GetFaceDofs(face, face_dof_row);

         int face_dof = face_dof_row[2];
         int index0 = face_dof_row[0];
         int index1 = face_dof_row[1];
         for (int i = 0; i < dim; i++)
         {
            int face_dof_index = face_dof + i * H1FESpace.GetNDofs();
            int node0_index = index0 + i * H1FESpace.GetNDofs();
            int node1_index = index1 + i * H1FESpace.GetNDofs();
            x_gf[face_dof_index] = 0.5 * (x_gf[node0_index] + x_gf[node1_index]);
         }
      }

      // Adjust all cell center to be averages of corners
      Vector cell_x(dim);
      Array<int> verts;
      for (int ci = 0; ci < L2FESpace.GetNDofs(); ci++)
      {
         int cell_vdof;
         cell_x = 0.;
         // Get center node dof, average, and update
         switch (dim)
         {
         case 1:
         {
            cell_vdof = L2FESpace.GetNF() + ci;

            pmesh->GetElementVertices(ci, verts);

            for (int j = 0; j < verts.Size(); j++)
            {
               cell_x[0] += x_gf[verts[j]];
            }
            cell_x /= verts.Size();
            x_gf[cell_vdof] = cell_x[0];
            break;
         }
         case 2:
         {
            cell_vdof = H1FESpace.GetNVDofs() + L2FESpace.GetNF() + ci;

            pmesh->GetElementVertices(ci, verts);

            for (int j = 0; j < verts.Size(); j++)
            {
               cell_x[0] += x_gf[verts[j]];
               int index = verts[j] + H1FESpace.GetNDofs();
               cell_x[1] += x_gf[index];
            }
            cell_x /= verts.Size();

            x_gf[cell_vdof] = cell_x[0];
            int index = cell_vdof + H1FESpace.GetNDofs();
            x_gf[index] = cell_x[1];
            break;
         }
         default:
         {
            MFEM_ABORT("Incorrect dim value provided.\n");
         }
         }
      } // end cell center average

      // Update pmesh reference grid function
      pmesh->NewNodes(x_gf, false);
      cout << "Mesh distorted\n";
   } // End distort mesh

   // Initialize mesh velocity
   ConstantCoefficient zero(0.0);
   mv_gf.ProjectCoefficient(zero);

   // Change class variables into static std::functions since virtual static member functions are not an option
   // and Coefficient class requires std::function arguments
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> sv0_static = 
      std::bind(&ProblemBase<dim>::sv0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<void(const Vector &, const double, Vector &)> v0_static = 
      std::bind(&ProblemBase<dim>::v0, problem_class, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
   std::function<double(const Vector &,const double)> ste0_static = 
      std::bind(&ProblemBase<dim>::ste0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> rho0_static = 
      std::bind(&ProblemBase<dim>::rho0, problem_class, std::placeholders::_1, std::placeholders::_2);
   std::function<double(const Vector &,const double)> p0_static = 
      std::bind(&ProblemBase<dim>::p0, problem_class, std::placeholders::_1, std::placeholders::_2);

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
   if (problem == 6)
   {
      double blast_energy = 0.25;
      double blast_position[] = {0.0, 0.0, 0.0};
      DeltaCoefficient e_coeff(blast_position[0], blast_position[1],
                               blast_position[2], blast_energy);
      ste_gf.ProjectCoefficient(e_coeff);
   }
   else
   {
      ste_gf.ProjectCoefficient(ste_coeff);
   }

   // PLF to build mass vector
   FunctionCoefficient rho_coeff(rho0_static); 
   rho_coeff.SetTime(t_init);
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(rho_coeff));
   m->Assemble();

   FunctionCoefficient p_coeff(p0_static);
   p_coeff.SetTime(t_init);

   /* Create Lagrangian Low Order Solver Object */
   LagrangianLOOperator<dim> hydro(H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, offset, use_viscosity, mm, CFL);

   /* Set parameters of the LagrangianLOOperator */
   hydro.SetMVOption(mv_option);
   hydro.SetMVLinOption(do_mv_linearization);
   hydro.SetFVOption(fv_option);
   hydro.SetViscOption(greedy);
   hydro.SetProblem(problem);
   hydro.SetMeshCheck(check_mesh);
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
   socketstream vis_rho_ex, vis_v_ex, vis_ste_ex, vis_p_ex;
   socketstream vis_rho_err, vis_v_err, vis_ste_err, vis_p_err;
   char vishost[] = "localhost";
   int  visport   = 19916;

   ParGridFunction rho_gf(&L2FESpace), rho_cont_gf(&H1cFESpace);
   ParGridFunction mc_gf(&L2FESpace); // Gridfunction to show mass conservation
   mc_gf = 0.;  // if a cells value is 0, mass is conserved

   /* Gridfunctions used in Triple Point */
   ParGridFunction press_gf(&L2FESpace), gamma_gf(&L2FESpace);

   if (problem == 6 || problem == 12)
   {
      Vector U(dim+2);
      for (int i = 0; i < press_gf.Size(); i++)
      {
         hydro.GetCellStateVector(S, i, U);
         double pressure = problem_class->pressure(U, pmesh->GetAttribute(i));
         press_gf[i] = pressure;
         gamma_gf[i] = problem_class->get_gamma(pmesh->GetAttribute(i));
      }
   }

   // Compute the density if we need to visualize it
   if (visualization || gfprint || visit || pview)
   {
      // Compute Density
      for (int i = 0; i < sv_gf.Size(); i++)
      {
         rho_gf[i] = 1./sv_gf[i];
      } 
      // Continuous projection
      GridFunctionCoefficient rho_gf_coeff(&rho_gf);
      rho_cont_gf.ProjectDiscCoefficient(rho_gf_coeff, mfem::ParGridFunction::AvgType::ARITHMETIC);
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

      if (problem == 6 || problem == 12 || problem == 15)
      {
         Wx += offx;
         VisualizeField(vis_press, vishost, visport, press_gf,
                        "Pressure", Wx, Wy, Ww, Wh);
      }
      if (problem == 12)
      {
         Wx += offx; 
         VisualizeField(vis_gamma, vishost, visport, gamma_gf,
                        "Gamma", Wx, Wy, Ww, Wh);
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
            rho_ex_gf.ProjectCoefficient(rho_coeff);
            vel_ex_gf.ProjectCoefficient(v_coeff);
            ste_ex_gf.ProjectCoefficient(ste_coeff);
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

            if (problem == 6 || problem == 12)
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
         
         if (problem == 6 || problem == 12 || problem == 15)
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
   std::ostringstream mesh_name, rho_name, v_name, ste_name, mv_name;
   mesh_name << gfprint_path 
               << setfill('0') 
               << setw(6)
               << 0
               << ".mesh";
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

   std::ofstream mv_ofs(mv_name.str().c_str());
   mv_ofs.precision(8);
   mv_gf.SaveAsOne(mv_ofs);
   mv_ofs.close();

   /* Print gamma/pressure grid function for Triple Point problem */
   if (problem == 6 || problem == 12 || problem == 15) 
   {
      std::ostringstream _press_name;
      _press_name << gfprint_path 
                  << "press.gf";

      std::ofstream press_ofs(_press_name.str().c_str());
      press_ofs.precision(8);
      press_gf.SaveAsOne(press_ofs);
      press_ofs.close();
   }
   if (problem == 12)
   {
      std::ostringstream gamma_name;
      gamma_name << gfprint_path 
                  << "gamma.gf";

      std::ofstream gamma_ofs(gamma_name.str().c_str());
      gamma_ofs.precision(8);
      gamma_gf.SaveAsOne(gamma_ofs);
      gamma_ofs.close();
   }

   VisItDataCollection visit_dc(_visit_basename, pmesh);
   visit_dc.RegisterField("Specific Volume", &sv_gf);
   visit_dc.RegisterField("Density", &rho_gf);
   visit_dc.RegisterField("Density c", &rho_cont_gf);
   visit_dc.RegisterField("Velocity", &v_gf);
   visit_dc.RegisterField("Mesh Velocity", &mv_gf);
   visit_dc.RegisterField("Specific Total Energy", &ste_gf);
   visit_dc.RegisterField("Mass Loss", &mc_gf);
   visit_dc.RegisterField("Pressure", &press_gf);
   visit_dc.RegisterField("Gamma", &gamma_gf);
   visit_dc.SetCycle(0);
   visit_dc.SetTime(0.0);
   visit_dc.Save();

   ParaViewDataCollection paraview_dc(_pview_basename, pmesh);
   paraview_dc.SetLevelsOfDetail(order_mv);
   paraview_dc.SetDataFormat(VTKFormat::ASCII); // BINARY also an option, easier to debug
   paraview_dc.RegisterField("Specific Volume", &sv_gf);
   paraview_dc.RegisterField("Density", &rho_gf);
   paraview_dc.RegisterField("Density c", &rho_cont_gf);
   paraview_dc.RegisterField("Velocity", &v_gf);
   paraview_dc.RegisterField("Mesh Velocity", &mv_gf);
   paraview_dc.RegisterField("Specific Total Energy", &ste_gf);
   paraview_dc.RegisterField("Mass Loss", &mc_gf);
   paraview_dc.RegisterField("Pressure", &press_gf);
   paraview_dc.RegisterField("Gamma", &gamma_gf);
   paraview_dc.SetCycle(0);
   paraview_dc.SetTime(0.0); 
   paraview_dc.Save();

   // Perform the time-integration by looping over time iterations
   // ti with a time step dt.  The main function call here is the
   // LagrangianLOOperator<dim>::MakeTimeStep() funcall.
   double t = t_init;

   if (convergence_testing)
   {
      dt = hmin; // Set timestep to smalled h val for convergence testing
   }

   bool last_step = false;
   BlockVector S_old(S);

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
      hydro.chrono_dij.Start();
      hydro.BuildDijMatrix(S);
      /* Check if we need to change CFL */
      if (problem_class->get_cfl_change() && t > problem_class->get_cfl_time_change() && hydro.GetCFL() != problem_class->get_cfl_second())
      {
         cout << "Changing CFL for Saltzman at time = " << t << endl;
         double CFL_new = problem_class->get_cfl_second();
         hydro.SetCFL(CFL_new);
      }

      if (optimize_timestep)
      {
         hydro.CalculateTimestep(S);
         dt = hydro.GetTimestep();
      }
      hydro.chrono_dij.Stop();

      if (t + dt >= t_final)
      {
         dt = t_final - t;
         last_step = true;
      }

      if (ti == max_tsteps) { last_step = true; }

      S_old = S;
      
      hydro.MakeTimeStep(S, t, dt, isCollapsed);
      t += dt;

      /* End iteration if the mesh has collapsed */
      if (isCollapsed)
      {
         last_step = true;
      }
      
      if (S_old.GetData() == S.GetData()) { cout << "\t State has not changed with step.\n"; return -1; }

      // Make sure that the mesh corresponds to the new solution state. This is
      // needed, because some time integrators use different S-type vectors
      // and the oper object might have redirected the mesh positions to those.
      pmesh->NewNodes(x_gf, false);

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
                 << ",\tt = " << std::setw(5) << std::setprecision(4) << t
                 << ",\tdt = " << std::setw(5) << std::setprecision(6) << dt
                 << ",\t|e| = " << std::setprecision(10) << std::scientific
                 << sqrt_norm
                 << endl;
         }
         // Fill grid function with mass information
         hydro.CheckMassConservation(S, mc_gf);

         // Turn back on suppression
         if (suppress_output)
         {
            // redirect cout to /dev/null
            cout.rdbuf(file.rdbuf());
         }

         // Compute the density if we need to visualize it
         if (visualization || gfprint || visit || pview)
         {
            // Compute Density
            for (int i = 0; i < sv_gf.Size(); i++)
            {
               rho_gf[i] = 1./sv_gf[i];
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

            if (problem == 6 || problem == 12 || problem == 15) // Visualize pressure
            {
               Vector U(dim+2);
               for (int i = 0; i < press_gf.Size(); i++)
               {
                  hydro.GetCellStateVector(S, i, U);
                  double pressure = problem_class->pressure(U, pmesh->GetAttribute(i));
                  press_gf[i] = pressure;
               }
               Wx += offx;

               // Pressure
               VisualizeField(vis_press, vishost, visport, press_gf,
                              "Pressure", Wx, Wy, Ww, Wh);
            }
            if (problem == 12) // Visualize gamma
            {
               vis_gamma.precision(8);
               ParGridFunction gamma_gf(&L2FESpace);
               for (int i = 0; i < gamma_gf.Size(); i++)
               {
                  gamma_gf[i] = problem_class->get_gamma(pmesh->GetAttribute(i));
               }
               // gamma
               Wx += offx;
               VisualizeField(vis_gamma, vishost, visport, gamma_gf,
                              "Gamma", Wx, Wy, Ww, Wh);
            }

            Wx += offx;
            VisualizeField(vis_mc, vishost, visport, mc_gf,
                           "Mass Conservation",
                           Wx, Wy, Ww, Wh);
            
            Wx = 0;
            Wy += offy;

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
                  rho_ex_gf.ProjectCoefficient(rho_coeff);
                  rho_ex_gf[0] = rho_gf[0];
                  vel_ex_gf.ProjectCoefficient(v_coeff);
                  ste_ex_gf.ProjectCoefficient(ste_coeff);
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
                  if (problem == 6 || problem == 12)
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
               
               if (problem == 6 || problem == 12 || problem == 15)
               {
                  Wx += offx;
                  VisualizeField(vis_p_err, vishost, visport, p_err,
                                 "Error: Pressure", Wx, Wy, Ww, Wh);
               }
               
            }
         }

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

            std::ofstream mv_ofs(mv_name.str().c_str());
            mv_ofs.precision(8);
            mv_gf.SaveAsOne(mv_ofs);
            mv_ofs.close();

            // Print continuous interpolation of density
            GridFunctionCoefficient rho_gf_coeff(&rho_gf);
            rho_cont_gf.ProjectDiscCoefficient(rho_gf_coeff, mfem::ParGridFunction::AvgType::ARITHMETIC);
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


            if (visit)
            {
               visit_dc.SetCycle(ti);
               visit_dc.SetTime(t);
               visit_dc.Save();
            }

            if (pview)
            {
               paraview_dc.SetCycle(ti);
               paraview_dc.SetTime(t);
               paraview_dc.Save();
            }     
         }
      }
   } // End time step iteration

   chrono.Stop();

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
                      << ".out";

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
      std::ostringstream mesh_name, rho_name, v_name, ste_name, massC_name, mv_name;
      mesh_name << gfprint_path 
                << "final.mesh";
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

      std::ofstream mv_ofs(mv_name.str().c_str());
      mv_ofs.precision(8);
      mv_gf.SaveAsOne(mv_ofs);
      mv_ofs.close();

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

      /* Always save last VisIt files */
      visit_dc.SetCycle(ti);
      visit_dc.SetTime(t);
      visit_dc.Save();

      /* Always save last ParaView files */
      paraview_dc.SetCycle(ti);
      paraview_dc.SetTime(t);
      paraview_dc.Save();
   } // End print final grid functions

   if (problem_class->has_exact_solution())
   {
      // Save exact to file as well
      BlockVector S_exact(offset, Device::GetMemoryType());
      ParGridFunction x_gf_exact, mv_gf_exact, sv_ex_gf, vel_ex_gf, ste_ex_gf;

      x_gf_exact.MakeRef(&H1FESpace, S_exact, offset[0]);
      mv_gf_exact.MakeRef(&H1FESpace, S_exact, offset[1]);
      sv_ex_gf.MakeRef(&L2FESpace, S_exact, offset[2]);
      vel_ex_gf.MakeRef(&L2VFESpace, S_exact, offset[3]);
      ste_ex_gf.MakeRef(&L2FESpace, S_exact, offset[4]);

      // Project exact solution
      if (problem_class->get_indicator() == "Vdw1")
      {
         problem_class->update(x_gf, t);
      }

      sv_coeff.SetTime(t);
      v_coeff.SetTime(t);
      ste_coeff.SetTime(t);

      sv_ex_gf.ProjectCoefficient(sv_coeff);
      vel_ex_gf.ProjectCoefficient(v_coeff);
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
         rho_ex_gf.ProjectCoefficient(rho_coeff);
         vel_ex_gf.ProjectCoefficient(v_coeff);
         ste_ex_gf.ProjectCoefficient(ste_coeff);
         sv_ex_gf.ProjectCoefficient(sv_coeff);

         // In the case of the Noh Problem, project 0 on the boundary of approx and exact
         if (problem_class->get_indicator() == "Noh")
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

         /* Compute relative errors */
         GridFunctionCoefficient rho_ex_coeff(&rho_ex_gf), vel_ex_coeff(&vel_ex_gf), ste_ex_coeff(&ste_ex_gf), sv_ex_coeff(&sv_ex_gf);

         // Velocity errors
         vel_L1_error_n = v_gf.ComputeL1Error(vel_ex_coeff) / vel_ex_gf.ComputeL1Error(zero);
         vel_L2_error_n = v_gf.ComputeL2Error(vel_ex_coeff) / vel_ex_gf.ComputeL2Error(zero);
         vel_Max_error_n = v_gf.ComputeMaxError(vel_ex_coeff) / vel_ex_gf.ComputeMaxError(zero);

         // Specific volume and specific total energy
         /* In the Sedov case, we do not get convergence of the specific total energy, but the total energy */
         if (problem == 6)
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
                        << "mass_loss " << hydro.CalcMassLoss(S) << "\n"
                        << "dt " << dt << "\n"
                        << "Endtime " << t << "\n";
                  
      convergence_file.close();
   } // Writing convergence file
   
   if (suppress_output)
   {
      // restore cout stream buffer
      cout.rdbuf(strm_buffer);
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
   delete pmesh0;
   delete CRFEC;
   delete problem_class;
   delete m;

   return 0;
}
