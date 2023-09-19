#include "mfem.hpp"
#include "laglos_solver.hpp"
#include "problem_template.h"
#include "sod.h"
#include "var-config.h"
#include <cassert>
#include <fstream>
#include <sstream>
#include <cmath>
#include "matplotlibcpp.h"

using namespace std;
using namespace mfem;
using namespace hydrodynamics;
namespace plt = matplotlibcpp;

/* ---------------- Parameters to be used for tests ---------------- */
// Linear velocity field
double a = 1.,
       b = 0.,
       c = 0.,
       d = 0.,
       e = 0,
       f = 0.,
       aq = 0.,
       bq = 0.,
       dq = 0.,
       eq = 0.;

// Problem
const int dim = 2;
const int problem = 1;     // This value doesn't matter (1-8)
int order_mv = 2;          // Order of mesh movement approximation space
int order_u = 0;
const string flag = "testing";
double tol = 1e-12;
const char *mesh_file_location = "/data/ref-square.mesh";
std::string result = std::string(LAGLOS_DIR) + std::string(mesh_file_location);
const char* mesh_file = result.c_str();
int mesh_refinements = 2;
bool use_viscosity = true; // Doesn't matter
bool _mm = false;          // Doesn't matter
double CFL = 0.5;          // Doesn't matter
static int num_procs = 0;
static int myid = -1;
bool suppress_test_output = false;
bool create_table = false;
int lower_refinement = 2;
int upper_refinement = 7;
/* ---------------- End Parameters ---------------- */

void velocity_exact(const Vector &x, const double & t, Vector &v);
int test_Ci_geo();

void plot_velocities();


int main(int argc, char *argv[])
{
   // Initialize MPI.
   Mpi::Init();
   num_procs = Mpi::WorldSize();
   myid = Mpi::WorldRank();

   // Set precision for all output values
   cout.precision(12);

   // Parse command line arguments
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&mesh_refinements, "-r", "--refinement", "Number of mesh refinements.");

   args.Parse();

   int d = 0;
   d += test_Ci_geo();

   plot_velocities();

   cout << "d: " << d << endl;

   return d;
}

void velocity_exact(const Vector &x, const double & t, Vector &v)
{
   if (dim != 2)
   {
      MFEM_ABORT("Dimension must be equal to 2 for mesh motion test.\n");
   }
   v[0] = aq * pow(x[0],2) + bq * pow(x[1],2) + a * x[0] + b * x[1] + c;
   v[1] = dq * pow(x[0],2) + eq * pow(x[1],2) + d * x[0] + e * x[1] + f;
}


/*
Purpose:
   In a linear velocity field, Cgeo should be the same at every node.  This 
   iterates over all nodes and validates this.
*/
int test_Ci_geo()
{
   cout << "test_Ci_geo\n";
   aq = 0., bq = 0., dq = 0., eq = 0.;
   a = 2., b = 4., c = 6., d = 8., e = 9., f = 10.;
   // Initialize mesh
   mfem::Mesh *mesh;
   mesh = new mfem::Mesh(mesh_file, true, true);

   // Refine the mesh
   for (int lev = 0; lev < mesh_refinements; lev++)
   {
      mesh->UniformRefinement();
   }

   // Construct parmesh object
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);            
   delete mesh;

   // Template FE stuff to construct hydro operator
   H1_FECollection H1FEC(order_mv, dim);
   H1_FECollection H1FEC_L(1, dim);
   L2_FECollection L2FEC(order_u, dim, BasisType::Positive);
   RT_FECollection CRFEC(0, dim);

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace H1FESpace_L(pmesh, &H1FEC_L, dim);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);
   ParFiniteElementSpace CRFESpace(pmesh, &CRFEC, dim);

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

   // Pair with corresponding gridfunctions
   // Only need position and velocity for testing
   ParGridFunction x_gf, mv_gf, sv_gf, v_gf, ste_gf;
   ParGridFunction v_cr_gf(&CRFESpace);
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   mv_gf.MakeRef(&H1FESpace, S, offset[1]);
   sv_gf.MakeRef(&L2FESpace, S, offset[2]);
   v_gf.MakeRef(&L2VFESpace, S, offset[3]);
   ste_gf.MakeRef(&L2FESpace, S, offset[4]);

   pmesh->SetNodalGridFunction(&x_gf);
   x_gf.SyncAliasMemory(S);

   Vector one(dim);
   one = 1.;

   VectorFunctionCoefficient v_exact_coeff(dim, &velocity_exact);
   mv_gf.ProjectCoefficient(v_exact_coeff);
   mv_gf.SyncAliasMemory(S);

   // Initialize specific volume, velocity, and specific total energy
   ConstantCoefficient one_const_coeff(1.0);
   sv_gf.ProjectCoefficient(one_const_coeff);
   sv_gf.SyncAliasMemory(S);

   VectorConstantCoefficient one_coeff(one);
   v_gf.ProjectCoefficient(one_coeff);
   v_gf.SyncAliasMemory(S);

   ste_gf.ProjectCoefficient(one_const_coeff);
   ste_gf.SyncAliasMemory(S);

   // Just leave templated for hydro construction
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(one_const_coeff));
   m->Assemble();

   ProblemBase<dim> * problem_class = new ProblemTemplate<dim>();

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   // Necessary parameters
   DenseMatrix Cgeo(dim), Cgeo_0(dim), C_error(dim);
   double _error = 0.;
   double dt = 1., t = 0.;

   // Must compute the intermediate face velocities before computing geometric velocity
   hydro.ComputeIntermediateFaceVelocities(S, t, "testing", &velocity_exact);

   // Iterate across all nodes and verify Cigeo is the same at all the geometric vertices.
   hydro.ComputeCiGeoRaviart(0, Cgeo_0);
   for (int node_it = 1; node_it < H1FESpace.GetNDofs() - L2FESpace.GetNDofs(); node_it++)
   {
      hydro.ComputeCiGeoRaviart(node_it, Cgeo);
      Add(Cgeo, Cgeo_0, -1., C_error);
      _error += C_error.FNorm();

      cout << "error on node: " << _error << endl;
      cout << "--- Ci on node: " << node_it << " ---\n";
      Cgeo.Print(cout);
   }

   if (abs(_error) < 1e-10)
   {
      return 0; // Test passed
   }
   else 
   {
      cout << "failed test_Ci_geo\n";
      cout << "error: " << _error << endl;
      return 1; // Test failed
   }
}


void plot_velocities()
{
   a = 1., b = 0., c = 0., d = 0., e = -0.5, f = 0.;
   aq = 0., bq = 0., dq = 0., eq = 0.;

   // Initialize mesh
   mfem::Mesh *mesh;
   mesh = new mfem::Mesh(mesh_file, true, true);

   // Refine the mesh
   for (int lev = 0; lev < mesh_refinements; lev++)
   {
      mesh->UniformRefinement();
   }

   // Construct parmesh object
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);            
   delete mesh;

   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-TestVDLFGI." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   // Template FE stuff to construct hydro operator
   H1_FECollection H1FEC(order_mv, dim);
   H1_FECollection H1FEC_L(1, dim);
   L2_FECollection L2FEC(order_u, dim, BasisType::Positive);
   RT_FECollection CRFEC(0, dim);

   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, dim);
   ParFiniteElementSpace H1FESpace_L(pmesh, &H1FEC_L, dim);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace L2VFESpace(pmesh, &L2FEC, dim);
   ParFiniteElementSpace CRFESpace(pmesh, &CRFEC, dim);

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

   // Pair with corresponding gridfunctions
   // Only need position and velocity for testing
   ParGridFunction x_gf, mv_gf, sv_gf, v_gf, ste_gf;
   
   x_gf.MakeRef(&H1FESpace, S, offset[0]);
   mv_gf.MakeRef(&H1FESpace, S, offset[1]);
   sv_gf.MakeRef(&L2FESpace, S, offset[2]);
   v_gf.MakeRef(&L2VFESpace, S, offset[3]);
   ste_gf.MakeRef(&L2FESpace, S, offset[4]);

   pmesh->SetNodalGridFunction(&x_gf);
   x_gf.SyncAliasMemory(S);

   Vector one(dim);
   one = 1.;
   VectorConstantCoefficient one_coeff(one);
   mv_gf.ProjectCoefficient(one_coeff);
   mv_gf.SyncAliasMemory(S);

   // Initialize specific volume, velocity, and specific total energy
   ConstantCoefficient one_const_coeff(1.0);
   sv_gf.ProjectCoefficient(one_const_coeff);
   sv_gf.SyncAliasMemory(S);

   v_gf.ProjectCoefficient(one_coeff);
   v_gf.SyncAliasMemory(S);

   ste_gf.ProjectCoefficient(one_const_coeff);
   ste_gf.SyncAliasMemory(S);

   // Just leave templated for hydro construction
   ParLinearForm *m = new ParLinearForm(&L2FESpace);
   m->AddDomainIntegrator(new DomainLFIntegrator(one_const_coeff));
   m->Assemble();

   ProblemBase<dim> * problem_class = new ProblemTemplate<dim>();

   mfem::hydrodynamics::LagrangianLOOperator<dim> hydro(H1FESpace, H1FESpace_L, L2FESpace, L2VFESpace, CRFESpace, m, problem_class, use_viscosity, _mm, CFL);

   Vector _vel(dim), true_vel(dim);
   double t = 0.;

   hydro.ComputeIntermediateFaceVelocities(S, t, "testing", &velocity_exact);
   hydro.ComputeGeoVRaviart(S);
   
   ParGridFunction v_cr_gf(&CRFESpace);
   ParGridFunction v_geo_gf(&H1FESpace_L);
   hydro.GetVCRgf(v_cr_gf);
   hydro.GetVGeogf(v_geo_gf);

   VectorFunctionCoefficient v_exact_coeff(dim, &velocity_exact);
   ParGridFunction v_exact_gf(&H1FESpace);
   v_exact_gf.ProjectCoefficient(v_exact_coeff);

   // Set int rule
   IntegrationRules IntRulesLo(0, Quadrature1D::GaussLobatto);
   const IntegrationRule gir = IntRules.Get(CRFESpace.GetFE(0)->GetGeomType(), 1);

   // Output integration rule information
   cout << "Order of integration rule: " << gir.GetOrder() << endl;
   cout << "number of integration points: " << gir.GetNPoints() << endl;
   cout << "weights:\n";
   const Array<double> weights = gir.GetWeights();
   for (int i = 0; i < gir.GetNPoints(); i++)
   {
      cout << "weight at i = " << i << ": " << weights[i] << endl;
      IntegrationPoint point = gir.IntPoint(i);
      cout << "x: " << point.x << ", y: " << point.y << endl;
   }
   ParFiniteElementSpace CRc(CRFESpace.GetParMesh(), CRFESpace.FEColl(), 1);
   const int size = CRc.GetVSize();
   cout << "Size: " << size << endl;

   cout << "Diagnostics:\n";
   cout << "CRFESpace.GetVSize(): " << CRFESpace.GetVSize() << endl;
   cout << "CRFESpace.TrueVSize(): " << CRFESpace.TrueVSize() << endl;
   cout << "CRFESpace.GetNDofs(): " << CRFESpace.GetNDofs() << endl;

   cout << "CRc.GetVSize(): " << CRc.GetVSize() << endl;
   cout << "CRc.TrueVSize(): " << CRc.TrueVSize() << endl;
   cout << "CRc.GetNDofs(): " << CRc.GetNDofs() << endl;

   Array<double> cell_areas(L2FESpace.GetNE());
   Array<double> bmn_sum_cells(L2FESpace.GetNE());

   // Store initial mass of cells
   for (int ci = 0; ci < L2FESpace.GetNE(); ci++)
   {
      double k = pmesh->GetElementVolume(ci);

      cell_areas[ci] = k;
      assert(k >= 0.);
   }

   DenseMatrix res(dim);

   Vector x(dim), vec_res(dim);
   double dt = .001;
   bool is_dt_changed = false;
   for (int node_it = 0; node_it < H1FESpace.GetNDofs() - L2FESpace.GetNDofs(); node_it++)
   {
      // hydro.GetNodePosition(S, node_it, x);
      // cout << "node position for node: " << node << endl;
      // x.Print(cout);

      hydro.ComputeNodeVelocityRaviart(node_it, dt, vec_res, is_dt_changed);
      hydro.UpdateNodeVelocity(S, node_it, vec_res);

      // restart nodal velocity computation if the timestep has been restricted
      if (is_dt_changed)
      {
         is_dt_changed = false;
         node_it = -1;
      }
   }

   // Fill center with average over corner vertices
   Array<int> verts;
   Vector vel_center(dim);
   for (int ci = 0; ci < L2FESpace.GetNDofs(); ci++)
   {
      double vel_center_x = 0., vel_center_y = 0.;
      pmesh->GetElementVertices(ci, verts);
      for (int j = 0; j < verts.Size(); j++)
      {
         int index = verts[j];
         vel_center_x += v_geo_gf[index];
         index = verts[j] + H1FESpace.GetNDofs();
         vel_center_y += v_geo_gf[index];
      }
      vel_center_x *= 1. / verts.Size();
      vel_center_y *= 1. / verts.Size();

      v_geo_gf[H1FESpace.GetNDofs() - L2FESpace.GetNDofs() + ci] = vel_center_x;
      v_geo_gf[2 * H1FESpace.GetNDofs() - L2FESpace.GetNDofs() + ci] = vel_center_y;
   }

   // hydro.ComputeCorrectiveFaceVelocities(S, t, dt);
   hydro.FillFaceVelocitiesWithAvg(S);
   hydro.FillCenterVelocitiesWithAvg(S);

   Array<int> fids, oris, row;
   mfem::Mesh::FaceInformation FI;
   double bmn_sum;
   int c, cp;
   Vector Uc(dim+2), Ucp(dim+2), c_vec(dim), Vf(dim), face_x(dim), vdof1_x(dim), vdof2_x(dim), n_int(dim), n_vec(dim);

   // for (int ci = 0; ci < L2FESpace.GetNE(); ci++)
   // {
   //    // Compute sum of bmn on cell
   //    pmesh->GetElementEdges(ci, fids, oris);

   //    // iterate over faces
   //    bmn_sum = 0.;
   //    for (int j = 0; j < fids.Size(); j++)
   //    {
   //       Vf = 0.;
   //       int face = fids[j];
   //       FI = pmesh->GetFaceInformation(face);

   //       // Compute intermediate face velocity on fly
   //       c = FI.element[0].index;
   //       cp = FI.element[1].index;
   //       hydro.GetCellStateVector(S, c, Uc);

   //       H1FESpace.GetFaceDofs(fids[j], row);
   //       int face_vdof1 = row[1], face_vdof2 = row[0], face_dof = row[2];
   //       hydro.GetNodePosition(S, face_vdof1, vdof1_x);
   //       hydro.GetNodePosition(S, face_vdof2, vdof2_x);
   //       hydro.GetNodePosition(S, face_dof, face_x);
         
   //       velocity_exact(face_x, t, Vf);

   //       hydro.CalcOutwardNormalInt(S, ci, face, n_int);
   //       n_vec = n_int;
   //       double F = n_vec.Norml2();
   //       n_vec /= F;

   //       double bmn = Vf * n_vec;
   //       bmn *= F;
   //       bmn_sum += bmn;
         
   //    } // End face iterator
   //    bmn_sum_cells[ci] = bmn_sum;
   // }

   /* ************************
   Move the mesh according to previously computed nodal velocities
   *************************** */ 
   cout << "x_gf pre move:\n";
   x_gf.Print(cout);
   add(x_gf, dt, mv_gf, x_gf);
   pmesh->NewNodes(x_gf, false);
   cout << "x_gf post move:\n";
   x_gf.Print(cout);
   t += dt;
   if (!suppress_test_output)
   {
      cout << "Done moving mesh.\n";
   }

   double num = 0., error_denom = 0., _error = 0.;
   /* ************************
   Compute cell area error
   *************************** */ 
   std::vector<double> cells_py(L2FESpace.GetNE());
   std::vector<double> cell_errors_py(L2FESpace.GetNE());
   // for (int ci = 0; ci < L2FESpace.GetNE(); ci++)
   // {
   //    if (!suppress_test_output)
   //    {
   //       cout << "Computing cell area error for cell: " << ci << endl;
   //    }
   //    // Compute cell area side
   //    double k = pmesh->GetElementVolume(ci);
   //    double _bmn_approx = (k - cell_areas[ci]) / dt;
   //    cout << "cell area pre move: " << cell_areas[ci] << endl;
   //    cout << "cell area post move: " << k << endl;
   //    cout << "dt: " << dt << endl;

   //    cout << "Cell wise comparison:\n";
   //    cout << "cell area difference: " << _bmn_approx << ", bmn_sum: " << bmn_sum_cells[ci] << endl;

   //    double val = abs(_bmn_approx - bmn_sum_cells[ci]);
   //    num += val;
   //    error_denom += abs(bmn_sum_cells[ci]);
   //    cells_py[ci] = ci;
   //    cell_errors_py[ci] = val;
   // }

   // if (error_denom < tol)
   // {
   //    _error = 0.;
   // }
   // else
   // {
   //    _error = abs(num / error_denom);
   // }

   /* ************************
   Visualization
   *************************** */ 
   // Output mesh to be visualized
   // Can be visualized with glvis -np # -m mesh-test-moved
   {
      int precision = 12;
      ostringstream mesh_name;
      mesh_name << "../results/mesh-TestVDLFGI-moved." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
   }

   socketstream vis_vh, vis_vgeo, vis_vexact;
   char vishost[] = "localhost";
   int visport = 19916;

   MPI_Barrier(pmesh->GetComm());
   vis_vgeo.precision(8);

   int Wx = 0, Wy = 0;
   const int Ww = 350, Wh = 350;
   int offx = Ww+10, offy = Wh+45;

   hydrodynamics::VisualizeField(vis_vgeo, vishost, visport, v_geo_gf, "Geometric velocity", Wx, Wy, Ww, Wh);

   Wx += offx;

   hydrodynamics::VisualizeField(vis_vh, vishost, visport, mv_gf, "Mesh Velocity", Wx, Wy, Ww, Wh);

   Wx += offx;

   hydrodynamics::VisualizeField(vis_vexact, vishost, visport, v_exact_gf, "Exact Velocity", Wx, Wy, Ww, Wh);

   // Plots cell errors using Python's matplotlib library
   if (!suppress_test_output)
   {
      plt::figure_size(1200, 780);
      plt::title("Error by cell");
      plt::plot(cells_py, cell_errors_py, "og");
      plt::xlabel("cell");
      plt::ylabel("error");
      plt::show();
   }
   
}
