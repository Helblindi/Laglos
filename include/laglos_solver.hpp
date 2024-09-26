#ifndef LAGLOS_SOLVER
#define LAGLOS_SOLVER

#include "mfem.hpp"
// #include "initial_vals.hpp"
#include "problem_base.h"
#include "geometry.hpp" // Mesh information
#include "lagrange_multiplier.hpp"
#include "lagrange_multiplier_dense.hpp" // TODO: Remove
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <unordered_set> // For cell face normal velocity function

using namespace std;

namespace mfem
{

namespace hydrodynamics
{

/// Visualize the given parallel grid function, using a GLVis server on the
/// specified host and port. Set the visualization window title, and optionally,
/// its geometry.
void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    ParGridFunction &gf, const char *title,
                    int x = 0, int y = 0, int w = 400, int h = 400,
                    bool vec = false);

struct TimingData
{
   // Total times for all major computations:
   // CG solves (H1 and L2) / force RHS assemblies / quadrature computations.
   StopWatch sw_cgH1, sw_cgL2, sw_force, sw_qdata;

   // Store the number of dofs of the corresponding local CG
   const HYPRE_Int L2dof;

   // These accumulate the total processed dofs or quad points:
   // #(CG iterations) for the L2 CG solve.
   // #quads * #(RK sub steps) for the quadrature data computations.
   HYPRE_Int H1iter, L2iter;
   HYPRE_Int quad_tstep;

   TimingData(const HYPRE_Int l2d) :
      L2dof(l2d), H1iter(0), L2iter(0), quad_tstep(0) { }
};

//
template <int dim>
class LagrangianLOOperator
{
protected:
   ParFiniteElementSpace &H1, &L2, &L2V, &CR, CRc;
   ParFiniteElementSpace &H1_L;
   ParFiniteElementSpace H1Lc;
   ParGridFunction v_CR_gf; // 5.7(b)
   ParGridFunction v_CR_gf_corrected; // Iteratively updated
   ParGridFunction v_CR_gf_fluxes;    // Iteratively updated
   ParGridFunction cell_bdr_flag_gf;  // Element indexing vector
   ParGridFunction LagrangeMultipliers;
   // Stores mesh velocities from the previous iteration 
   Vector lambda_max_vec; // TODO: remove, just for temp plotting
   ParGridFunction v_geo_gf; // 5.11
   ParMesh *pmesh;
   ParLinearForm *m_lf;
   HypreParVector *m_hpv;

   // Problem specific
   ProblemBase<dim> * pb;

   // Geometric
   Geometric<dim> geom;

   // Matrix to hold max wavespeed values dij
   SparseMatrix * dij_sparse;

   // FE spaces local and global sizes
   const int Vsize_H1;
   const int TVSize_H1;
   const HYPRE_Int GTVSize_H1;
   const int NDofs_H1;
   const int NVDofs_H1;
   const int NDofs_H1L;

   const int Vsize_L2;
   const int TVSize_L2;
   const HYPRE_Int GTVSize_L2;
   const int NDofs_L2;

   const int Vsize_L2V;
   const int TVSize_L2V;
   const HYPRE_Int GTVSize_L2V;
   const int NDofs_L2V;

   IntegrationRule RT_ir;
   const int RT_ir_order = 2;

   // Tables to relate cell to the contained faces
   // Ref: https://mfem.org/howto/nav-mesh-connectivity/
   Table element_face;
   Table vertex_edge;
   Table * vertex_element;
   Table * face_element;
   Table * edge_vertex;
   Array<int> block_offsets;
   Array<int> BdrElementIndexingArray; // Array to identify boundary faces
   Array<int> BdrVertexIndexingArray;  // Array to identify boundary vertices

   int el_num_faces;
   const int num_elements, num_vertices, num_faces, num_edges;

   double CFL;
   double timestep = 0.001;
   double timestep_first = 0.; // Set and used for activation function when prescribing left wall dirichlet BCs for Saltzman problem

   bool use_viscosity;
   bool mm;
   bool check_mesh;
   int mv_option = 0;
   bool do_mv_linearization;
   int fv_option = 0;
   int mv_it_option = 2;
   bool use_corner_velocity_MC_iteration = false;
   int corner_velocity_MC_num_iterations = 0;
   double mm_visc_face = 0., mm_cell = 0.;
   int problem = -1;

   Array<int> HiopHessIArr, HiopHessJArr;
   Array<int> HiopCGradIArr, HiopCGradJArr;
   Array<int> HiopDGradIArr, HiopDGradJArr;
   Array<double> HiopDGradData;

   Array<int> ess_bdr, dofs_list, ess_tdofs;
   Array<double> bdr_vals;

public:
   enum DofEntity {corner, face, cell};

   LagrangianLOOperator(ParFiniteElementSpace &h1,
                        ParFiniteElementSpace &h1_l,
                        ParFiniteElementSpace &l2,
                        ParFiniteElementSpace &l2v,
                        ParFiniteElementSpace &cr,
                        ParLinearForm *m,
                        ProblemBase<dim> *_pb,
                        Array<int> offset,
                        bool use_viscosity,
                        bool mm,
                        double CFL);
   ~LagrangianLOOperator();

   double GetCFL() { return this->CFL; }
   double GetTimestep() { return timestep; }
   void SetCFL(const double &_CFL) { this->CFL = _CFL; }

   void SetProblem(const int _problem) { this->problem = _problem; }
   void SetMeshCheck(const bool _check_mesh) { this->check_mesh = _check_mesh; }

   void GetEntityDof(const int GDof, DofEntity & entity, int & EDof);

   void CreateBdrElementIndexingArray();
   void CreateBdrVertexIndexingArray();
   void FillCellBdrFlag();
   void GetCellBdrFlagGF(ParGridFunction &_cell_bdr_flag_gf) { _cell_bdr_flag_gf = this->cell_bdr_flag_gf; }

   bool IsBdrVertex(const int & node) { return (BdrVertexIndexingArray[node] == 1); }

   void MakeTimeStep(Vector &S, const double & t, const double & dt, bool &isCollapsed);
   void ComputeStateUpdate(Vector &S_new, const double &t, const double dt);

   void InitializeDijMatrix();
   void BuildDijMatrix(const Vector &S);
   void CalculateTimestep(const Vector &S);

   void GetCellStateVector(const Vector &S, const int cell, Vector &U);
   void SetCellStateVector(Vector &S_new, const int cell, const Vector &U);

   /* cij comp */
   void CalcOutwardNormalInt(const Vector &S, const int cell, const int face, Vector & res);

   /* System timing */
   StopWatch chrono_mm, chrono_state, chrono_dij, chrono_mm_lin;

   /* Mesh movement */
   void SetMVOption(const int & option);
   void SetMVLinOption(const bool & option) { this->do_mv_linearization = option; }
   void SetFVOption(const int & option);
   void SetMVIterationOption(const int &option);
   void SetMVIteration(const int num_iterations);
   void SetMMViscFace(const double mm_visc);
   void SetMMCell(const double mm_consistency);
   void GetIntermediateFaceVelocity(const int & face, Vector & vel);
   void SetCorrectedFaceVelocity(const int & face, const Vector & vel); 
   void GetCorrectedFaceVelocity(const int & face, Vector & vel);       
   void SetCorrectedFaceFlux(const int & face, const Vector &   ); 
   void GetCorrectedFaceFlux(const int & face, Vector & flux);       

   void SetViGeo(const int &node, const Vector &vel);
   void GetViGeo(const int & node, Vector & vel);
   void GetLambdaMaxVec(Vector &lambda_max_vec) { lambda_max_vec = this->lambda_max_vec; }
   void GetVGeogf(ParGridFunction & _v_geo_gf) { _v_geo_gf = this->v_geo_gf; }
   void GetVCRgf(ParGridFunction & _v_CR_gf) { _v_CR_gf = this->v_CR_gf; }

   void UpdateMesh(const Vector &S) const;
   
   // Face velocity functions
   void ComputeIntermediateFaceVelocities(const Vector &S, 
                                             const double t,
                                             const string ="NA", 
                                             void (*test_vel)(const Vector&, const double&, Vector&) = NULL);
   void ComputeCorrectiveFaceVelocities(Vector &S, const double & t, const double & dt, const string ="NA", void (*test_vel)(const Vector&, const double&, Vector&) = NULL);
   void ComputeCorrectiveFaceFluxes(Vector &S, const double & t, const double & dt);
   void FillFaceVelocitiesWithAvg(Vector &S, const string ="NA", void (*test_vel)(const Vector&, const double&, Vector&) = NULL);

   // Fill mv_gf for cell centers
   void SetCellCenterAsCenter(Vector &S);
   void FillCenterVelocitiesWithL2(Vector &S);
   void FillCenterVelocitiesWithAvg(Vector &S);

   // Compute mesh velocities
   void ComputeMeshVelocities(Vector &S, const Vector &S_old, const double &t, const double &dt);

   // Check Jacobians to ensure mesh hasn't collapsed
   bool IsMeshCollapsed();

   // Normal vector mesh motion
   void tensor(const Vector & v1, const Vector & v2, DenseMatrix & dm);
   void ComputeGeoVNormal(const Vector &S, ParGridFunction &mv_gf_l);

   // Normal vector mesh motion with distributed viscosity (discussed on 09/05/2024)
   void ComputeGeoVNormalDistributedViscosity(Vector &S);
   
   // Raviart-Thomas mesh motion 
   void ComputeGeoVRaviart(Vector &S);
   void ComputeGeoVRaviart2(const Vector &S);

   // Cell Face Normal
   void ComputeGeoVCellFaceNormal(Vector &S);

   // CAVEAT
   void ComputeGeoVCAVEAT(Vector &S);

   // Combo of CAVEAT on boundary and Cell Face Normal on interior vertices
   void ComputeGeoVCAVEATCellFace(Vector &S);

   // The above, but using weighted least squares
   void ComputeGeoVCAVEATCellFaceWeighted(Vector &S);

   // Iterative method to update corner velocities, so less mass correction is needed
   void IterativeCornerVelocityMC(Vector &S, const double & dt);
   double ComputeIterationNormMC(Vector &S, const double & dt);

   // Iterative method using least squares
   void IterativeCornerVelocityLS(Vector &S, const double & dt);
   double compute_alpha(const Vector &Vadj, const Vector &aadj, 
                        const Vector &Vnode, const Vector &anode, 
                        const Vector &V3n, const Vector &a3n, 
                        const Vector &n_vec, const Vector &tau_vec,
                        const double &dt, const double &F);
   double ComputeIterativeLSGamma(Vector &S, const double & dt);
   double ComputeIterationNorm(const Vector &S, const ParGridFunction &mv_gf_prev_it, const double & dt);
   double ComputeFaceSecantNorm(Vector &S, const double & dt);

   // Iterative method using least squares (tangent and normal)
   void IterativeCornerVelocityTNLSnoncart(Vector &S, const double & dt);

   // Iterative method to compute the corner velocities using Least Squares
   // that minimizes the change in mass, plus some optional viscosity
   void IterativeCornerVelocityLSCellVolumeFaceVisc(Vector &S, const Vector &S_old, const double &dt);
   void IterativeCornerVelocityLSCellVolumeCellVisc(Vector &S, const Vector &S_old, const double &dt);
   void IterativeCornerVelocityLSCellVolumeMv2Visc(Vector &S, const Vector &S_old, const ParGridFunction &mv2_gf, const double &dt);
   void IterativeCornerVelocityLSCellVolumeMv2FaceVisc(Vector &S, const Vector &S_old, const ParGridFunction &mv2_gf, const double &dt);
   void VerifyContributions(const Vector &S, const Vector &S_old, const ParGridFunction &mv2_gf, const double &dt, const int &it);
   double ComputeCellVolume(const Vector &S, const int &cell);
   double ComputeCellVolumeNorm(const Vector &S, const Vector &S_old, const double &dt);
   void compare_gamma2(const Vector &S, const Vector &S_old, const double &dt, const int &it);

   // Average Velocities
   void ComputeAverageVelocities(Vector &S);

   // Trying Flux on LHS
   void IterativeCornerVelocityFLUXLS(Vector &S, const double & dt);

   // HiOp Lagrange Multipliers implementation
   void CalcMassVolumeVector(const Vector &S, const Vector &S_old, const double &dt, Vector &massvec);
   void CalcCellAveragedCornerVelocityVector(const Vector &S, const Vector &S_old, const bool &is_weighted, int td_flag, ParGridFunction &mv_gf_l);
   void ComputeCellAverageVelocityAtNode(const Vector &S, const Vector &S_old, const int node, const bool &is_weighted, int td_flag, Vector &node_v);
   void DistributeFaceViscosityToVelocity(const Vector &S, Vector &mv_gf);
   void SolveHiOp(const Vector &S, const Vector &S_old, const int & target_option, const double &t, const double &dt, ParGridFunction &mv_gf_l);
   void SolveHiOpDense(const Vector &S, const Vector &S_old, const int & target_option, const double &t, const double &dt, ParGridFunction &mv_gf_l);
   
   /* Functions to linearize the mesh velocity */
   void ComputeDeterminant(const DenseMatrix &C, const double &dt, double & d, int obj_index);
   void ComputeCiGeo(const ParGridFunction &mv_gf_l, const int & node, DenseMatrix & res);
   void IntGrad(const ParGridFunction &mv_gf_l, const int cell, DenseMatrix & res);
   void ComputeLinearizedNodeVelocities(const ParGridFunction &mv_gf_l, ParGridFunction &mv_gf_linearized, const double &t, const  double &dt, const string ="NA", void (*test_vel)(const Vector&, const double&, Vector&) = NULL);
   void ComputeLinearizedNodeVelocity(const ParGridFunction &mv_gf_l, const int & node, const double & dt, Vector &node_v, bool &is_dt_changed);

   // Enforce Boundary Conditions
   void EnforceExactBCOnCell(const Vector &S, const int & cell, const double &t, 
                             const double &dt, Vector & state_val);
   void EnforceMVBoundaryConditions(Vector &S, const double &t, const double &dt);

   // Validate mass conservation
   double CalcMassLoss(const Vector &S);
   void CheckMassConservation(const Vector &S, ParGridFunction & mc_gf);

   // Various print functions
   void SaveStateVecsToFile(const Vector &S, const string &output_file_prefix, const string &output_file_suffix);
};

} // end ns hydrodynamics

} // end ns mfem

#endif // LAGLOS_SOLVER