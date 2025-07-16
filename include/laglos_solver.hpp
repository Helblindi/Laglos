#ifndef LAGLOS_SOLVER
#define LAGLOS_SOLVER

#include "mfem.hpp"
#include "test_problems_include.h"
#include "geometry.hpp" // Mesh information
#include "lagrange_multiplier.hpp"
#include "lagrange_multiplier_dense.hpp" // TODO: Remove
#include "laglos_assembly.hpp"
#include "elastic.hpp" //NF//MS
#include "laglos_tools.hpp"
#include "mfem/linalg/dtensor.hpp" // For Reshape
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <unordered_set> // For cell face normal velocity function

using namespace std;

namespace mfem
{

namespace hydroLO
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


class LagrangianLOOperator : public LimitedTimeDependentOperator
{
protected:
   const int dim;
   
   ParFiniteElementSpace &H1, &L2, &L2V, &CR, CRc;
   ParFiniteElementSpace &H1_L;
   ParFiniteElementSpace smesh_H1L;
   ParFiniteElementSpace H1Lc;
   const ParGridFunction &rho0_gf;
   mutable ParGridFunction x_gf;
   mutable ParGridFunction smesh_x_gf;
   mutable ParGridFunction mv_gf;
   mutable ParGridFunction v_CR_gf; // 5.7(b)
   ParGridFunction v_CR_gf_corrected; // Iteratively updated
   ParGridFunction v_CR_gf_fluxes;    // Iteratively updated
   ParGridFunction cell_bdr_flag_gf;  // Element indexing vector
   ParGridFunction LagrangeMultipliers;
   // Stores mesh velocities from the previous iteration 
   Vector lambda_max_vec; // TODO: remove, just for temp plotting
   ParGridFunction v_geo_gf; // 5.11
   ParMesh *pmesh, *smesh;
   ParLinearForm *m_lf;
   HypreParVector *m_hpv;

   mutable HypreParVector *k_hpv;

   // Problem specific
   ProblemBase * pb;

   // Geometric
   Geometric geom;

   // Matrix to hold max wavespeed values dij
   Array<int> dij_I, dij_J;
   Array<double> dij_data;
   SparseMatrix * dij_sparse = NULL;

   // Matrices to hold geometric vectors cij
   SparseMatrix * cij_sparse_x = NULL;
   SparseMatrix * cij_sparse_y = NULL;
   SparseMatrix * cij_sparse_z = NULL;
   void FreeCij();

   // FE spaces local and global sizes
   const int Vsize_H1;
   const int TVSize_H1;
   const HYPRE_Int GTVSize_H1;
   const int NDofs_H1;
   const int NVDofs_H1;
   const int NDofs_H1L;
   const int Vsize_H1L;

   const int Vsize_L2;
   const int TVSize_L2;
   const HYPRE_Int GTVSize_L2;
   const int NDofs_L2;
   const int order_t;

   const int Vsize_L2V;
   const int TVSize_L2V;
   const HYPRE_Int GTVSize_L2V;
   const int NDofs_L2V;

   // Mesh information
   const int NE, smesh_NBE;
   IntegrationRule RT_ir;
   const int RT_ir_order = 2;

   // Elasticity //NF//MS
   bool use_elasticity = false;
   const IntegrationRule ir;

   // HO
   Table dof_h1_dof_l2; // Table to relate H1 and L2 dofs
   void BuildDofH1DofL2Table();
   mutable QuadratureData qdata; // Data associated with each quadrature point in the mesh
   const int l2dofs_cnt, l2vdofs_cnt;
   void ComputeMassConservativeDensity(ParGridFunction &rho) const;
   Vector initial_masses, initial_volumes;
   void MassesAndVolumesAtPosition(const ParGridFunction &u, const GridFunction &x,
                                   Vector &el_mass, Vector &el_vol) const;
   mutable DenseTensor Me, Me_inv; // Energy mass matrix and its inverse
   MassIntegrator * mi; // Mass integrator for mass matrix
   void ComputeHydroLocRHS(const Vector &S, const int &el, Vector &loc_tau_rhs, Vector &loc_e_rhs, DenseMatrix &loc_v_rhs) const;
   int FindFaceIndexPmesh(const int &el1, const int &el2) const;
   int FindFaceIndexSmesh(const int &el1, const int &el2) const;

   // Tables to relate cell to the contained faces
   // Ref: https://mfem.org/howto/nav-mesh-connectivity/
   Table element_face;
   Table vertex_edge;
   Table * vertex_element;
   Table * face_element;
   Table * edge_vertex;

   Table * smesh_vertex_element;
   Table * smesh_face_element;
   Table * smesh_edge_vertex;
   Table smesh_vertex_edge;
   Array<int> block_offsets;
   Array<int> smesh_BdrElementIndexingArray; // Array to identify boundary faces
   Array<int> smesh_BdrVertexIndexingArray;  // Array to identify boundary vertices

   int el_num_faces;
   const int num_vertices, num_faces;
   const int smesh_num_faces;

   double CFL;
   double timestep = 0.001;
   mutable double timestep_first = 0.; // Set and used for activation function when prescribing left wall dirichlet BCs for Saltzman problem

   int visc;
   bool mm;
   bool compute_mv = true;
   bool use_greedy_viscosity;
   bool post_process_density;
   int mv_option = 0;
   bool do_mv_linearization;
   int fv_option = 0;
   int mv_it_option = 2;
   bool use_corner_velocity_MC_iteration = false;
   int corner_velocity_MC_num_iterations = 0;
   double mm_visc_face = 0., mm_cell = 0.;
   double mv_target_visc_coeff = 0.;
   int problem = -1;

   Array<int> HiopHessIArr, HiopHessJArr;
   Array<int> HiopCGradIArr, HiopCGradJArr;
   Array<int> HiopDGradIArr, HiopDGradJArr;
   Array<double> HiopDGradData;

   bool is_L2_connectivity_built = false;
   int L2ConnectivitySize = 0;
   Table L2Connectivity;

   int ess_tdofs_cart_size;
   Array<int> ess_bdr, dofs_list, ess_tdofs;
   mutable Array<double> bdr_vals;
   Array<int> add_ess_tdofs;
   mutable Array<double> add_bdr_vals;

   /* Time series data */
   Array<double> ts_timestep, ts_t, ts_dijmax, ts_dijavg, ts_min_detJ, ts_min_detJ_cell;
   Array<double> ts_kidder_avg_rad_ext, ts_kidder_avg_rad_int, ts_kidder_avg_density, ts_kidder_avg_entropy;

public:
   Elastic *elastic; //NF//MS
   enum DofEntity {corner, face, cell};

   LagrangianLOOperator(const int &_dim,
                        const int size,
                        ParFiniteElementSpace &h1,
                        ParFiniteElementSpace &h1_l,
                        ParFiniteElementSpace &l2,
                        ParFiniteElementSpace &l2v,
                        ParFiniteElementSpace &cr,
                        Coefficient &rho0_coeff,
                        const ParGridFunction &rho0_gf,
                        ParLinearForm *m,
                        const IntegrationRule &_ir,
                        ProblemBase *_pb,
                        Array<int> offset,
                        int visc,
                        int elastic_eos,
                        bool mm,
                        double CFL);
   ~LagrangianLOOperator();

   virtual void MultUnlimited(const Vector &S, Vector &dS_dt) const;
   virtual void LimitMult(const Vector &S, Vector &dS_dt) const;

   const Array<int> &GetBlockOffsets() const { return block_offsets; }

   /* Limiter */
   ParFiniteElementSpace & GetL2FE() { return L2; }
   ParFiniteElementSpace & GetL2VFE() { return L2V; }
   ParFiniteElementSpace & GetH1FE() { return H1; }

   double GetCFL() { return this->CFL; }
   double GetTimestep() const { return timestep; }
   void SetCFL(const double &_CFL) { this->CFL = _CFL; }

   void SetProblem(const int _problem) { this->problem = _problem; }
   void SetDensityPP(const bool _post_process_density) { this->post_process_density = _post_process_density; }
   bool GetDensityPP() const { return this->post_process_density; }
   void SetViscOption(const bool _use_greedy_viscosity) { this->use_greedy_viscosity = _use_greedy_viscosity; }
   void GetEntityDof(const int GDof, DofEntity & entity, int & EDof);

   void CreateBdrElementIndexingArray();
   void CreateBdrVertexIndexingArray();
   void FillCellBdrFlag();
   void GetCellBdrFlagGF(ParGridFunction &_cell_bdr_flag_gf) { _cell_bdr_flag_gf = this->cell_bdr_flag_gf; }

   bool IsBdrVertex(const int & node) { return (smesh_BdrVertexIndexingArray[node] == 1); }

   void SolveHydro(const Vector &S, Vector &dS_dt) const;
   void EnforceL2BC(Vector &S, const double &t, const double &dt);

   void InitializeDijMatrix();
   void BuildDijMatrix(const Vector &S);
   void CalculateTimestep(const Vector &S);

   void GetStateVector(const Vector &S, const int &index, Vector &U) const;
   void SetStateVector(Vector &S_new, const int &index, const Vector &U) const;

   /* cij comp */
   void CheckSkewSymmetry() const;
   void BuildCijMatrices();
   void GetLocalCij(const int &i, const int &j, Vector &cij) const;
   void BuildL2ConnectivityTable();
   void GetL2ConnectivityTable(Table &_L2Connectivity) const;
   void CalcOutwardNormalInt(const Vector &S, const int cell, const int face, Vector & res) const;

   /* System timing */
   mutable StopWatch chrono_mm, chrono_state, chrono_dij, chrono_mm_lin, chrono_hiop;

   //NF//MS
   void SetElasticity(const bool _use_elasticity) { this->use_elasticity = _use_elasticity; }
   void ComputeSigmaDComp(const Vector &S, const int &e, DenseMatrix &sigma_e) const;
   void ComputeSigmaGF(const Vector &S, ParGridFunction &sigma_gf) const;
   void ComputeFGF(ParGridFunction &f_gf) const;
   void ComputeESheerGF(ParGridFunction &e_sheer_gf) const;
   void SetShearModulus(const double &_mu) { elastic->set_shear_modulus(_mu); }

   /* Mesh movement */
   void UpdateMeshVelocityBCs(const double &t, const double &dt);
   void SolveMeshVelocities(const Vector &S, Vector &dS_dt) const;
   void SetMVTargetViscCoeff(const double & coeff);
   void SetMVOption(const int & option);
   void SetComputeMV(const bool & option) { this->compute_mv = option; }
   void SetMV(const ParGridFunction &_mv) { this->mv_gf = _mv; }
   void GetMV(ParGridFunction & _mv) const { _mv = this->mv_gf; }
   void SetMVLinOption(const bool & option) { this->do_mv_linearization = option; }
   void SetFVOption(const int & option);
   void SetMVIterationOption(const int &option);
   void SetMVIteration(const int num_iterations);
   void SetMMViscFace(const double mm_visc);
   void SetMMCell(const double mm_consistency);
   void GetIntermediateFaceVelocity(const int & face, Vector & vel) const;
   void SetCorrectedFaceVelocity(const int & face, const Vector & vel); 
   void GetCorrectedFaceVelocity(const int & face, Vector & vel);       
   void SetCorrectedFaceFlux(const int & face, const Vector &   ); 
   void GetCorrectedFaceFlux(const int & face, Vector & flux);   
   // HO
   void ComputeMVHOfirst(const Vector &S, ParGridFunction &dxdt_gf) const;

   void SetViGeo(const int &node, const Vector &vel);
   void GetViGeo(const int & node, Vector & vel);
   void GetLambdaMaxVec(Vector &lambda_max_vec) { lambda_max_vec = this->lambda_max_vec; }
   void GetVGeogf(ParGridFunction & _v_geo_gf) { _v_geo_gf = this->v_geo_gf; }
   void GetVCRgf(ParGridFunction & _v_CR_gf) { _v_CR_gf = this->v_CR_gf; }

   void UpdateMesh(const Vector &S) const;
   
   // Face velocity functions
   void ComputeIntermediateFaceVelocities(const Vector &S) const;
   void ComputeCorrectiveFaceVelocities(Vector &S, const double & t, const double & dt, const string ="NA", void (*test_vel)(const Vector&, const double&, Vector&) = NULL);
   void ComputeCorrectiveFaceFluxes(Vector &S, const double & t, const double & dt);
   void FillFaceVelocitiesWithAvg(ParGridFunction &dxdt) const;
   void FillFaceVelocitiesWithButterfly(ParGridFunction &dxdt) const;

   // Fill mv_gf for cell centers
   void SetCellCenterAsCenter(Vector &S);
   void FillCenterVelocitiesWithL2(const Vector &S, Vector &dSdt) const;
   void FillCenterVelocitiesWithAvg(Vector &dxdt) const;

   // Normal vector mesh motion
   void tensor(const Vector & v1, const Vector & v2, DenseMatrix & dm) const;
   void ComputeGeoVNormal(const Vector &S, ParGridFunction &dxdt_gf) const;

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
   void CalcMassVolumeVector(const Vector &S, const double &dt, Vector &massvec);
   void CalcCellAveragedCornerVelocityVector(const Vector &S, const Vector &S_old, const bool &is_weighted, int td_flag, ParGridFunction &mv_gf_l);
   void ComputeCellAverageVelocityAtNode(const Vector &S, const Vector &S_old, const int node, const bool &is_weighted, int td_flag, Vector &node_v);
   void DistributeFaceViscosityToVelocity(const Vector &S, Vector &mv_gf);
   void ComputeVelocityLumpedMass(Vector & row_sums);
   void ComputeVelocityLumpedMassByHand(const Vector &S, Vector & row_sums);
   void ComputeAverageMassAtGeo(const Vector &S, Vector &vec_weights);
   void SolveHiOp(const Vector &S, const Vector &S_old, const int &lm_option, const int &target_option, const double &t, const double &dt, ParGridFunction &mv_gf_l);
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

   // Enforce Mass Conservation
   void SetMassConservativeDensity(Vector &S) const;
   void ComputeDensity(const Vector &S, ParGridFunction &rho_gf) const;

   // Validate mass conservation
   void ValidateMassConservation(const Vector &S, ParGridFunction & mc_gf, double &mass_loss) const;
   void SetInitialMassesAndVolumes(const Vector &S);

   // Compute various time series data
   void ComputeMinDetJ(int &cell, double &minDetJ);
   bool IsMeshCollapsedGeom(const Vector &S);
   bool ComputeTimeSeriesData(const Vector &S, const double &t, const double &dt);

   // Various print functions
   void SaveStateVecsToFile(const Vector &S, const string &output_file_prefix, const string &output_file_suffix);
   void SaveTimeSeriesArraysToFile(const string &output_file_prefix, const string &output_file_suffix);

   // Kidder specific function
   void ComputeKidderAvgIntExtRadii(const Vector &S, double &avg_rad_int, double &avg_rad_ext);
   void ComputeKidderAvgDensityAndEntropy(const Vector &S, double &avg_density, double &avg_entropy);

   // Debugging
   mutable bool l2_dof_x_set = false;
   mutable ParGridFunction _l2x_gf;
   void SetL2DofX() const
   {
      pmesh->GetNodes(_l2x_gf);
      l2_dof_x_set = true;
   }
   void GetL2DofX(const int &dof, Vector &x) const;
};

class HydroODESolver : public ODESolver
{
protected:
   LagrangianLOOperator *hydro_oper;
   BlockVector dS_dt, S0;
public:
   HydroODESolver() : hydro_oper(NULL) { }
   virtual void Init(TimeDependentOperator&);
   virtual void Step(Vector&, double&, double&)
   { MFEM_ABORT("Time stepping is undefined."); }
};

class RK2AvgSolver : public HydroODESolver
{
public:
   RK2AvgSolver() { }
   virtual void Init(TimeDependentOperator &_f);
   virtual void Step(Vector &S, double &t, double &dt);
};

class HydroRK3SSPSolver : public HydroODESolver
{
protected:
   BlockVector k, y;
public:
   HydroRK3SSPSolver() { }
   virtual void Init(TimeDependentOperator &_f);
   virtual void Step(Vector &S, double &t, double &dt);
};

} // end ns hydroLO

} // end ns mfem

#endif // LAGLOS_SOLVER