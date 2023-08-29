#ifndef LAGLOS_SOLVER
#define LAGLOS_SOLVER

#include "mfem.hpp"
// #include "initial_vals.hpp"
#include "problem_base.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>

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
   ParGridFunction v_CR_gf; // 5.7(b)
   ParGridFunction v_geo_gf; // 5.11
   ParMesh *pmesh;
   ParLinearForm *m_lf;
   HypreParVector *m_hpv;

   // Problem specific
   ProblemBase<dim> * pb;

   // Matrix to hold max wavespeed values dij
   SparseMatrix * dij_sparse;

   // FE spaces local and global sizes
   const int Vsize_H1;
   const int TVSize_H1;
   const HYPRE_Int GTVSize_H1;
   const int NDofs_H1;
   const int NVDofs_H1;

   const int Vsize_L2;
   const int TVSize_L2;
   const HYPRE_Int GTVSize_L2;
   const int NDofs_L2;

   const int Vsize_L2V;
   const int TVSize_L2V;
   const HYPRE_Int GTVSize_L2V;
   const int NDofs_L2V;

   IntegrationRule RT_ir;

   // Tables to relate cell to the contained faces
   Table element_face;
   Table * vertex_element;
   Table * face_element;
   Array<int> block_offsets;
   Array<int> BdrElementIndexingArray;
   Array<int> BdrVertexIndexingArray;  // Array to identify boundary vertices

   int el_num_faces;
   const int num_elements, num_vertices, num_faces, num_edges;

   double CFL;
   double timestep = 0.001;
   double timestep_first = 0.; // Set and used for activation function when prescribing left wall dirichlet BCs for Saltzman problem

   bool use_viscosity;
   bool mm;

public:
   enum DofEntity {corner, face, cell};

   LagrangianLOOperator(ParFiniteElementSpace &h1,
                        ParFiniteElementSpace &l2,
                        ParFiniteElementSpace &l2v,
                        ParFiniteElementSpace &cr,
                        ParLinearForm *m,
                        ProblemBase<dim> *_pb,
                        bool use_viscosity,
                        bool mm,
                        double CFL);
   ~LagrangianLOOperator();

   double GetCFL() { return this->CFL; }
   double GetTimestep() { return timestep; }
   void SetCFL(const double &_CFL) { this->CFL = _CFL; }

   void GetEntityDof(const int GDof, DofEntity & entity, int & EDof);

   void CreateBdrElementIndexingArray();

   void CreateBdrVertexIndexingArray();

   bool IsBdrVertex(const int & node) { return (BdrVertexIndexingArray[node] == 1); }

   void MakeTimeStep(Vector &S, const double & t, double & dt);

   void ComputeStateUpdate(Vector &S_new, const double &t, const double dt);

   void InitializeDijMatrix();
   void BuildDijMatrix(const Vector &S);
   void CalculateTimestep(const Vector &S);

   void GetCellStateVector(const Vector &S, const int cell, Vector &U);

   void SetCellStateVector(Vector &S_new, const int cell, const Vector &U);

   /* cij comp */
   void CalcOutwardNormalInt(const Vector &S, const int cell, const int face, Vector & res);

   void Orthogonal(Vector &v);

   /* Mesh movement */
   void ComputeIntermediateFaceVelocities(const Vector &S, 
                                             const double t,
                                             const string ="NA", 
                                             void (*test_vel)(const Vector&, const double&, Vector&) = NULL);
   void GetIntermediateFaceVelocity(const int & face, Vector & vel);

   void ComputeNodeVelocityRT(const int & node, double & dt, Vector &node_v, bool &is_dt_changed);
   void IntGradRT(const int cell, DenseMatrix & res);
   void ComputeGeoV();
   void GetViGeo(const int & node, Vector & vel);
   void ComputeCiGeo(const int & node, DenseMatrix & res);
   void GetVCRgf(ParGridFunction & _v_CR_gf) { _v_CR_gf = this->v_CR_gf; }
   void GetVGeogf(ParGridFunction & _v_geo_gf) { _v_geo_gf = this->v_geo_gf; }
   
   void ComputeDeterminant(const DenseMatrix &C, const double &dt, double & d);
   void ComputeNodeVelocities(Vector &S, const double & t, double & dt, const string ="NA", void (*test_vel)(const Vector&, const double&, Vector&) = NULL);
   void ComputeCorrectiveFaceVelocities(Vector &S, const double & t, const double & dt, const string ="NA", void (*test_vel)(const Vector&, const double&, Vector&) = NULL);

   void FillCenterVelocitiesWithAvg(Vector &S);
   void FillFaceVelocitiesWithAvg(Vector &S, const string ="NA", void (*test_vel)(const Vector&, const double&, Vector&) = NULL);

   void ComputeMeshVelocities(Vector &S, const double & t, double & dt, const string ="NA", 
                              void (*test_vel)(const Vector&, const double&, Vector&) = NULL);
   void UpdateNodeVelocity(Vector &S, const int & node, const Vector & vel);
   void GetNodeVelocity(const Vector &S, const int & node, Vector & vel);
   void GetNodePosition(const Vector &S, const int & node, Vector & x);

   void EnforceExactBCOnCell(const Vector &S, const int & cell, const double &t, 
                             const double &dt, Vector & state_val);
   
   void EnforceMVBoundaryConditions(Vector &S, const double &t, const double &dt);
   
   double CalcMassLoss(const Vector &S);
   void CheckMassConservation(const Vector &S, ParGridFunction & mc_gf);

   // Various print functions
   void SaveStateVecsToFile(const Vector &S, const string &output_file_prefix, const string &output_file_suffix);
};

} // end ns hydrodynamics

} // end ns mfem

#endif // LAGLOS_SOLVER