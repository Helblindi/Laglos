#ifndef LAGLOS_SOLVER
#define LAGLOS_SOLVER

#include "mfem.hpp"
#include "laglos_assembly.hpp"
#include "initial_vals.hpp"
#include <iostream>
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

//
template <int dim, int problem>
class LagrangianLOOperator
{
protected:
   ParFiniteElementSpace &H1, &L2, &L2V;
   ParMesh *pmesh;
   ParLinearForm *m_lf;
   HypreParVector *m_hpv;

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

   Array<int> block_offsets;
   Array<int> BdrElementIndexingArray;
   Array<int> BdrVertexIndexingArray;  // Array to identify boundary vertices
   Array<double> vstar_arr;  // Array to store vstar for each face
   mutable Vector v_face_intermediate; // (5.7b)
   // const Array<int> &ess_tdofs;
   const int NE, l2dofs_cnt, h1dofs_cnt, l2vdofs_cnt;

   const int num_faces, num_vertices;

   double CFL;
   double timestep = 0.001;
   double timestep_first = 0.; // Set and used for activation function when prescribing left wall dirichlet BCs for Saltzman problem

   bool use_viscosity;
   bool mm;

public:
   LagrangianLOOperator(ParFiniteElementSpace &h1,
                        ParFiniteElementSpace &l2,
                        ParFiniteElementSpace &l2v,
                        ParLinearForm *m,
                        bool use_viscosity,
                        bool mm,
                        double CFL);
   ~LagrangianLOOperator();

   double GetCFL();

   void SetCFL(const double &_CFL); // STOPPED HERE.

   void IterateOverCells(); // TODO: Delete

   void CreateBdrElementIndexingArray();

   void CreateBdrVertexIndexingArray();

   bool IsBdrVertex(const int & node);

   void CreateVStarArr(const Vector &S);

   void MakeTimeStep(Vector &S, double & t, double dt);

   void ComputeStateUpdate(Vector &S_new, const double &t, const double dt);

   double GetTimeStepEstimate(const Vector &S);

   void CalculateTimestep(const Vector &S);

   double GetTimestep();

   void GetCellStateVector(const Vector &S, const int cell, Vector &U);

   void SetCellStateVector(Vector &S_new, const int cell, const Vector &U);

   /* cij comp */
   void CalcOutwardNormalInt(const Vector &S, const int cell, const int face, Vector & res);

   void Orthogonal(Vector &v);

   Vector GetIntDerRefShapeFunctions();

   /* RP Functions */
   double cosineSimilarity(const Vector &v1, const Vector &v2);

   /* Mesh movement */
   // void form_intermediate_velocity(const Vector &S, const double dt);

   void get_intermediate_face_velocity(const int & face, Vector & vel);

   void tensor(const Vector & v1, const Vector & v2, DenseMatrix & dm);

   // tests from 2023-01-05
   void compute_node_velocity_cwa(Vector &S, const double & t, const double & dt);
   void compute_node_velocity_RP(Vector &S, const double & t, const double & dt);
   void compute_node_velocity_LS(const Vector &S, 
                                 const Table &vertex_edge,  
                                 const int &vertex,
                                 const double & t, 
                                 const double & dt, 
                                 Vector &vertex_v,
                                 DenseMatrix &C,
                                 Vector &D,
                                 const string flag="NA", 
                                 void (*test_vel)(const Vector&, const double&, Vector&) = NULL);
   
   // Functions representing development on April 2023
   void compute_A(const DenseMatrix & C, const double d, const double &dt, DenseMatrix &A);
   void compute_B(const DenseMatrix &C, const Vector & D, const double d, const double &dt, Vector &B);
   void compute_determinant(const DenseMatrix &C, const double &dt, double & d);
   void compute_corrected_node_velocity(const DenseMatrix &C, const Vector &D, const double & dt, const Vector &vertex_x, Vector &vertex_v, const string flag="NA", void (*test_vel)(const Vector&, const double&, Vector&) = NULL);

   void compute_node_velocities(Vector &S, const double & t, const double & dt, const string flag="NA", void (*test_vel)(const Vector&, const double&, Vector&) = NULL);
   
   void compute_corrective_face_velocities(Vector &S, const double & t, const double & dt, const string flag="NA", void (*test_vel)(const Vector&, const double&, Vector&) = NULL);

   void fill_center_velocities_with_average(Vector &S, const string flag="NA", void (*test_vel)(const Vector&, const double&, Vector&) = NULL);
   void fill_face_velocities_with_average(Vector &S, const string flag="NA", void (*test_vel)(const Vector&, const double&, Vector&) = NULL);

   void update_node_velocity(Vector &S, const int & node, const Vector & vel);

   void get_node_velocity(const Vector &S, const int & node, Vector & vel);

   void get_node_position(const Vector &S, const int & node, Vector & x);

   void EnforceBoundaryConditions(Vector &S);

   void MoveMesh(Vector &S, GridFunction & x_gf, GridFunction & mv_gf_new, const double & t, const double & dt);

   void CheckMassConservation(const Vector &S);

   double CheckMassLoss(const Vector &S);

   /* Problem Description Functions */
   static double internal_energy(const Vector &U);

   static double specific_internal_energy(const Vector &U);
   static double pressure(const Vector &U);

   static double compute_lambda_max(const Vector & U_i,
                                    const Vector & U_j,
                                    const Vector & n_ij,
                                    const string flag="NA");
   
   static Vector velocity(const Vector & U);

   DenseMatrix flux(const Vector &U);

   // double ComputeViscosity(const Vector & UL, const Vector & UR);

};

} // end ns hydrodynamics

} // end ns mfem

#endif // LAGLOS_SOLVER