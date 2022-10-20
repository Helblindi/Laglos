#ifndef LAGLOS_SOLVER
#define LAGLOS_SOLVER

#include "mfem.hpp"
#include "laglos_assembly.hpp"
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
template <int dim>
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

   const int Vsize_L2;
   const int TVSize_L2;
   const HYPRE_Int GTVSize_L2;
   const int NDofs_L2;

   const int Vsize_L2V;
   const int TVSize_L2V;
   const HYPRE_Int GTVSize_L2V;
   const int NDofs_L2V;

   Array<int> block_offsets;
   // Reference to the current mesh configuration.
   mutable ParGridFunction x_gf;
   // const Array<int> &ess_tdofs;
   const int NE, l2dofs_cnt, h1dofs_cnt, l2vdofs_cnt;

   double timestep = 0.001;

public:
   static constexpr double gamma = 7./5.;
   LagrangianLOOperator(ParFiniteElementSpace &h1,
                        ParFiniteElementSpace &l2,
                        ParFiniteElementSpace &l2v,
                        ParLinearForm *m);
   ~LagrangianLOOperator();

   void MakeTimeStep(Vector &S, double & t, double dt);

   double GetTimeStepEstimate(const Vector &S) const;

   void CalculateTimestep(const Vector &S);

   double GetTimestep();

   void GetCellStateVector(const Vector &S, const int cell, Vector &U);

   void SetCellStateVector(Vector &S_new, const int cell, const Vector &U);

   /* cij comp */
   void build_C(const Vector &S, double dt);

   void CalcOutwardNormalInt(const Vector &S, const int cell, const int face, Vector & res);

   void Orthogonal(Vector &v);

   Vector GetIntDerRefShapeFunctions();

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