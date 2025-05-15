/*********************************************************
* This file serves as a template for setting up a problem.
* Steps to implement a problem:
*  1) Copy this file and rename in $PROJECT_SRC/include/ to 
*     correspond to your problem
*  2) Modify problem specific constants a, b, gamma
*  3) Override pressure function
*  4) Specify initial conditions given by rho0, v0, and sie0
*  5) Add corresponding problem option to laglos.cpp
*  6) Add to include header file problems/test_problems_include.h
*  7) Enforce boundary conditions in laglos_solver.cpp
*     a) This may require edits to the following functions 
*        depending on the type of BCs that will be 
*        implemented
*        - LagrangianLOOperator::CreateBdrElementIndexingArray()
*        - LagrangianLOOperator::CreateBdrVertexIndexingArray()
*        - LagrangianLOOperator::FillCellBdrFlag()
*        - LagrangianLOOperator::UpdateMeshVelocityBCs()
*        - LagrangianLOOperator::EnforceL2BC()
*        - ProblemBase::get_mv_bcs_need_updating()
*        - ProblemBase::update_additional_BCs()
*     b) To note about boundary conditions in the mesh file:
*        1: vx = 0
*        2: vy = 0
*        3: vz = 0
*        4: vr = 0 (radial boundary condition)
*        5: Misc (possibly dirichlet condition)
*  8) Edit SaveStateVecsToFile if the problem is a radial one.
*********************************************************/


#include "mfem.hpp"
#include "problem_base.h"
#include "riemann1D.hpp"
#include <cmath>

using namespace mfem;
using namespace std;

namespace mfem
{
namespace hydroLO
{

class ElasticTwist: public ProblemBase
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _gamma = 3.4;
   bool _known_exact_solution = false;
   bool _thbcs = false; // Indicator for thermal boundary conditions
   bool _mvbcs = false; // Indicator for mv boundary conditions
   string _indicator = "ElasticTwist";

   /* Material parameters */
   double rho = 2.7E3, p = 1.E5;
   const double p_inf = 2.15E10;
   const double _mu = 2.6E10;

   /* Problem specific parameters, can be tweaked */
   const double omega = 40000; // revolutions per second
   const double R = 0.05;      // Internal radius
   const double R2 = R*R;

public:
   ElasticTwist(const int &_dim) : ProblemBase(_dim)
   {
      if (dim != 2)
      {
         /*
         I do not know how to implement 2d velocity for a 1d problem,
         hence we must be compiled in 2D.
         */
         MFEM_ABORT("Dimension must be 2 for the elastic twist problem.\n");
      }
      this->set_gamma(_gamma);
      this->set_pinf(p_inf);
      this->set_shear_modulus(_mu);
      this->set_indicator(_indicator);
      this->set_thbcs_indicator(_thbcs);
      this->set_mvbcs_indicator(_mvbcs);
      this->set_exact_solution(_known_exact_solution);

      // Set Equation of state
      this->eos = std::unique_ptr<EquationOfState>(new NobleAbelStiffenedGasEOS(this->b,this->q,this->p_inf));
   }
   
   /* Override specific update functions */
   void lm_update(const double b_covolume) override 
   {
      this->set_b(b_covolume);
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double p0(const Vector &x, const double & t) const override
   {
      return p;
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double rho0(const Vector &x, const double & t) const override
   {
      if (t < 1e-12)
      {
         return rho;
      }
      else
      {
         MFEM_ABORT("No exact solution.\n");
         return -1.;
      }
   }
   void v0(const Vector &x, const double & t, Vector &v) const override
   {
      v = 0.;
      if (t < 1e-12)
      {
         double val = x[0] * x[0] + x[1] * x[1];
         if (val < R2)
         {
            v[0] = -1. * omega * x[1];
            v[1] = omega * x[0];
         }
      }

      return;
   }

   double sie0(const Vector &x, const double & t) const override
   {
      double _rho = rho0(x,t);
      double _p = p0(x,t);
      return this->eos->energy(_p, _rho, this->get_gamma());
   }
}; // End class

} // ns hydroLO
} // ns mfem
