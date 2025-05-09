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
#include <cmath>
#include <string>

using namespace mfem;
using namespace std;

namespace mfem
{
namespace hydroLO
{

class TaylorCoefficient : public Coefficient
{
public:
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
   {
      Vector x(2);
      T.Transform(ip, x);
      return 3.0 / 8.0 * M_PI * ( cos(3.0*M_PI*x(0)) * cos(M_PI*x(1)) -
                                    cos(M_PI*x(0))     * cos(3.0*M_PI*x(1)) );
   }
};

class TaylorGreenProblem: public ProblemBase
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 5./3.;
   bool _distort_mesh = false;
   bool _known_exact_solution = true;
   bool _thbcs = true; // Indicator for thermal boundary conditions
   bool _mvbcs = true; // Indicator for mv boundary conditions
   string _indicator = "TaylorGreen"; // Possible: saltzmann

public:
   TaylorGreenProblem(const int &_dim) : ProblemBase(_dim)
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_thbcs_indicator(_thbcs);
      this->set_mvbcs_indicator(_mvbcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);

      // Set Equation of state
      this->eos = std::unique_ptr<EquationOfState>(new IdealGasEOS());
   }

   /* Optionally overridden, or removed */
   double get_gamma(const int &cell_attr = 0) const override { return _gamma; }
   void lm_update(const double b_covolume) override {}
   void update(Vector vec, double t = 0.) override {}

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double p0(const Vector &x, const double & t) const override
   {
      double _rho = rho0(x,t);
      double _sie = sie0(x,t);
      return this->eos->pressure(_rho, _sie, this->get_gamma());
   }
   double rho0(const Vector &x, const double & t) const override
   {
      return 1.;
   }
   void v0(const Vector &x, const double & t, Vector &v) const override
   {
      v[0] = std::sin(M_PI * x[0]) * std::cos(M_PI * x[1]);
      v[1] = -1. * std::cos(M_PI * x[0]) * std::sin(M_PI * x[1]);
      return;
   }
   double sie0(const Vector &x, const double & t) const override
   {
      const double denom = 2.0 / 3.0;  // (5/3 - 1) * density.
      double val;
      if (x.Size() == 2)
      {
         val = 1.0 + (cos(2*M_PI*x(0)) + cos(2*M_PI*x(1))) / 4.0;
      }
      else
      {
         val = 100.0 + ((cos(2*M_PI*x(2)) + 2) *
                        (cos(2*M_PI*x(0)) + cos(2*M_PI*x(1))) - 2) / 16.0;
      }
      return val/denom;
   }

}; // End class

} // ns hydroLO
} // ns mfem