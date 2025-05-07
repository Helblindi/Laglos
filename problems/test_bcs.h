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

template<int dim>
class TestBCs: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 7/5.;
   bool _distort_mesh = false;
   bool _known_exact_solution = false;
   bool _bcs = true; // Indicator for boundary conditions
   string _indicator = "TestBCs"; // Possible: saltzmann

   double rhoL = 1., rhoR = 0.125;
   double pL = 1., pR = .1;
   double vL = 0., vR = 0.;
   double x_center = 0.5;

   // CFL change
   bool _change_cfl = false;

public:
   TestBCs()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_bcs_indicator(_bcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);

      // Set Equation of state
      this->eos = std::unique_ptr<EquationOfState>(new IdealGasEOS());
   }

   /* Override specific update functions */
   void lm_update(const double b_covolume) override 
   {
      this->set_b(b_covolume);
   }

   double p0(const Vector &x, const double & t) const override
   {
      if (t < 1e-12)
      {
         return (x(0) < x_center) ? pL : pR;
      }
      assert(false);
      return 0.1;
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double rho0(const Vector &x, const double & t) const override
   {
      if (x[0] < x_center)
      {
         return rhoL;
      }
      else { return rhoR; }
   }
   void v0(const Vector &x, const double & t, Vector &v) const override
   {
      v = 0.;
      if (x[0] < x_center)
      {
         v[0] = vL;
      }
      else { v[0] = vR; }
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