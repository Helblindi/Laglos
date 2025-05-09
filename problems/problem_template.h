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
class ProblemTemplate : public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 1.4;
   bool _distort_mesh = false;
   bool _known_exact_solution = false;
   bool _thbcs = false; // Indicator for thermal boundary conditions
   bool _mvbcs = false; // Indicator for mv boundary conditions
   string _indicator = ""; // Possible: saltzmann

   // CFL change, can remove if not needed
   bool _change_cfl = false;
   double _cfl_first = 0.5;
   double _cfl_second = 0.5;
   double _cfl_time_change = 0.01;

public:
   ProblemTemplate()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_thbcs_indicator(_thbcs);
      this->set_mvbcs_indicator(_mvbcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);
      // CFL change
      this->set_cfl_change(_change_cfl);
      this->set_cfl_first(_cfl_first);
      this->set_cfl_second(_cfl_second);
      this->set_cfl_time_change(_cfl_time_change);

      // Set Equation of state
      this->eos = std::unique_ptr<EquationOfState>(new IdealGasEOS()); // Options are IdealGasEOS, VanDerWaalsEOS, NobleAbelStiffenedGasEOS, and PolytropicEOS
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
      /*
      Must Override
      */
      return 0.;
   }
   double rho0(const Vector &x, const double & t) const override
   {
      /*
      Must Override
      */
      return 0.;
   }
   void v0(const Vector &x, const double & t, Vector &v) const override
   {
      /*
      Must Override
      */
      return;
   }
   double sie0(const Vector &x, const double & t) const override
   {
      /*
      Must Override
      */
      return 0.;
   }

}; // End class

} // ns hydroLO
} // ns mfem