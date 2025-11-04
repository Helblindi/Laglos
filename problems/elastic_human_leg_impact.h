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

class ElasticHumanLegImpact: public ProblemBase
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _gamma = 4.22;
   bool _known_exact_solution = false;
   bool _thbcs = false; // Indicator for thermal boundary conditions
   bool _mvbcs = false; // Indicator for mv boundary conditions
   string _indicator = "ElasticHumanLegImpact";

   //https://www.sciencedirect.com/science/article/pii/S0021999109002654?fr=RR-2&ref=pdf_download&rr=928faaf93aca69c5
   // 5.2 elastic projectile plate section 7.1
   double rho_bone = 1986., rho_flesh = 1000., rho_bar = 2400.; // kg/m^3
   double v_proj = -13.4112; // m/s, equivalent to 30 mph impact speed
   const double p_inf = 3.42E8; // Pa
   /* Different shear moduli for projectile plate */
   const double _mu = 9.2E8; // Pa


   /* helper function to determine the region based on x,y coords */
   bool is_bone_region(const Vector &x) const
   {
      double _x = x[0];
      double _y = x[1];
      if (_x >= 0. && _x <= 0.03 && _y >= 0.0 && _y <= 0.24) { return true; }
      else { return false; }
   }
   bool is_flesh_region(const Vector &x) const
   {
      double _x = x[0];
      double _y = x[1];
      if (_x >= 0.03 && _x <= 0.09 && _y >= 0.0 && _y <= 0.24) { return true; }
      else { return false; }
   }
   bool is_bar_region(const Vector &x) const
   {
      double _x = x[0];
      double _y = x[1];
      if (_x >= 0.09 && _x <= 0.19 && _y >= 0.07 && _y <= 0.17) { return true; }
      else { return false; }
   }

public:
   ElasticHumanLegImpact(const int &_dim) : ProblemBase(_dim)
   {
      if (dim != 2)
      {
         MFEM_ABORT("Dimension must be 2 for the elastic human leg impact problem.\n");
      }
      this->set_indicator(_indicator);
      this->set_pinf(p_inf);
      this->set_shear_modulus(_mu);
      this->set_thbcs_indicator(_thbcs);
      this->set_mvbcs_indicator(_mvbcs);
      this->set_exact_solution(_known_exact_solution);
      this->set_gamma(_gamma);

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
      if (t < 1e-12)
      {
         return 0.;
      }
      else
      {
         MFEM_ABORT("No exact solution.\n");
         return -1.;
      }
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double rho0(const Vector &x, const double & t) const override
   {
      assert(t < 1e-12); // No exact solution
      if (is_bone_region(x)) {
         return rho_bone;
      } else if (is_flesh_region(x)) {
         return rho_flesh;
      } else if (is_bar_region(x)) {
         return rho_bar;
      } else {
         MFEM_ABORT("Point not in any defined region.\n");
         return -1.;
      }
   }
   void v0(const Vector &x, const double & t, Vector &v) const override
   {
      assert(t < 1e-12);
      v = 0.;
      if (is_bar_region(x)) { v[0] = v_proj; }

      return;
   }

   double sie0(const Vector &x, const double & t) const override
   {
      double _rho = rho0(x,t);
      double _p = p0(x,t);
      return this->eos->energy(_p, _rho, _gamma);
   }
}; // End class

} // ns hydroLO
} // ns mfem
