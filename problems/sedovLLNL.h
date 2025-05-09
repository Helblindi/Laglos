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
#include "sedov_exact.hpp"
#include <cmath>
#include <string>

using namespace mfem;
using namespace std;

namespace mfem
{
namespace hydroLO
{

template<int dim>
class SedovLLNLProblem: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 1.4;
   bool _distort_mesh = false;
   bool _known_exact_solution = true;
   bool _thbcs = true; // Indicator for thermal boundary conditions
   bool _mvbcs = true; // Indicator for mv boundary conditions
   string _indicator = "Sedov"; // Possible: saltzmann

public:
   SedovLLNLProblem()
   {
      if (dim != 2)
      {
         MFEM_ABORT("Dimension != 2 is not supported.\n");
      }
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

   /* Override specific update functions */
   void lm_update(const double b_covolume) override 
   {
      this->set_b(b_covolume);
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double rho0(const Vector &x, const double & t) const override
   {
      double xyt[3];
      xyt[0] = x[0], xyt[1] = x[1], xyt[2] = t;

      double rho = sedov_rho(xyt);
      return rho;
   }
   void v0(const Vector &x, const double & t, Vector &v) const override
   {
      v = 0.;

      double xyt[3];
      xyt[0] = x[0], xyt[1] = x[1], xyt[2] = t;

      v[0] = sedov_vx(xyt);
      v[1] = sedov_vy(xyt);

      return;
   }
   double sie0(const Vector &x, const double & t) const override
   {
      if (t < 0.)
      {
         std::cout << "SedovLLNLProblem::sie0: t < 0.0\n";
         std::cout << "Should not be here.\n";
      }
      double xyt[3];
      xyt[0] = x[0], xyt[1] = x[1], xyt[2] = t;

      Vector v(dim);
      v0(x, t, v);
      double p = p0(x,t);
      double rho = rho0(x,t);
      double e = p / ((_gamma - 1.) * rho);
      double se = sedov_e(xyt);
      // e and se and the same, no bug in sedov_exact
      // if (x.Norml2() < .1)
      // {
      //    cout << "pressure: " << p << ", density: " << rho << ", sie: " << e << ", se: " << se << endl;
      // }
      
      // return sedov_e(xyt);
      return e;
   }
   double p0(const Vector &x, const double &t) const override
   {
      double xyt[3];
      xyt[0] = x[0], xyt[1] = x[1], xyt[2] = t;
      return sedov_p(xyt);
   }

   /*********************************************************
    * Sedov specific functions
    *********************************************************/
}; // End class

} // ns hydroLO
} // ns mfem