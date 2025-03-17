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
class NohProblem: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 5./3.;
   bool _distort_mesh = false;
   bool _known_exact_solution = true;
   bool _bcs = false; // Indicator for boundary conditions
   string _indicator = "Noh"; // Possible: saltzmann

public:
   NohProblem()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_bcs_indicator(_bcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);
   }

   /* Override specific update functions */
   void lm_update(const double b_covolume) override 
   {
      this->set_b(b_covolume);
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double pressure(const double &rho, const double &sie, const int &cell_attr=0) override
   {
      return (this->get_gamma() - 1.) * rho * sie;
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double p0(const Vector &x, const double &t) override
   {
      return (_gamma - 1.) * rho0(x,t) * sie0(x,t);
   }

   double rho0(const Vector &x, const double & t) override
   {
      /* Initial condition */
      double norm = x.Norml2();

      if (t < 1.e-16) {
         return 1.;
      }
      /* Exact solution */
      else 
      {
         if (t < 3. * norm) {
            return 1.0 + t / norm;
         } else { // (t / 3. >= norm)
            return 16.0;
         }
      }
   }
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      double norm = x.Norml2();
      v = 0.;

      /* Initial condition */
      if (t < 1.e-16) {
         if (norm > 1.e-16) {
            v[0] = -x[0] / norm, v[1] = -x[1] / norm;
         }
      } // Check here

      /* Exact solution */
      else if (t < 3. * norm) {
         v[0] = -x[0] / norm, v[1] = -x[1] / norm;
      } 

      return;
   }
   double sie0(const Vector &x, const double & t) override
   {
      double norm = x.Norml2();

      /* Initial condition */
      if (t < 1.e-16) {
         // if (norm > 1.e-16) {
         return 1.e-12 / (this->get_gamma() - 1.);
         // }
      }

      /* Exact solution */
      else {
         if (t < 3. * norm) {
            return 1.e-12 / (this->get_gamma() - 1.);
         } 
         else  {
            // All KE converted to IE on the backside of the shock
            return 0.5 + 1.e-12 / (this->get_gamma() - 1.);
         }
      }
      cout << "norm: " << norm << ", t: " << t << "\n";
      MFEM_ABORT("Noh problem, this spot should not have been encountered.\n");
      return 0.;
   }

}; // End class

} // ns hydroLO
} // ns mfem