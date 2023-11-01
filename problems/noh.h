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
*********************************************************/


#include "mfem.hpp"
#include "problem_base.h"
#include <cmath>
#include <string>

using namespace mfem;
using namespace std;

namespace mfem
{
namespace hydrodynamics
{

template<int dim>
class NohProblem: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 5./3.;
   bool distort_mesh = false;
   bool known_exact_solution = true;
   bool bcs = false; // Indicator for boundary conditions
   string indicator = "Noh"; // Possible: saltzmann

   // CFL change
   bool _change_cfl = false;
   // constexpr static double CFL_first = 0.5;
   // constexpr static double CFL_second = 0.5;
   // constexpr static double CFL_time_change = 0.01;

public:
   NohProblem()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
   }

   /* Override getters */
   string get_indicator() override { return indicator; }
   bool get_distort_mesh() override { return distort_mesh; }
   bool has_exact_solution() override { return known_exact_solution; }
   bool has_boundary_conditions() override { return bcs; }

   bool change_cfl() override { return _change_cfl; }
   // double get_cfl_first() override { return CFL_first; }
   // double get_cfl_second() override { return CFL_second; }
   // double get_cfl_time_change() override { return CFL_time_change; }

   /* Override specific update functions */
   void lm_update(const double b_covolume) override 
   {
      this->set_b(b_covolume);
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double pressure(const Vector &U) override
   {
      return (this->get_gamma() - 1.) * this->internal_energy(U);
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
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

} // ns hydrodynamics
} // ns mfem