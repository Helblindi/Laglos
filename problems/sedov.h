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
class SedovProblem: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 7/5.;
   bool distort_mesh = false;
   bool known_exact_solution = false;
   bool bcs = false; // Indicator for boundary conditions
   string indicator = ""; // Possible: saltzmann

   // CFL change
   bool _change_cfl = false;
   // constexpr static double CFL_first = 0.5;
   // constexpr static double CFL_second = 0.5;
   // constexpr static double CFL_time_change = 0.01;

public:
   SedovProblem()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
   }

   /* Override getters */
   bool get_distort_mesh() override { return distort_mesh; }
   bool has_exact_solution() override { return known_exact_solution; }

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
      if (t < 1.e-16)
      {
         // Initial condition 
         return 1.;
      }
      else 
      {
         return 0.;
      }
   }
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      v = 0.;

      if (t < 1.e-16)
      {
         // Initial condition
         return;
      }
      return;
   }
   double sie0(const Vector &x, const double & t) override
   {
      double norm = x.Norml2();
      if (t < 1.e-16)
      {
         // Initial condition
         if (norm <= 0.05) {
            return 0.979264;
         }
         else
         {
            return 0.;
         }
      }
      return 0.;
   }

}; // End class

} // ns hydrodynamics
} // ns mfem