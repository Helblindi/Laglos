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
class ProblemTemplate: public ProblemBase<dim>
{
public:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double a = 0., b = 0., gamma = 7/5.;
   bool distort_mesh = false;
   bool known_exact_solution = false;
   bool bcs = false; // Indicator for boundary conditions
   string indicator = ""; // Possible: saltzmann

   /* Override getters */
   virtual double get_a() override { return a; }
   virtual double get_b() override { return b; }
   virtual double get_gamma() override { return gamma; }
   virtual bool get_distort_mesh() override { return distort_mesh; }
   virtual bool has_exact_solution() override { return known_exact_solution; }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   virtual double pressure(const Vector &U) override
   {
      /*
      Must Override
      */
      return 0.;
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   virtual double rho0(const Vector &x, const double & t) override
   {
      /*
      Must Override
      */
      return 0.;
   }
   virtual void v0(const Vector &x, const double & t, Vector &v) override
   {
      /*
      Must Override
      */
      return;
   }
   virtual double sie0(const Vector &x, const double & t) override
   {
      /*
      Must Override
      */
      return 0.;
   }

}; // End class

} // ns hydrodynamics
} // ns mfem