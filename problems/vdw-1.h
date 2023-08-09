/*********************************************************
* This file serves as a template for setting up a problem.
* Steps to implement a problem:
*  1) Copy this file and rename in $PROJECT_SRC/include/ to 
*     correspond to your problem
*  2) Modify problem specific constants a, b, gamma
*  3) Override pressure function
*  4) Specify initial conditions given by rho0, v0, and sie0
*  5) Add corresponding problem option to laglos.cpp #TODO WIP
*********************************************************/


#include "mfem.hpp"
#include "problem_base.h"
#include <cmath>

using namespace mfem;
using namespace std;

namespace mfem
{
namespace hydrodynamics
{

template<int dim>
class VdwTest1: public ProblemBase<dim>
{
public:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double a = 1.;     
   double b = 1.;
   double gamma = 1.02;
   bool distort_mesh = false;
   bool known_exact_solution = false;

   /* Override getters */
   virtual double get_a() override { return a; }
   virtual double get_b() override { return b; }
   virtual double get_gamma() override { return gamma; }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   virtual double pressure(const Vector &U) override
   {
       // Use van der Waals
      double rho = 1. / U[0];
      double sie = this->specific_internal_energy(U);

      double val = (this->get_gamma() - 1.) * (rho * sie + this->get_a() * pow(rho, 2)) / (1. - this->get_b() * rho) - this->get_a() * pow(rho,2);
      return val;
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   virtual double rho0(const Vector &x, const double & t) override
   {
      if (t < 1.e-16) {
         if (x[0] <= 0.)
         {
            return 0.10;
         }
         else
         {
            assert(x[0] <= 1.);
            return 0.39;
         }
      }
      else {
         return 0.5; // TODO: Exact representation of sie0
      }
   }
   virtual void v0(const Vector &x, const double & t, Vector &v) override
   {
      if (t < 1.e-16) {
         if (x[0] <= 0.)
         {
            v[0] = -0.475504638574729;
            return;
         }
         else
         {
            assert(x[0]<=1.);
            v[0] = -0.121375781741349;
            return;
         }
      }
      else {
         v[0] = 14.;
         return; 
      }
   }
   virtual double sie0(const Vector &x, const double & t) override
   {
      if (t < 1.e-16) {
         if (x[0] <= 0.)
         {
            return 14.337916411885988;
         }
         else
         {
            assert(x[0]<=1.);
            return 14.560722040683306;
         }
      }
      else {
         return 14.; // TODO: Exact representation of sie0
      }
   }

}; // End class

} // ns hydrodynamics
} // ns mfem