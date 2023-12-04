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
class VdwTest4: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 1., _b = 1., _gamma = 1.02;
   bool distort_mesh = false;
   bool known_exact_solution = false;
   bool _bcs = false;

public:
   VdwTest4()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_bcs_indicator(_bcs);
   }

   /* Override getters */
   virtual bool get_distort_mesh() override { return distort_mesh; }
   virtual bool has_exact_solution() override { return known_exact_solution; }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   virtual double pressure(const Vector &U, const int &cell_attr=0) override
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
            return 0.9932;
         }
         else
         {
            assert(x[0] <= 1.);
            return 0.9500;
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
            v[0] = 3.;
            return;
         }
         else
         {
            v[0] = -3.;
            return;
         }
      }
      else {
         v = 0.;
         return;
      }
   }
   virtual double sie0(const Vector &x, const double & t) override
   {
      if (t < 1.e-16) {
         if (x[0] <= 0.)
         {
            return 0.029143658477667977;
         }
         else
         {
            return 6.688157894736825;
         }
      }
      else {
         return .5; // TODO: Exact representation of sie0
      }
   }

}; // End class

} // ns hydrodynamics
} // ns mfem