/*********************************************************
* The one-dimensional smooth wave test case comes from 
* section 5.2 of the following paper:
*
* Authors: Guermond, Nazarov, Popov, Tomas
* Title: SECOND-ORDER INVARIANT DOMAIN PRESERVING
*        APPROXIMATION OF THE EULER EQUATIONS
*        USING CONVEX LIMITING
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
class Smooth1D: public ProblemBase<dim>
{
public:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double a = 0.;     
   double b = 0.;
   double gamma = 7/5.;
   bool distort_mesh = false;
   bool known_exact_solution = true;

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
      // pL = pR = 1.
      return (this->get_gamma() - 1.) * this->internal_energy(U);
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   virtual double rho0(const Vector &x, const double & t) override
   {
      assert(dim == 1);
      double x0 = 0.1, x1 = 0.3;
      if (x0 <= x[0] - t && x[0] - t < x1)
      {
         double _rho = 1. + pow(2.,6) * (1. / pow(x1-x0,6)) * pow(x[0]-t-x0,3) * pow(x1-x[0]+t,3);
         return _rho;
      }
      else 
      {
         return 1.;
      }
   }
   virtual void v0(const Vector &x, const double & t, Vector &v) override
   {
      assert(dim==1);
      v = 1.;
      return;
   }
   virtual double sie0(const Vector &x, const double & t) override
   {
      return 1.0 / this->rho0(x, t) / (this->get_gamma() - 1.0);
   }

}; // End class

} // ns hydrodynamics
} // ns mfem