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
class SmoothWave: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 7/5.;
   bool distort_mesh = false;
   bool known_exact_solution = true;

public:
   SmoothWave()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
   }

   /* Override getters */
   bool get_distort_mesh() override { return distort_mesh; }
   bool has_exact_solution() override { return known_exact_solution; }

   /* Override specific update functions */
   void lm_update(const double b_covolume) override 
   {
      this->set_b(b_covolume);
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double pressure(const Vector &U, const int &cell_attr=0) override
   {
      // pL = pR = 1.
      return (this->get_gamma() - 1.) * this->internal_energy(U);
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double rho0(const Vector &x, const double & t) override
   {
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
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      v[0] = 1.;
      return;
   }
   double sie0(const Vector &x, const double & t) override
   {
      return 1.0 / this->rho0(x, t) / (this->get_gamma() - 1.0);
   }

}; // End class

} // ns hydrodynamics
} // ns mfem