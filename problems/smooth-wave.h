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
namespace hydroLO
{

class SmoothWave: public ProblemBase
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 7/5.;
   bool _distort_mesh = false;
   bool _known_exact_solution = true;
   bool _thbcs = false;
   bool _mvbcs = false; /* For use with extrapolated segment with indicator 99 */
   string _indicator = "SmoothWave";

public:
   SmoothWave(const int &_dim) : ProblemBase(_dim)
   {
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

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double p0(const Vector &x, const double &t) const override
   {
      // pL = pR = 1
      return 1.;
   }

   double rho0(const Vector &x, const double & t) const override
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
   void v0(const Vector &x, const double & t, Vector &v) const override
   {
      v[0] = 1.;
      return;
   }
   double sie0(const Vector &x, const double & t) const override
   {
      double _p = p0(x,t);
      double _rho = rho0(x,t);
      return this->eos->energy(_p, _rho, this->get_gamma());
   }

}; // End class

} // ns hydroLO
} // ns mfem