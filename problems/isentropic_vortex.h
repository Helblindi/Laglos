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
class IsentropicVortex: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0.;     
   double _b = 0.;
   double _gamma = 7./5.;
   bool distort_mesh = false;
   bool known_exact_solution = true;
   bool bcs = true;
   string _indicator = "IsentropicVortex"; // Possible: saltzmann

   // Free steam conditions
   const double rho_inf = 1., p_inf = 1., T_inf = 1.;
   const double beta = 5.;
   double xc_0 = 0., xc_1 = 0.;
   double vc_0 = 0., vc_1 = 0.; // u_inf

public:
   IsentropicVortex()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
   }

   /* Override getters */
   bool get_distort_mesh() override { return distort_mesh; }
   bool has_exact_solution() override { return known_exact_solution; }
   bool has_boundary_conditions() override { return bcs; }

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
      double density = 1. / U[0];
      return pow(density, this->get_gamma());
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double pressure(const Vector &x, const double &t)
   {
      double density = rho0(x,t);
      return pow(density, this->get_gamma());
   }
   double rho0(const Vector &x, const double & t) override
   {
      Vector center(2);
      center[0] = xc_0 + vc_0 * t;
      center[1] = xc_1 + vc_1 * t;

      Vector x_bar = x;
      x_bar -= center;
      const double r = x_bar.Norml2();

      const double dT = -1. * (exp(1 - pow(r,2)) * pow(beta, 2) * (this->get_gamma() - 1.)) / (8. * this->get_gamma() * pow(M_PI, 2));
      const double rho = pow(T_inf + dT, 1./(this->get_gamma() - 1.));

      return rho;
   }
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      Vector center(2);
      center[0] = xc_0 + vc_0 * t;
      center[1] = xc_1 + vc_1 * t;

      Vector x_bar = x;
      x_bar -= center;
      const double r = x_bar.Norml2();

      const double coeff = (exp((1. - pow(r,2)) / 2.) * beta) / (2. * M_PI);

      v[0] = vc_0 - coeff * x_bar[1];
      v[1] = vc_1 + coeff * x_bar[0];
      
      return;
   }
   double sie0(const Vector &x, const double & t) override
   {
      return pressure(x,t) / this->rho0(x, t) / (this->get_gamma() - 1.0);
   }

}; // End class

} // ns hydrodynamics
} // ns mfem