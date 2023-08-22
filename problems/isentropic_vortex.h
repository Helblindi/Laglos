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
public:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double a = 0.;     
   double b = 0.;
   double gamma = 7./5.;
   bool distort_mesh = false;
   bool known_exact_solution = true;
   bool bcs = true;
   string indicator = "IsentropicVortex"; // Possible: saltzman

   // Free steam conditions
   const double rho_inf = 1., p_inf = 1., T_inf = 1.;
   // Vector u_inf(2);
   // // u_inf[0] = 2.; u_inf[1] = 0.;
   const double beta = 5.;
   double xc_0 = 0., xc_1 = 0.;

   /* Override getters */
   virtual double get_a() override { return a; }
   virtual double get_b() override { return b; }
   virtual double get_gamma() override { return gamma; }
   virtual string get_indicator() override { return indicator; }
   virtual bool get_distort_mesh() override { return distort_mesh; }
   virtual bool has_exact_solution() override { return known_exact_solution; }
   virtual bool has_boundary_conditions() override { return bcs; }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   virtual double pressure(const Vector &U) override
   {
      cout << "IV::Pressure\n";
      double density = 1. / U[0];
      return pow(density, this->get_gamma());
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   virtual double pressure(const Vector &x, const double &t)
   {
      cout << "IV::Pressure\n";
      double density = rho0(x,t);
      return pow(density, this->get_gamma());
   }
   virtual double rho0(const Vector &x, const double & t) override
   {
      cout << "IV::rho\n";
      Vector center(2);
      center[0] = xc_0 + 2. * t;
      center[1] = xc_1;

      Vector x_bar = x;
      x_bar -= center;
      const double r = x_bar.Norml2();

      const double dT = -1. * (exp(1 - pow(r,2)) * pow(beta, 2) * (this->get_gamma() - 1.)) / (8. * this->get_gamma() * pow(M_PI, 2));
      const double rho = pow(T_inf + dT, 1./(this->get_gamma() - 1.));
      cout << "rho: " << rho << endl;

      return rho;
   }
   virtual void v0(const Vector &x, const double & t, Vector &v) override
   {
      cout << "IV::v0\n";
      Vector center(2);
      center[0] = xc_0 + 2. * t;
      center[1] = xc_1;

      Vector x_bar = x;
      x_bar -= center;
      const double r = x_bar.Norml2();

      const double coeff = (exp((1. - pow(r,2)) / 2.) * beta) / (2. * M_PI);

      // v[0] = u_inf[0] + coeff * x_bar[1] * -1.;
      // v[1] = u_inf[1] + coeff * x_bar[0];
      v[0] = 2. + coeff * x_bar[1] * -1.;
      v[1] = 0. + coeff * x_bar[0];
      return;
   }
   virtual double sie0(const Vector &x, const double & t) override
   {
      cout << "IV::sie0\n";
      return pressure(x,t) / this->rho0(x, t) / (this->get_gamma() - 1.0);
   }

}; // End class

} // ns hydrodynamics
} // ns mfem