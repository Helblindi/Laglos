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
   bool known_exact_solution = true; // exact is known, information in supplementary material

   double rhoL = 0.1, rhoR = 0.39;
   double vL = -0.475504638574729, vR = -0.121375781741349;
   double pL = 0.022084258693080, pR = 0.039073167077590;
   double sieL = 14.337916411885988, sieR = 14.560722040683306;

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
            return rhoL;
         }
         else
         {
            assert(x[0] <= 1.);
            return rhoR;
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
            v[0] = vL;
            return;
         }
         else
         {
            assert(x[0]<=1.);
            v[0] = vR;
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
            return sieL;
         }
         else
         {
            assert(x[0]<=1.);
            return sieR;
         }
      }
      else {
         return 14.; // TODO: Exact representation of sie0
      }
   }

private:
   double rk_rho(const double & _rho, const double & _S)
   {
      double num = _S * this->get_gamma() * pow(_rho, this->get_gamma() - 1.);
      num *= (this->get_gamma() + 1.) * pow(1. - this->get_b() * _rho, -2. - this->get_gamma());
      num -= 6. * this->get_a() * _rho;

      double denom = pow(1 - this->get_b() * _rho, -1. - this->get_gamma());
      denom *= (_S * this->get_gamma() * pow(_rho, this->get_gamma()) - 2 * this->get_a() * pow(_rho,2) * pow(1.-this->get_b() * _rho, this->get_gamma() + 1.)) * _rho;

      return 2. * sqrt(denom) / num;
   }

   void sol_rho(double & _rhoz, 
                double & _Sz, 
                double (*phi)(const double &, const double &), // &rk_rho
                Vector & _xx,
                Vector & _rho)
   {
      double dx = _xx[0] - 0.;
      double k1 = phi(_rhoz, _Sz);
      double k2 = phi(_rhoz+dx*k1/2.,_Sz);
      double k3 = phi(_rhoz+dx*k2/2.,_Sz);
      double k4 = phi(_rhoz+dx*k3,_Sz);
      _rho[0] = _rhoz + (dx/6.)*(k1+2.*k2+2.*k3+k4);

      for (int i = 0; i < _xx.Size() - 1; i++)
      {
         dx = _xx[i+1] - _xx[i];
         k1 = phi(_rho[i], _Sz);
         k2 = phi(_rho[i]+dx*k1/2.,_Sz);
         k3 = phi(_rho[i]+dx*k2/2.,_Sz);
         k4 = phi(_rho[i]+dx*k3,_Sz);
         _rho[i+1] = _rho[i] + (dx/6.)*(k1+2.*k2+2.*k3+k4);
      }

      return;
   }

}; // End class

} // ns hydrodynamics
} // ns mfem

/***
 * Exact Solution process outlined in supplementary.pdf
 * 
 * 
 * 
 * 
 * 
 * 
 * 
*/