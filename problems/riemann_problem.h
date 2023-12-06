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
class RiemannProblem: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 1.4;
   bool _distort_mesh = false;
   bool _known_exact_solution = false;
   bool _bcs = false; // Indicator for boundary conditions
   string _indicator = ""; // Possible: saltzmann

   double rhoL = 1., rhoR = 1.;
   double vL = 0., vR = 0.;
   double pL = .1, pR = .1;
   double x0 = 0.;

public:
   RiemannProblem()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_bcs_indicator(_bcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);
   }

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
         if (x[0] <= x0)
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
         if (x[0] <= x0)
         {
            v[0] = vL;
            return;
         }
         else
         {
            v[0] = vR;
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
         double rho = rho0(x,t);
         double pressure = p0(x);

         return ((pressure + this->get_a() * pow(rho,2)) * (1. - this->get_b()*rho)  / (rho * (this->get_gamma() - 1.))) - this->get_a() * rho;
      }
      else {
         MFEM_ABORT("Exact solution for vdw3 not programmed.\n");
         return -1.;
      }
   }

   // Initial values are in terms of pressure
   double p0(const Vector &x)
   {
      if (x[0] <= x0)
      {
         return pL;
      }
      else 
      {
         assert(x[0] <= 1.);
         return pR;
      }
   }

}; // End class

} // ns hydrodynamics
} // ns mfem