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
class SaltzmannProblem: public ProblemBase<dim>
{
public:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double a = 0., b = 0., gamma = 5./3.;
   bool distort_mesh = true;
   bool known_exact_solution = true;
   bool bcs = true;
   string indicator = "saltzmann";

   // CFL change
   bool _change_cfl = true;
   constexpr static double CFL_first = 0.01;
   constexpr static double CFL_second = 0.25;
   constexpr static double CFL_time_change = 0.01; // From Boscheri's paper

   virtual bool change_cfl() override { return _change_cfl; }
   virtual double get_cfl_first() override { return CFL_first; }
   virtual double get_cfl_second() override { return CFL_second; }
   virtual double get_cfl_time_change() override { return CFL_time_change; }

   double rotation_angle = 0.; // 0 - 1D horizontal velocity
   double rhoL = 1.0, rhoR = 1.0;;
   double rhoLstar = 3.9992502342988532;
   double rhoRstar = 3.9992502342988532;
   double vL = 2., vR = 0., vstar = 1;
   double pL = (this->get_gamma() - 1.) * pow(10., -4); // this is actually sie
   double pR= (this->get_gamma() - 1.) * pow(10., -4);  // also sie
   double pstar = 1.3334833281256511;
   double lambda1m = 0.66658333854101581;
   double lambda1p = 0.66658333854101581; 
   double lambda3 = 1.3334166614589844;

   double cL = sqrt(this->get_gamma() * pL / rhoL);
   double x0 = 0.0000000001; // Initial shock position

   /* Override getters */
   virtual double get_a() override { return a; }
   virtual double get_b() override { return b; }
   virtual double get_gamma() override { return gamma; }
   virtual string get_indicator() override { return indicator; }
   virtual bool get_distort_mesh() override { return distort_mesh; }
   virtual bool has_exact_solution() override { return known_exact_solution; }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   virtual double pressure(const Vector &U) override
   {
      return (this->get_gamma() - 1.) * this->internal_energy(U);
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   virtual double pressure(const Vector &x, const double &t) 
   {
      if (t == 0) { return pow((this->get_gamma() - 1), 2) * 10e-4; } // TODO: Change pressure
      else {
         Vector _n(dim);
         _n[0] = cos(rotation_angle);
         if (dim > 1)
         {
            _n[1] = sin(rotation_angle);
         }
         double _x_tilde = x * _n;

         double _xi = (_x_tilde - x0) / t;

         // Saltzman TODO
         if (_x_tilde <= lambda1p * t)
         {
            return pL;
         }
         else if (_x_tilde <= lambda3 * t)
         {
            return pstar;
         }
         else
         {
            return pR;
         }
      }
   }

   virtual double rho0(const Vector &x, const double & t) override
   {
      assert(dim < 3); //
      if (t == 0) { return 1.; }
      else 
      {
         Vector _n(dim);
         _n[0] = cos(rotation_angle);
         if (dim > 1)
         {
            _n[1] = sin(rotation_angle);
         }
         double _x_tilde = x * _n;

         double _xi = (_x_tilde - x0) / t;

         // Saltzman TODO
         if (_x_tilde <= lambda1p * t)
         {
            return rhoL;
         }
         else if (_x_tilde <= lambda3 * t)
         {
            return rhoRstar;
         }
         else
         {
            return rhoR;
         }
      }
      return 0.;
   }
   virtual void v0(const Vector &x, const double & t, Vector &v) override
   {
      if (t == 0) { v = 0.; }
      else 
      {
         Vector _n(dim);
         _n[0] = cos(rotation_angle);
         if (dim > 1)
         {
            _n[1] = sin(rotation_angle);
         }
         double _x_tilde = x * _n;

         double _xi = (_x_tilde - x0) / t;

         // our velocity will be some scaling of the normal vector
         v = _n;

         // Saltzman TODO
         if (_x_tilde <= lambda1p * t)
         {
            v *= vL;
         }
         else if (_x_tilde <= lambda3 * t)
         {
            v *= vstar;
         }
         else
         {
            v *= vR;
         }
      }
      return;
   }
   virtual double sie0(const Vector &x, const double & t) override
   {
      if (t == 0) { return pow(10, -4); }
      else {
         double _p = pressure(x,t);
         double _rho = rho0(x,t);
         return _p / (_rho * (this->get_gamma() - 1.));
      }
   }

}; // End class

} // ns hydrodynamics
} // ns mfem