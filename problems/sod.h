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
#include "riemann1D.hpp"
#include <cmath>

using namespace mfem;
using namespace std;

namespace mfem
{
namespace hydrodynamics
{

template<int dim>
class SodProblem: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 1.4;
   bool distort_mesh = false;
   bool known_exact_solution = true;
   bool bcs = true; // Indicator for boundary conditions
   string _indicator = "Sod"; // Possible: saltzmann

   double rhoL = 1.0, rhoR = 0.125, pL = 1.0, pR = 0.1, vL = 0., vR = 0.;
   double x_center = 0.5;

public:
   SodProblem()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_bcs_indicator(bcs);
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
      // Initial
      // _pL = 1.0;
      // _pR= 0.1;
      return (this->get_gamma() - 1.) * this->internal_energy(U);
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double pressure(const Vector &x, const double & t)
   {
      if (t < 1e-12)
      {
         return (x(0) < x_center) ? pL : pR;
      }
      else
      {
         double params[8];
         params[0] = rhoL; params[3] = rhoR; // rho
         params[1] = pL; params[4] = pR;     // p
         params[2] = params[5] = vL;         // u
         params[6] = this->get_gamma();      // gamma
         params[7] = x_center;               // x_center
         riemann1D::init(params);

         double _p[2];
         _p[0] = x[0]; _p[1] = t;
         return riemann1D::p(_p);
      }
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double rho0(const Vector &x, const double & t) override
   {
      if (t < 1e-12)
      {
         switch (dim)
         {
            case 1:
            case 2:
            {
               return (x(0) < x_center) ? rhoL : rhoR;
            }
            default:
            {
               MFEM_ABORT("Invalid dimension provided.\n");
            }
         }
      }
      else // use exact solution in riemann1D.hpp
      {
         double params[8];
         params[0] = rhoL; params[3] = rhoR; // rho
         params[1] = pL; params[4] = pR;     // p
         params[2] = params[5] = vL;         // u
         params[6] = this->get_gamma();      // gamma
         params[7] = x_center;               // x_center
         riemann1D::init(params);

         double p[2];
         p[0] = x[0]; p[1] = t;
         return riemann1D::rho(p);
      }
   }
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      if (t < 1e-12)
      {
         v = vL;
      }
      else
      {
         double params[8];
         params[0] = rhoL; params[3] = rhoR; // rho
         params[1] = pL; params[4] = pR;     // p
         params[2] = params[5] = vL;         // u
         params[6] = this->get_gamma();      // gamma
         params[7] = x_center;               // x_center
         riemann1D::init(params);

         double p[2];
         p[0] = x[0]; p[1] = t;
         v[0] = riemann1D::v(p);
      }
      return;
   }

   double sie0(const Vector &x, const double & t) override
   {
      return (x(0) < x_center) ? pressure(x,t) / this->rho0(x, t) / (this->get_gamma() - 1.0) // Sod
                        : pressure(x,t) / this->rho0(x, t) / (this->get_gamma() - 1.0);
   }

}; // End class

} // ns hydrodynamics
} // ns mfem