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
class SodProblem: public ProblemBase<dim>
{
public:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double a = 0.;     
   double b = 0.;
   double gamma = 1.4;
   bool distort_mesh = false;
   bool known_exact_solution = false;

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double pressure(const Vector &U) override
   {
      // Initial
      // _pL = 1.0;
      // _pR= 0.1;
      return (gamma - 1.) * this->internal_energy(U);
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double rho0(const Vector &x, const double & t) override
   {
      switch (dim)
      {
         case 1:
         case 2:
         {
            return (x(0) < 0.5) ? 1.0 : 0.125;
         }
         default:
         {
            MFEM_ABORT("Invalid dimension provided.\n");
         }
      }
      return 1.;
   }
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      v = 0.;
      return;
   }

   double sie0(const Vector &x, const double & t) override
   {
      return (x(0) < 0.5) ? 1.0 / this->rho0(x, t) / (gamma - 1.0) // Sod
                        : 0.1 / this->rho0(x, t) / (gamma - 1.0);
   }

}; // End class

} // ns hydrodynamics
} // ns mfem