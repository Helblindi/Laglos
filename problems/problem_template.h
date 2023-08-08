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
class ProblemTemplate: public ProblemBase
{
public:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   a = 0.;     
   b = 0.;
   gamma = 0.;
   distort_mesh = false;
   known_exact_solution = false;

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double pressure(const Vector &U)
   {
      /*
      Must Override
      */
      return 0;
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double rho0(const Vector &x, const double & t)
   {
      /*
      Must Override
      */
      return 0;
   }
   void v0(const Vector &x, const double & t, Vector &v)
   {
      /*
      Must Override
      */
      return 0;
   }
   void sie0(const Vector &x, const double & t)
   {
      /*
      Must Override
      */
      return 0;
   }

}; // End class

} // ns hydrodynamics
} // ns mfem