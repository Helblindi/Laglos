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
class TriplePoint: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma_1 = 1.5, _gamma_2 = 1.4;
   bool _distort_mesh = false;
   bool _known_exact_solution = false;
   bool _bcs = true; // Indicator for boundary conditions
   string _indicator = "TriplePoint"; // Possible: saltzmann

   // CFL change
   bool _change_cfl = false;

public:
   TriplePoint()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_indicator(_indicator);
      this->set_bcs_indicator(_bcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);
   }

   /* Override getters */
   double get_gamma(const int &cell_attr) override {
      assert(cell_attr != 0 && "Must pass in a cell_attr to any ProblemBase::get_gamma funcalls.\n");
      return (cell_attr == 1) ? _gamma_1 : _gamma_2;
   }

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
      // TODO: Must fix to use cell specific gamma
      double _g;
      if (cell_attr == 1) { 
         _g = _gamma_1;
      } else {
         assert(cell_attr != 0 && "Must pass in a cell_attr to any ProblemBase::pressure funcalls.\n");
         _g = _gamma_2;
      }
      return (_g - 1.) * this->internal_energy(U);
      
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double rho0(const Vector &x, const double & t) override
   {
      return (dim == 2) ? (x(0) > 1.0 && x(1) > 1.5) ? 0.125 : 1.0
         /* dim = 3 */  : x(0) > 1.0 && ((x(1) < 1.5 && x(2) < 1.5) ||
                                         (x(1) > 1.5 && x(2) > 1.5)) ? 0.125 : 1.0;
      // return (dim == 2) ? (x(0) > 1.0 && x(1) < 1.5) ? 1. : .125; // same pressure
   }
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      v = 0.;
      return;
   }
   double sie0(const Vector &x, const double & t) override
   {
      if (x[0] <= 1)
      {
         // D_1
         return 2.;
         // return 1.6; // same pressure
      }
      else
      {
         if (x[1] <= 1.5)
         {
            // D_2
            return 0.25;
         }
         else 
         {
            // D_3
            assert(x[1] <= 3.);
            return 1.6;
         }
      }
   }

}; // End class

} // ns hydrodynamics
} // ns mfem