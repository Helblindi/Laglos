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
#include "sedov_exact.hpp"
#include <cmath>
#include <string>

using namespace mfem;
using namespace std;

namespace mfem
{
namespace hydrodynamics
{

template<int dim>
class SedovProblem: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 1.4;
   bool _distort_mesh = false;
   bool _known_exact_solution = true;
   bool _bcs = true; // Indicator for boundary conditions
   string _indicator = "Sedov"; // Possible: saltzmann

   // Constants specific to Sedov problem
   double h = 1., cell_vol = 1.;
   // double sedov_energy_initial = 0.979264;
   double sedov_energy_initial = 0.3;
   // double sedov_energy_initial = 0.025;

public:
   SedovProblem()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_bcs_indicator(_bcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);
   }

   /* Override specific update functions */
   void lm_update(const double b_covolume) override 
   {
      this->set_b(b_covolume);
   }
   virtual void update(Vector params, double t) override {
      // params is a vector [hmax, cell_vol]
      if (t <= 1.e-16)
      {
         this->h = params[0];
         this->cell_vol = params[1];
      }
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double pressure(const Vector &U, const int &cell_attr=0) override
   {
      return (this->get_gamma() - 1.) * this->internal_energy(U);
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double rho0(const Vector &x, const double & t) override
   {
      if (t < 1.e-16)
      {
         // Initial condition 
         return 1.;
      }
      else 
      {
         if (dim != 2)
         {
            MFEM_ABORT("Dimension != 2 is not supported.\n");
         }
         double xyt[3];
         xyt[0] = x[0], xyt[1] = x[1], xyt[2] = t;
         return sedov_rho(xyt);
      }
   }
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      v = 0.;

      if (t < 1.e-16)
      {
         // Initial condition
         return;
      }
      else
      {
         double xyt[3];
         xyt[0] = x[0], xyt[1] = x[1], xyt[2] = t;

         v[0] = sedov_vx(xyt);
         v[1] = sedov_vy(xyt);
      }

      return;
   }
   double sie0(const Vector &x, const double & t) override
   {
      assert(h != 0.);
      assert(cell_vol != 0.);
      double norm = x.Norml2();
      if (t < 1.e-16)
      {
         // Initial condition
         if (norm <= h) {
            assert(cell_vol > 0.);
            return sedov_energy_initial / cell_vol;
         }
         else
         {
            return 0.;
         }
      }
      else
      {
         double xyt[3];
         xyt[0] = x[0], xyt[1] = x[1], xyt[2] = t;

         return sedov_e(xyt);
      }
   }

   /*********************************************************
    * Sedov specific functions
    *********************************************************/
}; // End class

} // ns hydrodynamics
} // ns mfem