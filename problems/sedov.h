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
   bool distort_mesh = false;
   bool known_exact_solution = true;
   bool bcs = false; // Indicator for boundary conditions
   string indicator = ""; // Possible: saltzmann

   // CFL change
   bool _change_cfl = false;

   // Constants specific to Sedov problem
   double h = 0., cell_vol = 0.;

public:
   SedovProblem()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
   }

   /* Override getters */
   bool get_distort_mesh() override { return distort_mesh; }
   bool has_exact_solution() override { return known_exact_solution; }

   bool change_cfl() override { return _change_cfl; }
   // double get_cfl_first() override { return CFL_first; }
   // double get_cfl_second() override { return CFL_second; }
   // double get_cfl_time_change() override { return CFL_time_change; }

   /* Override specific update functions */
   void lm_update(const double b_covolume) override 
   {
      this->set_b(b_covolume);
   }
   virtual void update(Vector x_gf, double t) override {
      if (t <= 1.e-16)
      {  
         cout << "setting internal energy params.\n";
         this->h = x_gf[0];
         this->cell_vol = x_gf[1];
      }
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double pressure(const Vector &U) override
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
            return 0.979264 / (4. * cell_vol);
            // return 2. / (4. * cell_vol);
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
   void set_internal_energy_params(double _h, double _cell_vol)
   {
      cout << "setting internal energy params.\n";
      this->h = h;
      this->cell_vol = _cell_vol;
   }
}; // End class

} // ns hydrodynamics
} // ns mfem