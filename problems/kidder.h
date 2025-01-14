/*********************************************************
* This file serves as a template for setting up a problem.
* Steps to implement a problem:
*  1) Copy this file and rename in $PROJECT_SRC/include/ to 
*     correspond to your problem
*  2) Modify problem specific constants a, b, this->get_gamma()
*  3) Override pressure function
*  4) Specify initial conditions given by rho0, v0, and sie0
*  5) Add corresponding problem option to laglos.cpp
*  6) Add to include header file problems/test_problems_include.h
*  7) Enforce boundary conditions in laglos_solver.cpp
*     a) This may require edits to the following functions 
*        depending on the type of BCs that will be 
*        implemented
*        - CreateBdrElementIndexingArray()
*        - CreateBdrVertexIndexingArray()
*        - FillCellBdrFlag()
*        - ComputeStateUpdate()
*        - EnforceMVBoundaryConditions()
*  8) Edit SaveStateVecsToFile if the problem is a radial one.
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
class KidderProblem: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   const double r1 = 0.9, r2 = 1.0, P1 = 0.1, P2 = 10., rho2 = .01;
   const double _a = 0., _b = 0., _gamma = 5./3.;
   bool _distort_mesh = false;
   bool _known_exact_solution = false;
   bool _bcs = false; // Indicator for boundary conditions
   string _indicator = "Kidder"; // Possible: saltzmann

   // CFL change, can remove if not needed
   bool _change_cfl = false;
   double _cfl_first = 0.5;
   double _cfl_second = 0.5;
   double _cfl_time_change = 0.01;

public:
   KidderProblem()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_bcs_indicator(_bcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);
      // CFL change
      this->set_cfl_change(_change_cfl);
      this->set_cfl_first(_cfl_first);
      this->set_cfl_second(_cfl_second);
      this->set_cfl_time_change(_cfl_time_change);
   }

   /* Optionally overridden, or removed */
   double get_gamma(const int &cell_attr = 0) override { return _gamma; }
   void lm_update(const double b_covolume) override {}
   void update(Vector vec, double t = 0.) override {}

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
   double p0(const Vector &x, const double & t) override
   {
      double r = x.Norml2();
      if (t < 1e-12)
      {
         double s = P2 / (pow(rho2, this->get_gamma()));
         if (r < r1)
         {
            return P1;
         }
         else if (r < r2)
         {
            return s * pow(rho0(x,t), this->get_gamma());
         }
         else 
         {
            return P2;
         }
      }
      else 
      {
         MFEM_ABORT("P0 must be overridden\n");
      }
      return 0.;
   }
   double rho0(const Vector &x, const double & t) override
   {
      double r = x.Norml2();
      if (t < 1e-12)
      { 
         double rho1 = rho2 * pow(P1/P2, 1./this->get_gamma());
         if (r < r1)
         {
            return rho1;
         }
         else if (r < r2)
         {
            double val = (pow(r2,2) - pow(r,2)) / (pow(r2,2) - pow(r1,2)) * pow(rho1, this->get_gamma()-1);
            val += (pow(r,2) - pow(r1,2)) / (pow(r2,2) - pow(r1,2)) * pow(rho2, this->get_gamma()-1);
            val = pow(val, 1./(this->get_gamma()-1.));
            return val;
         }
         else 
         {
            return rho2;
         }
         
      }
      else 
      {
         MFEM_ABORT("rho0 must be overridden\n");
      }

      return 0.;
   }
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      if (t < 1e-12)
      {
         v = 0.;
      }
      else 
      {
         MFEM_ABORT("v0 must be overridden\n");
      }
      return;
   }
   double sie0(const Vector &x, const double & t) override
   {
      return p0(x,t) / this->rho0(x, t) / (this->get_gamma() - 1.0);
   }

   /*********************************************************
    * Problem specific functions
    *********************************************************/
   double compute_focusing_time()
   {
      double rho1 = rho2 * pow(P1/P2, 1./this->get_gamma());
      double s = P2 / (pow(rho2, this->get_gamma()));
      double c1 = sqrt(s*this->get_gamma()*pow(rho1, this->get_gamma()-1.));
      double c2 = sqrt(s*this->get_gamma()*pow(rho2, this->get_gamma()-1.));

      double num = (this->get_gamma() - 1.) * (pow(r2,2) - pow(r1,2));
      double denom = 2. * (pow(c2,2) - pow(c1,2));
      double tau = sqrt(num / denom);
      return tau;
   }
   double h(const double t)
   {
      double tau = compute_focusing_time();
      return sqrt(1. - pow(t/tau, 2));
   }

   void GetBoundaryState(const double &t, const int &bdr_attr, Vector &state) override
   { 
      state.SetSize(dim + 2);
      state = 0.;

      double hmult = pow(h(t), 2. * this->get_gamma() / (1. - this->get_gamma()));
      double P1t = P1 * hmult;
      double P2t = P2 * hmult;
      double s = P2 / (pow(rho2, this->get_gamma()));
      double rho2t = pow(P2t / s, 1. / this->get_gamma());
      double rho1t = rho2t * pow(P1t/P2t, 1./this->get_gamma());

      switch (bdr_attr)
      {
      case 4: // inner radius
      {
         double rho1 = rho2 * pow(P1/P2, 1./this->get_gamma());
         state[0] = rho1;
         state[3] = P1t / rho1 / (this->get_gamma() - 1.);
      }
      case 5: // outer radius
      {
         state[0] = rho2;
         state[3] = P2t / rho2 / (this->get_gamma() - 1.);
      }
      } 
   }

   virtual void get_additional_BCs(const FiniteElementSpace &fes, Array<int> ess_bdr, Array<int> &add_ess_tdofs, Array<double> &add_bdr_vals) override { }
   virtual void update_additional_BCs(const double &t, const double timestep_first, Array<double> &add_bdr_vals) override { }
}; // End class

} // ns hydrodynamics
} // ns mfem