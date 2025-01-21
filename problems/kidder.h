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
   const double _a = 0., _b = 0., _gamma = 2.;
   const double r1 = 0.9, r2 = 1.0, rho1 = 1., rho2 = 2., s = 1.;
   const double P1 = pow(rho1, _gamma), P2 = pow(rho2, _gamma);
   const double tau;
   bool _distort_mesh = false;
   bool _known_exact_solution = true;
   bool _bcs = false; // Indicator for boundary conditions
   string _indicator = "Kidder"; // Possible: saltzmann

   /* Initialize focusing time */
   double compute_focusing_time()
   {

      double c1 = sqrt(this->get_gamma()*P1/rho1);
      double c2 = sqrt(this->get_gamma()*P2/rho2);

      double num = (this->get_gamma() - 1.) * (pow(r2,2) - pow(r1,2));
      double denom = 2. * (pow(c2,2) - pow(c1,2));
      double tau = sqrt(num / denom);
      return tau;
   }

public:
   KidderProblem() : tau(compute_focusing_time())
   {
      cout << "rho1: " << rho1 << endl;
      cout << "rho2: " << rho2 << endl;
      cout << "r1: " << r1 << endl;
      cout << "r2: " << r2 << endl;
      cout << "P1: " << P1 << endl;
      cout << "P2: " << P2 << endl;
      cout << "focusing time: " << tau << endl;
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_bcs_indicator(_bcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);
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
    * NOTE: the exact solution for pressure, density, and 
    * velocity are in terms of initial particle locations.
    * Hence these functions assume that the first parameter 
    * is the initial particle location, rather than the 
    * particle location at time t.
    *********************************************************/
   double p0(const Vector &x, const double & t) override
   {
      /* Get initial pressure */
      double _p0 = pow(rho0(x,0), this->get_gamma());

      if (t < 1e-12)
      {
         return _p0;
      }
      else 
     {
         double _h = h(t);
         return _p0 * pow(_h, 2. * this->get_gamma() / (1. - this->get_gamma()));
      }
      return 0.;
   }
   double rho0(const Vector &x, const double & t) override
   {
      double r = x.Norml2();
      /* Get initial density */
      double _rho0 = (pow(r2,2) - pow(r,2)) / (pow(r2,2) - pow(r1,2)) * pow(rho1, this->get_gamma()-1);
      _rho0 += (pow(r,2) - pow(r1,2)) / (pow(r2,2) - pow(r1,2)) * pow(rho2, this->get_gamma()-1);
      _rho0 = pow(_rho0, 1./(this->get_gamma()-1.));
      if (t < 1e-12)
      { 
         return _rho0;
      }
      else 
      {
         double _h = h(t);
         return _rho0 * pow(_h, 2. / (1. - this->get_gamma()));
      }
   }
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      if (t < 1e-12)
      {
         v = 0.;
      }
      else 
      {
         // double r = x.Norml2();
         double _dhdt = dhdt(t);
         v = x;
         v *= -_dhdt; // ||v|| is r * dhdt, however must divide x by r to normalize. So r cancels
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
   double h(const double t)
   {
      return sqrt(1. - pow(t/tau, 2));
   }
   double dhdt(const double t)
   {
      return t / pow(tau, 2) / h(t);     
   }

   void GetBoundaryState(const Vector &x, const double &t, const int &bdr_attr, Vector &state) override
   { 
      // cout << "Kidder::GetBoundaryState\n";
      state.SetSize(dim + 2);
      state = 0.;

      Vector vel(dim);
      double norm_x = x.Norml2();
      vel = x;
      vel *= -1.;
      vel /= norm_x;
      double _dhdt = dhdt(t);
      double hmult = pow(h(t), 2. * this->get_gamma() / (1. - this->get_gamma()));
      double P1t = P1 * hmult;
      double P2t = P2 * hmult;
      double rho1t = rho1 * hmult;
      double rho2t = rho2 * hmult;

      switch (bdr_attr)
      {
      case 4: // inner radius
      {
         // Specific volume
         state[0] = 1. / rho1t;
         // Velocity
         vel *= r1 * _dhdt;
         state[1] = vel[0];
         state[2] = vel[1];
         // Specific total energy
         const double ke = vel * vel;
         state[3] = P1t / rho1t / (this->get_gamma() - 1.) + 0.5 * ke;
         break;
      }
      case 5: // outer radius
      {
         // Specific volume
         state[0] = 1. / rho2t;
         // Velocity
         vel *= r2 * _dhdt;
         state[1] = vel[0];
         state[2] = vel[1];
         // Specific total energy
         const double ke = vel * vel;
         state[3] = P2t / rho2t / (this->get_gamma() - 1.) + 0.5 * ke;
         break;
      }
      default:
      {
         MFEM_ABORT("Invalid boundary condition for Kidder problem.\n");
      }
      } 
   }

   virtual void get_additional_BCs(const FiniteElementSpace &fes, Array<int> ess_bdr, Array<int> &add_ess_tdofs, Array<double> &add_bdr_vals) override { }
   virtual void update_additional_BCs(const double &t, const double timestep_first, Array<double> &add_bdr_vals) override { }
}; // End class

} // ns hydrodynamics
} // ns mfem