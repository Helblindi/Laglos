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
*     a) This may require edits to the following functions 
*        depending on the type of BCs that will be 
*        implemented
*        - LagrangianLOOperator::CreateBdrElementIndexingArray()
*        - LagrangianLOOperator::CreateBdrVertexIndexingArray()
*        - LagrangianLOOperator::FillCellBdrFlag()
*        - LagrangianLOOperator::UpdateMeshVelocityBCs()
*        - LagrangianLOOperator::EnforceL2BC()
*        - ProblemBase::get_mv_bcs_need_updating()
*        - ProblemBase::update_additional_BCs()
*     b) To note about boundary conditions in the mesh file:
*        1: vx = 0
*        2: vy = 0
*        3: vz = 0
*        4: vr = 0 (radial boundary condition)
*        5: Misc (possibly dirichlet condition)
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
namespace hydroLO
{

template<int dim>
class RiemannProblem: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0;
   bool _distort_mesh = false;
   bool _known_exact_solution = false;
   bool _thbcs = false; // Indicator for thermal boundary conditions
   bool _mvbcs = true; // Indicator for mv boundary conditions
   string _indicator = "riemann"; // Possible: saltzmann

   /* Top states of TP */
   // double _gamma_1 = 1.5, _gamma_2 = 1.5;
   // double rhoL = 1., rhoR = .125;
   // double vL = 0., vR = 0.;
   // double eL = 2., eR = 1.6;

   /* Bottom states of TP */
   // double _gamma_1 = 1.5, _gamma_2 = 1.4;
   // double rhoL = 1., rhoR = 1.;
   // double vL = 0., vR = 0.;
   // double eL = 2., eR = 0.25;

   /* bottom = Left, top = Right*/
   double _gamma_1 = 1.4, _gamma_2 = 1.5;
   double rhoL = 1., rhoR = 0.125;
   double vL = 0., vR = 0.;
   double eL = 0.25, eR = 1.6;

   double x0 = 1.;

public:
   RiemannProblem()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_indicator(_indicator);
      this->set_thbcs_indicator(_thbcs);
      this->set_mvbcs_indicator(_mvbcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);
   }

   double get_gamma(const int &cell_attr) override {
      assert(cell_attr != 0 && "Must pass in a cell_attr to any ProblemBase::get_gamma funcalls.\n");
      return (cell_attr == 1) ? _gamma_1 : _gamma_2;
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   virtual double pressure(const double &rho, const double &sie, const int &cell_attr=0) override
   {
      // Use van der Waals
      // double rho = 1. / U[0];
      // double sie = this->specific_internal_energy(U);

      // double val = (this->get_gamma(cell_attr) - 1.) * (rho * sie + this->get_a() * pow(rho, 2)) / (1. - this->get_b() * rho) - this->get_a() * pow(rho,2);
      // return val;
      
      /* How the TriplePoint problem handles pressure */
      double _g;
      if (cell_attr == 1) { 
         _g = _gamma_1;
      } else {
         assert(cell_attr != 0 && "Must pass in a cell_attr to any ProblemBase::pressure funcalls.\n");
         _g = _gamma_2;
      }
      return (_g - 1.) * rho * sie;
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double p0(const Vector &x, const double &t) override
   {
      double _gamma = 0.;
      if (x[0] <= 1)
      {
         _gamma = _gamma_1;
      }
      else
      {
         _gamma = _gamma_2;
      }
      return (_gamma - 1.) * rho0(x,t) * sie0(x,t);
   }

   virtual double rho0(const Vector &x, const double & t) override
   {
      if (t < 1.e-16) {
         if (x[0] <= x0)
         {
            return rhoL;
         }
         else
         {
            // assert(x[0] <= 1.);
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

         if (x[0] <= 1)
         {
            return eL;
         }
         else
         {
            return eR;
         }
      }
      MFEM_ABORT("Should not be calling sie0 unless initializing.\n");
   }

}; // End class

} // ns hydroLO
} // ns mfem