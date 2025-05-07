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
   bool _bcs = true; // Indicator for boundary conditions
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
      this->set_bcs_indicator(_bcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);

      // Set Equation of state
      this->eos = std::unique_ptr<EquationOfState>(new IdealGasEOS());
   }

   double get_gamma(const int &cell_attr) const override {
      assert(cell_attr != 0 && "Must pass in a cell_attr to any ProblemBase::get_gamma funcalls.\n");
      return (cell_attr == 1) ? _gamma_1 : _gamma_2;
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double p0(const Vector &x, const double &t) const override
   {
      double rho = rho0(x,t);
      double sie = sie0(x,t);
      double _gamma = 0.;
      if (x[0] <= 1)
      {
         _gamma = _gamma_1;
      }
      else
      {
         _gamma = _gamma_2;
      }
      return this->eos->pressure(rho, sie, _gamma);
   }

   virtual double rho0(const Vector &x, const double & t) const override
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
   virtual void v0(const Vector &x, const double & t, Vector &v) const override
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
   virtual double sie0(const Vector &x, const double & t) const override 
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