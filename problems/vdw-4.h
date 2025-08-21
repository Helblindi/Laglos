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

using namespace mfem;
using namespace std;

namespace mfem
{
namespace hydroLO
{

class VdwTest4: public ProblemBase
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 1., _b = 1., _gamma = 1.02;
   bool _distort_mesh = false;
   bool _known_exact_solution = false;
   bool _thbcs = false;
   bool _mvbcs = false;
   string _indicator = "Vdw4";

   // Problem specifics
   double initial_shock = 0.0;
   double rhoL = 0.9932, rhoR = 0.95;
   double vL = 3., vR = -3.;
   double pL = 2.;
   double pR = 2.;
   // double eL = 0.029143658477667977;
   // double eR = 6.688157894736825

public:
   VdwTest4(const int &_dim) : ProblemBase(_dim)
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_thbcs_indicator(_thbcs);
      this->set_mvbcs_indicator(_mvbcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);

      // Set Equation of state
      this->eos = std::unique_ptr<EquationOfState>(new VanDerWaalsEOS(this->a, this->b));
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   virtual double rho0(const Vector &x, const double & t) const override
   {
      if (t < 1.e-16) {
         if (x[0] <= initial_shock)
         {
            return rhoL;
         }
         else
         {
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
         if (x[0] <= initial_shock)
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
         double rho = rho0(x,t);
         double pressure = p0(x,t);
         return this->eos->energy(pressure, rho, this->get_gamma());
      }
      else {
         MFEM_ABORT("Exact solution for vdw3 not programmed.\n");
         return -1.;
      }
   }

   // Initial values are in terms of pressure
   virtual double p0(const Vector &x, const double & t) const override
   {
      if (t < 1.E-12)
      {
         if (x[0] <= initial_shock)
         {
            return pL;
         }
         else 
         {
            return pR;
         }
      }
      else {
         MFEM_ABORT("Exact solution for vdw4 not programmed.\n");
         return -1.;
      }
   }

}; // End class

} // ns hydroLO
} // ns mfem