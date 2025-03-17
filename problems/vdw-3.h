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

template<int dim>
class VdwTest3: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 1., _b = 1., _gamma = 1.02;
   bool _distort_mesh = false;
   bool _known_exact_solution = false;
   bool _bcs = false;
   string _indicator = "Vdw3";

   // Problem specifics
   double initial_shock = 0.0;
   double rhoL = 0.25, rhoR = 0.000049;
   double vL = 0., vR = 0.;
   double pL = 3.e-2;
   double pR = 5.e-8;

public:
   VdwTest3()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_bcs_indicator(_bcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   virtual double pressure(const double &rho, const double &sie, const int &cell_attr=0) override
   {
      // Use van der Waals
      double val = (this->get_gamma() - 1.) * (rho * sie + this->get_a() * pow(rho, 2)) / (1. - this->get_b() * rho) - this->get_a() * pow(rho,2);
      return val;
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   virtual double rho0(const Vector &x, const double & t) override
   {
      if (t < 1.e-16) {
         if (x[0] <= initial_shock)
         {
            return rhoL;
         }
         else
         {
            assert(x[0] <= 1.);
            return rhoR;
         }
      }
      else {
         MFEM_ABORT("Exact solution for vdw3 not programmed.\n");
         return -1.;
      }
   }
   virtual void v0(const Vector &x, const double & t, Vector &v) override
   {
      v = 0.;
      if (t < 1.e-16) {
         if (x[0] <= initial_shock)
         {
            v[0] = vL;
         }
         else
         {
            assert(x[0] <= 1.);
            v[0] = vR;
         }
      }
   }

   virtual double sie0(const Vector &x, const double & t) override 
   {
      if (t < 1.e-16) {
         double rho = rho0(x,t);
         double pressure = p0(x);

         return ((pressure + this->get_a() * pow(rho,2)) * (1. - this->get_b()*rho)  / (rho * (this->get_gamma() - 1.))) - this->get_a() * rho;
      }
      else {
         MFEM_ABORT("Exact solution for vdw3 not programmed.\n");
         return -1.;
      }
   }

   // Initial values are in terms of pressure
   double p0(const Vector &x)
   {
      if (x[0] <= initial_shock)
      {
         return pL;
      }
      else 
      {
         assert(x[0] <= 1.);
         return pR;
      }
   }

}; // End class

} // ns hydroLO
} // ns mfem