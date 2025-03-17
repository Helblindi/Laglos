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
class SodRadial: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 1.4;
   bool _distort_mesh = false;
   bool _known_exact_solution = true;
   bool _bcs = true; // Indicator for boundary conditions
   string _indicator = "SodRadial"; // Possible: saltzmann

   // Initial conditions
   // Here, L is inside of shock bubble, R is outside
   double rhoL = 1.0, rhoR = 0.125, pL = 1.0, pR = 0.1, vL = 0., vR = 0.;
   double x_center = 0.5;


public:
   SodRadial()
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

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double pressure(const double &rho, const double &sie, const int &cell_attr=0) override
   {
      return (this->get_gamma() - 1.) * rho * sie;
   }

   double p0(const Vector &x, const double & t) override
   {
      double norm_x = x.Norml2();
      if (t < 1e-12)
      {
         return (norm_x < x_center) ? pL : pR;
      }
      else
      {
         double params[8];
         params[0] = rhoL; params[3] = rhoR; // rho
         params[1] = pL; params[4] = pR;     // p
         params[2] = params[5] = vL;         // u
         params[6] = this->get_gamma();      // gamma
         params[7] = x_center;               // x_center
         riemann1D::init(params);

         double _p[2];
         _p[0] = norm_x; _p[1] = t;
         return riemann1D::p(_p);
      }
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double rho0(const Vector &x, const double & t) override
   {
      double norm_x = x.Norml2();
      if (t < 1e-12)
      {
         switch (dim)
         {
            case 1:
            case 2:
            {
               return (norm_x < x_center) ? rhoL : rhoR;
            }
            default:
            {
               MFEM_ABORT("Invalid dimension provided.\n");
            }
         }
      }
      else // use exact solution in riemann1D.hpp
      {
         double params[8];
         params[0] = rhoL; params[3] = rhoR; // rho
         params[1] = pL; params[4] = pR;     // p
         params[2] = params[5] = vL;         // u
         params[6] = this->get_gamma();      // gamma
         params[7] = x_center;               // x_center
         riemann1D::init(params);

         double p[2];
         p[0] = norm_x; p[1] = t;
         return riemann1D::rho(p);
      }
   }
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      if (t < 1e-12)
      {
         v = vL;
      }
      else
      {
         double norm_x = x.Norml2();

         double params[8];
         params[0] = rhoL; params[3] = rhoR; // rho
         params[1] = pL; params[4] = pR;     // p
         params[2] = params[5] = vL;         // u
         params[6] = this->get_gamma();      // gamma
         params[7] = x_center;               // x_center
         riemann1D::init(params);

         double p[2];
         p[0] = norm_x; p[1] = t;
         double v_coeff = riemann1D::v(p);
         v = x;
         v *= v_coeff / norm_x;
      }
      return;
   }

   double sie0(const Vector &x, const double & t) override
   {
      double norm_x = x.Norml2();

      return p0(x,t) / this->rho0(x, t) / (this->get_gamma() - 1.0);
   }

}; // End class

} // ns hydroLO
} // ns mfem