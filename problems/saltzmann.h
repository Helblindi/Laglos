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
class SaltzmannProblem: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   mutable double timestep_first = 0.;
   double _a = 0., _b = 0., _gamma = 5./3.;
   bool _distort_mesh = true;
   bool _known_exact_solution = true;
   bool _bcs = true;
   bool _mv_bcs_need_updating = true;
   string _indicator = "saltzmann";

   // CFL change
   bool _change_cfl = true;
   double _cfl_first = 0.01;
   double _cfl_second = 0.5;
   double _cfl_time_change = 0.01; // From Boscheri's paper

   double rotation_angle = 0.; // 0 - 1D horizontal velocity
   double rhoL = 1.0, rhoR = 1.0;;
   double rhoLstar = 3.9992502342988532;
   double rhoRstar = 3.9992502342988532;
   double vL = 2., vR = 0., vstar = 1;
   double pL = (this->get_gamma() - 1.) * pow(10., -4); // this is actually sie
   double pR= (this->get_gamma() - 1.) * pow(10., -4);  // also sie
   double pstar = 1.3334833281256511;
   double lambda1m = 0.66658333854101581;
   double lambda1p = 0.66658333854101581; 
   double lambda3 = 1.3334166614589844;

   double cL = sqrt(this->get_gamma() * pL / rhoL);
   double x0 = 0.0000000001; // Initial shock position

public:
   SaltzmannProblem()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_bcs_indicator(_bcs);
      this->set_mv_bcs_need_updating_indicator(_mv_bcs_need_updating);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);
      // CFL change
      this->set_cfl_change(_change_cfl);
      this->set_cfl_first(_cfl_first);
      this->set_cfl_second(_cfl_second);
      this->set_cfl_time_change(_cfl_time_change);
   }

   /* Override specific update functions */
   void lm_update(const double b_covolume) override 
   {
      this->set_b(b_covolume);
   }

   void get_additional_BCs(const FiniteElementSpace &fes, Array<int> ess_bdr, Array<int> &add_ess_tdofs, Array<double> &add_bdr_vals, const Geometric<dim> &geom=NULL) override 
   {
      std::cout << "saltzman::get_additional_BCs\n";

      Array<int> dofs_list;
      ess_bdr = 0;
      ess_bdr[4] = 1;

      /* Boundary conditions: This test enforces nonzero dirichlet conditions on left wall - marker 5 */
      fes.GetEssentialTrueDofs(ess_bdr, dofs_list, 0); // only x coordinate is nonzero
      add_ess_tdofs.Append(dofs_list);

      /* Fill bdr vals with 0 in this case */
      int size_add_bdr_dofs = add_ess_tdofs.Size();
      add_bdr_vals.SetSize(size_add_bdr_dofs);

      /* since current time is t=0, this val is 0. */
      add_bdr_vals = 0.;
   }

   /* The dirichlet condition that is enforced on the left hand side changes in time */
   void update_additional_BCs(const double &t, const double timestep_first, Array<double> &add_bdr_vals, const Geometric<dim> &geom=NULL, const ParGridFunction &x_gf=NULL) override 
   {
      assert(timestep_first > 0.);
      /* Solve for Dirichlet velocity at left wall */
      double _xi = t / (2*timestep_first);
      double _psi = (4. - (_xi + 1.) * (_xi - 2.) * ((_xi - 2.) - (abs(_xi-2.) + (_xi-2.)) / 2.)) / 4.;
      add_bdr_vals = _psi;
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double pressure(const double &rho, const double &sie, const int &cell_attr=0) override
   {
      return (this->get_gamma() - 1.) * rho * sie;
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double p0(const Vector &x, const double &t) override
   {
      if (t == 0) { return pow((this->get_gamma() - 1), 2) * 10e-4; } // TODO: Change pressure
      else {
         Vector _n(dim);
         _n[0] = cos(rotation_angle);
         if (dim > 1)
         {
            _n[1] = sin(rotation_angle);
         }
         double _x_tilde = x * _n;

         double _xi = (_x_tilde - x0) / t;

         // Saltzman TODO
         if (_x_tilde <= lambda1p * t)
         {
            return pL;
         }
         else if (_x_tilde <= lambda3 * t)
         {
            return pstar;
         }
         else
         {
            return pR;
         }
      }
   }

   double rho0(const Vector &x, const double & t) override
   {
      assert(dim < 3); //
      if (t == 0) { return 1.; }
      else 
      {
         Vector _n(dim);
         _n[0] = cos(rotation_angle);
         if (dim > 1)
         {
            _n[1] = sin(rotation_angle);
         }
         double _x_tilde = x * _n;

         double _xi = (_x_tilde - x0) / t;

         // Saltzman TODO
         if (_x_tilde <= lambda1p * t)
         {
            return rhoL;
         }
         else if (_x_tilde <= lambda3 * t)
         {
            return rhoRstar;
         }
         else
         {
            return rhoR;
         }
      }
      return 0.;
   }
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      if (t == 0) { v = 0.; }
      else 
      {
         Vector _n(dim);
         _n[0] = cos(rotation_angle);
         if (dim > 1)
         {
            _n[1] = sin(rotation_angle);
         }
         double _x_tilde = x * _n;

         double _xi = (_x_tilde - x0) / t;

         // our velocity will be some scaling of the normal vector
         v = _n;

         // Saltzman TODO
         if (_x_tilde <= lambda1p * t)
         {
            v *= vL;
         }
         else if (_x_tilde <= lambda3 * t)
         {
            v *= vstar;
         }
         else
         {
            v *= vR;
         }
      }
      return;
   }
   double sie0(const Vector &x, const double & t) override
   {
      if (t == 0) { return pow(10, -5); }
      else {
         double _p = p0(x,t);
         double _rho = rho0(x,t);
         return _p / (_rho * (this->get_gamma() - 1.));
      }
   }

}; // End class

} // ns hydroLO
} // ns mfem