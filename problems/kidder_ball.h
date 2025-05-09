/*********************************************************
* KidderBall Problem
*
* See Burton, Morgan, Carney, Kenamond 
* Reduction of dissipation in Lagrange cell-centered hydrodynamics 
* (CCH) through corner gradient reconstruction(CGR)
*
* Notes:
* 1) Exact velocity is enforced on all boundary nodes.
* 2) The exact solution is enforced on all boundary cells.
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

class KidderBallProblem: public ProblemBase
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0., _b = 0., _gamma = 5./3.;
   bool _distort_mesh = false;
   bool _known_exact_solution = true;
   bool _thbcs = true; // Indicator for thermal boundary conditions
   bool _mvbcs = true; // Indicator for mv boundary conditions
   bool _mv_bcs_need_updating = true;
   int size_add_bdr_dofs;
   Array<int> bdr_dofs_list;
   string _indicator = "KidderBall"; // Possible: saltzmann

public:
   KidderBallProblem(const int &_dim) : ProblemBase(_dim)
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_thbcs_indicator(_thbcs);
      this->set_mvbcs_indicator(_mvbcs);
      this->set_mv_bcs_need_updating_indicator(_mv_bcs_need_updating);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);

      // Set Equation of state
      this->eos = std::unique_ptr<EquationOfState>(new IdealGasEOS());
   }

   /* Optionally overridden, or removed */
   double get_gamma(const int &cell_attr = 0) const override { return _gamma; }
   void lm_update(const double b_covolume) override {}
   void update(Vector vec, double t = 0.) override {}

   void get_additional_BCs(const FiniteElementSpace &fes, Array<int> ess_bdr, Array<int> &add_ess_tdofs, Array<double> &add_bdr_vals, const Geometric *geom=NULL) override 
   {
      std::cout << "KidderBall::get_additional_BCs\n";

      ess_bdr = 0;
      ess_bdr[4] = 1;

      /* Boundary conditions: This test enforces nonzero dirichlet conditions on all boundary vertices */
      fes.GetEssentialTrueDofs(ess_bdr, bdr_dofs_list);
      add_ess_tdofs.Append(bdr_dofs_list);
      size_add_bdr_dofs = add_ess_tdofs.Size();
      add_bdr_vals.SetSize(size_add_bdr_dofs);

      // Let value of exact solution be set by first call to update BCs.
      assert(add_bdr_vals.Size() == size_add_bdr_dofs);
   }

   /* The dirichlet condition that is enforced on the left hand side changes in time */
   void update_additional_BCs(const double &t, const double timestep_first, Array<double> &add_bdr_vals, const Geometric *geom=NULL, const ParGridFunction *x_gf=NULL) override 
   {
      /* Validate we do not divide by 0 and than the array of dofs is the right size */
      assert(add_bdr_vals.Size() == size_add_bdr_dofs);
      Vector x(dim);
      Vector v(dim);
      /* Solve for exact velocity at the boundary points */
      for (int bdr_v_it = 0; bdr_v_it < size_add_bdr_dofs / dim; bdr_v_it++)
      {
         int dof = bdr_dofs_list[bdr_v_it];
         geom->GetNodePosition(*x_gf, dof, x);
         this->v0(x, t, v);
         add_bdr_vals[bdr_v_it] = v[0];
         add_bdr_vals[bdr_v_it + size_add_bdr_dofs / dim] = v[1];
      }
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double p0(const Vector &x, const double & t) const override
   {
      double _rho = rho0(x,t);
      double _sie = sie0(x,t);
      return this->eos->pressure(_rho, _sie, this->get_gamma());
   }
   double rho0(const Vector &x, const double & t) const override
   {
      double _r = x.Norml2();
      double val = -1. * pow(_r,2) / (1. + pow(t-1.,2));
      val = exp(val);
      val *= 2. * pow(1 + pow(t-1.,2), -1.5);
      return val;
   }
   void v0(const Vector &x, const double & t, Vector &v) const override
   {
      double _r = x.Norml2();
      double val = (t-1.) / (1. + pow(t-1.,2));
      v = x;
      v *= val;
      
      return;
   }
   double sie0(const Vector &x, const double & t) const override
   {
      double val = 0.75 / (1. + pow(t-1,2));
      return val;
   }

}; // End class

} // ns hydroLO
} // ns mfem