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

class KidderProblem: public ProblemBase
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   const double _a = 0., _b = 0., _gamma = 2.;
   const double r1 = 0.9, r2 = 1.0, rho1 = 1., rho2 = 2., s = 1.;
   const double P1 = pow(rho1, _gamma), P2 = pow(rho2, _gamma);
   const double tau;
   double t_exact;
   bool _distort_mesh = false;
   bool _known_exact_solution = true;
   bool _thbcs = false; // Indicator for thermal boundary conditions
   bool _mvbcs = false; // Indicator for mv boundary conditions
   bool _mv_bcs_need_updating = true;
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
   KidderProblem(const int &_dim) : ProblemBase(_dim), tau(compute_focusing_time()), t_exact(0.)
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

   /*********************************************************
    * Initial State functions
    * NOTE: the exact solution for pressure, density, and 
    * velocity are in terms of initial particle locations.
    * Hence these functions assume that the first parameter 
    * is the initial particle location, rather than the 
    * particle location at time t.
    *********************************************************/
   double p0(const Vector &x, const double & t) const override
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
   double rho0(const Vector &x, const double & t) const override
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
   void v0(const Vector &x, const double & t, Vector &v) const override
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
   double sie0(const Vector &x, const double & t) const override
   {
      double _p = p0(x,t);
      double _rho = rho0(x,t);
      return this->eos->energy(_p, _rho, this->get_gamma());
   }

   /*********************************************************
    * Problem specific functions
    *********************************************************/
   double h(const double t) const
   {
      return sqrt(1. - pow(t/tau, 2));
   }
   double dhdt(const double t) const
   {
      return t / pow(tau, 2) / h(t);     
   }

   void GetBoundaryState(const Vector &x, const int &bdr_attr, Vector &state, const double &t) override
   { 
      double _t = 0.;
      if (t == 0.) {
         _t = t_exact;
      } else {
         _t = t;
      }

      state.SetSize(dim + 2);
      state = 0.;

      Vector vel(dim);
      double norm_x = x.Norml2();
      vel = x;
      vel *= -1.;
      vel /= norm_x;
      double _dhdt = dhdt(_t);
      double hmult = pow(h(_t), 2. * this->get_gamma() / (1. - this->get_gamma()));
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

   void get_additional_BCs(const FiniteElementSpace &fes, Array<int> ess_bdr, Array<int> &add_ess_tdofs, Array<double> &add_bdr_vals, const Geometric *geom=NULL) override { }
   void update_additional_BCs(const double &t, const double timestep_first, Array<double> &add_bdr_vals, const Geometric *geom=NULL, const ParGridFunction *x_gf=NULL) override
   {
      // Need to update the current time since the function in LaglosSolver that calls GetBoundaryState does not have access to the current time
      t_exact = t;
   }
}; // End class

} // ns hydroLO
} // ns mfem