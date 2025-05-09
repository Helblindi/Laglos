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

class IsentropicVortex: public ProblemBase
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 0.;     
   double _b = 0.;
   double _gamma = 7./5.;
   bool _distort_mesh = false;
   bool _known_exact_solution = true;
   bool _thbcs = true;
   bool _mvbcs = true;
   string _indicator = "IsentropicVortex"; // Possible: saltzmann

   // Free steam conditions
   const double rho_inf = 1., p_inf = 1., T_inf = 1.;
   const double beta = 5.;
   double xc_0 = 0., xc_1 = 0.;
   double vc_0 = 0., vc_1 = 0.; // u_inf

public:
   IsentropicVortex(const int &_dim) : ProblemBase(_dim)
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
      this->eos = std::unique_ptr<EquationOfState>(new PolytropicEOS(1.));
   }

   /* Override specific update functions */
   void lm_update(const double b_covolume) override 
   {
      this->set_b(b_covolume);
   }

   void get_additional_BCs(const FiniteElementSpace &fes, Array<int> ess_bdr, Array<int> &add_ess_tdofs, Array<double> &add_bdr_vals, const Geometric *geom=NULL) override
   {
      Array<int> dofs_list;
      ess_bdr = 0;
      ess_bdr[4] = 1;

      /* Boundary conditions: This test freezes all nodes on the boundaries - marker 5 */
      for (int d = 0; d < dim; d++)
      {
         fes.GetEssentialTrueDofs(ess_bdr, dofs_list, d);
         add_ess_tdofs.Append(dofs_list);
      }
      /* remove possible duplicates */
      add_ess_tdofs.Unique();

      /* Fill bdr vals with 0 in this case */
      add_bdr_vals.SetSize(add_ess_tdofs.Size());
      add_bdr_vals = 0.;

      assert(add_bdr_vals.Size() == add_ess_tdofs.Size());
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double p0(const Vector &x, const double &t) const override
   {
      double _rho = rho0(x,t);
      double _sie = sie0(x,t);
      return this->eos->pressure(_rho, _sie, this->get_gamma());
   }
   double rho0(const Vector &x, const double & t) const override
   {
      Vector center(2);
      center[0] = xc_0 + vc_0 * t;
      center[1] = xc_1 + vc_1 * t;

      Vector x_bar = x;
      x_bar -= center;
      const double r = x_bar.Norml2();

      const double dT = -1. * (exp(1 - pow(r,2)) * pow(beta, 2) * (this->get_gamma() - 1.)) / (8. * this->get_gamma() * pow(M_PI, 2));
      const double rho = pow(T_inf + dT, 1./(this->get_gamma() - 1.));

      return rho;
   }
   void v0(const Vector &x, const double & t, Vector &v) const override
   {
      Vector center(2);
      center[0] = xc_0 + vc_0 * t;
      center[1] = xc_1 + vc_1 * t;

      Vector x_bar = x;
      x_bar -= center;
      const double r = x_bar.Norml2();

      const double coeff = (exp((1. - pow(r,2)) / 2.) * beta) / (2. * M_PI);

      v[0] = vc_0 - coeff * x_bar[1];
      v[1] = vc_1 + coeff * x_bar[0];
      
      return;
   }
   double sie0(const Vector &x, const double & t) const override
   {
      double _p = p0(x,t);
      double _rho = rho0(x,t);
      return this->eos->energy(_p, _rho, this->get_gamma());
   }

}; // End class

} // ns hydroLO
} // ns mfem