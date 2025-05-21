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

#include "eos.h"
#include "mfem.hpp"
#include "problem_base.h"
#include <cmath>

using namespace mfem;
using namespace std;

namespace mfem
{
namespace hydroLO
{

class VDWIsentropicVortex: public ProblemBase
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 1.;     
   double _b = 0.;
   double _gamma = 1.5;
   bool _distort_mesh = false;
   bool _known_exact_solution = true;
   bool _thbcs = true;
   bool _mvbcs = true;
   string _indicator = "VDWIsentropicVortex"; // Possible: saltzmann

   // Free steam conditions
   const double rho_inf = 0.1, p_inf = 1.;
   const double beta = 20.;
   Vector x0, vinf;
   const double r0 = 1.0;

   /* Do not change */
   // double C = (p_inf + _a * rho_inf * rho_inf) / std::pow(rho_inf, 3./2.);
   static constexpr double F = -30.1;
   static constexpr double C = 31.9390; // 101. / sqrt(10);

public:
   VDWIsentropicVortex(const int &_dim) : x0(_dim), vinf(_dim), ProblemBase(_dim)
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_thbcs_indicator(_thbcs);
      this->set_mvbcs_indicator(_mvbcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);
      
      /* Set x0 and vinf */
      x0 = 0.;
      vinf = 0.;

      // Set Equation of state
      this->eos = std::unique_ptr<EquationOfState>(new VanDerWaalsEOS(_a,_b));
   }

   // void get_additional_BCs(const FiniteElementSpace &fes, Array<int> ess_bdr, Array<int> &add_ess_tdofs, Array<double> &add_bdr_vals, const Geometric *geom=NULL) override
   // {
   //    Array<int> dofs_list;
   //    ess_bdr = 0;
   //    ess_bdr[4] = 1;

   //    /* Boundary conditions: This test freezes all nodes on the boundaries - marker 5 */
   //    for (int d = 0; d < dim; d++)
   //    {
   //       fes.GetEssentialTrueDofs(ess_bdr, dofs_list, d);
   //       add_ess_tdofs.Append(dofs_list);
   //    }
   //    /* remove possible duplicates */
   //    add_ess_tdofs.Unique();

   //    /* Fill bdr vals with 0 in this case */
   //    add_bdr_vals.SetSize(add_ess_tdofs.Size());
   //    add_bdr_vals = 0.;

   //    assert(add_bdr_vals.Size() == add_ess_tdofs.Size());
   // }
   double psi(const Vector &x) const
   {
      double _r02 = r0 * r0;
      double _val = x * x;
      _val *= -1./ _r02;
      _val += 1.;
      _val *= 0.5;
      _val = std::exp(_val);
      _val *= beta / 2. / M_PI;
      return _val;
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double p0(const Vector &x, const double &t) const override
   {
      double _rho = rho0(x,t);
      return C * std::pow(_rho, this->get_gamma()) - this->get_a() * _rho * _rho;
   }
   double rho0(const Vector &x, const double & t) const override
   {
      Vector xbar(2);
      subtract(x, x0, xbar);
      xbar.Add(-1.*t, vinf);
      // MFEM_WARNING("CHECK ME");

      double _psi = psi(xbar);
      double _val = _psi * _psi / 2. / r0 / r0 + F;
      _val *= 2./this->get_a();
      _val += 9.*C*C / 16. / this->get_a() / this->get_a();
      _val = sqrt(_val);
      _val *= -0.5;
      _val += (3. * C) / 8. / this->get_a();
      return _val;
   }
   void v0(const Vector &x, const double & t, Vector &v) const override
   {
      Vector xbar(2);
      subtract(x, x0, xbar);
      xbar.Add(-1.*t, vinf);
      // MFEM_WARNING("CHECK ME");

      double _psi = psi(xbar);

      v.SetSize(2);
      v = vinf;
      v[0] -= _psi * xbar[1];
      v[1] += _psi * xbar[0];

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