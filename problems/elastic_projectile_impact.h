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
#include "riemann1D.hpp"
#include <cmath>

using namespace mfem;
using namespace std;

namespace mfem
{
namespace hydroLO
{

class ElasticProjectileImpact: public ProblemBase
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   // double _gamma_g = 1.4, _gamma_s = 4.22;
   double _gamma_g = 1.4, _gamma_s = 3.4; // aortic
   const int cell_attr_g = 1, cell_attr_s = 50;
   bool _known_exact_solution = false;
   bool _thbcs = false; // Indicator for thermal boundary conditions
   bool _mvbcs = true; // Indicator for mv boundary conditions
   bool _distort_mesh = false;
   string _indicator = "ElasticProjectileImpact";

   //https://www.sciencedirect.com/science/article/pii/S0021999109002654?fr=RR-2&ref=pdf_download&rr=928faaf93aca69c5
   // Elastic projectile impact - See vilar-shu-maire-2D
   double rho_g = 1., rho_s = 2.7E3; // kg/m^3, aortic
   double v_proj = -150.; // m/s
   const double p_inf = 2.15E10; // Pa, aortic 
   const double _mu = 2.6E10; // For aortic model



   /* helper function to determine the region based on x,y coords */
   bool is_projectile_region(const Vector &x) const
   {
      double _x = x[0];
      double _y = x[1];
      /* projectile location */
      if (_x >= 0. && _x <= 5. && _y >= 0. && _y <= 1.) { return true; }
      else { 
         cout << "_x: " << _x << ", _y: " << _y << endl;
         return false; }
   }

public:
   ElasticProjectileImpact(const int &_dim) : ProblemBase(_dim)
   {
      if (dim != 2)
      {
         /*
         I do not know how to implement 2d velocity for a 1d problem,
         hence we must be compiled in 2D.
         */
         MFEM_ABORT("Dimension must be 2 for the elastic shear problem.\n");
      }
      this->set_indicator(_indicator);
      this->set_pinf(p_inf);
      this->set_shear_modulus(_mu);
      this->set_thbcs_indicator(_thbcs);
      this->set_mvbcs_indicator(_mvbcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);

      // Set Equation of state
      this->eos = std::unique_ptr<EquationOfState>(new NobleAbelStiffenedGasEOS(this->b,this->q,this->p_inf));
   }

   double get_gamma(const int &cell_attr) const override {
      /*
      cell_attr of 50 represents the solid material in the center 
      while 1 represents the gas outside
      */
      assert((cell_attr == cell_attr_s || cell_attr == cell_attr_g) && "Must pass in a cell_attr to any ProblemBase::get_gamma funcalls.\n");
      return (cell_attr == cell_attr_s) ? _gamma_s : _gamma_g;
   }

   void get_additional_BCs(const FiniteElementSpace &fes, Array<int> ess_bdr, Array<int> &add_ess_tdofs, Array<double> &add_bdr_vals, const Geometric *geom=NULL) override
   {
      Array<int> dofs_list;
      ess_bdr = 0;
      ess_bdr[4] = 1;

      /* 
      Boundary conditions: This test enforces dirichlet at right boundary on
      the x-coord only. [marker 5] 
      The remaining boundary is treated with do-nothing bcs. 
      */
      int d = 0;
      fes.GetEssentialTrueDofs(ess_bdr, dofs_list, d);
      add_ess_tdofs.Append(dofs_list);

      /* remove possible duplicates */
      add_ess_tdofs.Unique();

      /* Fill bdr vals with 0 in this case */
      add_bdr_vals.SetSize(add_ess_tdofs.Size());
      add_bdr_vals = v_proj;

      assert(add_bdr_vals.Size() == add_ess_tdofs.Size());
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double p0(const Vector &x, const double & t) const override
   {
      if (t < 1e-12)
      {
         if (is_projectile_region(x)) {
            // const double _gam = _gamma_s;
            // return this->rho0(x, t) * this->sie0(x,t) / (_gam - 1.0) - p_inf * _gam;
            return 0.;

         } else
         {
            return 101325.; // 1 atm = 101.325 kPa
         }
      }
      else
      {
         MFEM_ABORT("No exact solution.\n");
         return -1.;
      }
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double rho0(const Vector &x, const double & t) const override
   {
      assert(t < 1e-12); // No exact solution
      if (is_projectile_region(x)) {
         return rho_s;
      } else {
         return rho_g;
      }
   }
   void v0(const Vector &x, const double & t, Vector &v) const override
   {
      assert(t < 1e-12);
      v = 0.;
      if (is_projectile_region(x)) { v[0] = v_proj; }

      return;
   }

   double sie0(const Vector &x, const double & t) const override
   {
      cout << "sie0\n";
      double _rho = rho0(x,t);
      double _p = p0(x,t);
      if (is_projectile_region(x))
      {
         return this->eos->energy(_p, _rho, _gamma_s);
      }
      else
      {
         // Value corresponds to 1 atm = 101.325 kPa
         // return 2.5331125E5; // J / kg
         MFEM_ABORT("Not successfully implemented");
         // Need to look at how the pressure is also computed.
         return _p * (_gamma_g - 1.) / _rho;
      }
   }
}; // End class

} // ns hydroLO
} // ns mfem
