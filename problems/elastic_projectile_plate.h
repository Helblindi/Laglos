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

template<int dim>
class ElasticProjectilePlate: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _gamma_g = 1.4, _gamma_s = 4.22;
   const int cell_attr_g = 1, cell_attr_s = 50;
   bool _known_exact_solution = false;
   bool _thbcs = false; // Indicator for thermal boundary conditions
   bool _mvbcs = false; // Indicator for mv boundary conditions
   string _indicator = "ElasticProjectilePlate";

   //https://www.sciencedirect.com/science/article/pii/S0021999109002654?fr=RR-2&ref=pdf_download&rr=928faaf93aca69c5
   // 5.2 elastic projectile plate section 7.1
   double rho_g = 1., rho_s = 8.9E3; // kg/m^3
   double v_proj = 800.; // m/s
   const double p_inf = 3.42E10; // Pa
   /* Different shear moduli for projectile plate */
   const double _mu = 9.2E10; // .002
   // const double _mu = 9.2E9; // .002 for bounceback
   // const double _mu = 9.2E8; // 0.002 re-stiffens 
   // const double _mu = 0.; // bound to crash


   /* helper function to determine the region based on x,y coords */
   bool is_solid_region(const Vector &x)
   {
      return (is_projectile_region(x) || is_plate_region(x));
   }
   bool is_projectile_region(const Vector &x)
   {
      double _x = x[0];
      double _y = x[1];
      /* projectile location */
      if (_x >= 0.2 && _x <= .3 && _y >= 0.45 && _y <= 0.55) { return true; }
      else { return false; }
   }
   bool is_plate_region(const Vector &x)
   {
      double _x = x[0];
      double _y = x[1];
      /* plate location */
      if (_x >= 0.3 && _x <= .4 && _y >= 0.25 && _y <= 0.75) { return true; }
      return false;
   }

public:
   ElasticProjectilePlate()
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
      this->set_exact_solution(_known_exact_solution);
   }

   double get_gamma(const int &cell_attr) override {
      /*
      cell_attr of 50 represents the solid material in the center 
      while 1 represents the gas outside
      */
      assert((cell_attr == cell_attr_s || cell_attr == cell_attr_g) && "Must pass in a cell_attr to any ProblemBase::get_gamma funcalls.\n");
      return (cell_attr == cell_attr_s) ? _gamma_s : _gamma_g;
   }
   
   /* Override specific update functions */
   void lm_update(const double b_covolume) override 
   {
      this->set_b(b_covolume);
   }

   void get_additional_BCs(const FiniteElementSpace &fes, Array<int> ess_bdr, Array<int> &add_ess_tdofs, Array<double> &add_bdr_vals, const Geometric<dim> &geom=NULL) override
   {
      Array<int> dofs_list;
      ess_bdr = 0;
      ess_bdr[4] = 1;

      /* 
      Boundary conditions: This test enforces dirichlet at left boundary on
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
   double pressure(const double &rho, const double &sie, const int &cell_attr=0) override
   {
	   return (this->get_gamma(cell_attr) - 1.) * rho * sie - this->get_gamma(cell_attr) * p_inf;
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double p0(const Vector &x, const double & t) override
   {
      if (t < 1e-12)
      {
         if (is_solid_region(x)) {
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
   double rho0(const Vector &x, const double & t) override
   {
      assert(t < 1e-12); // No exact solution
      if (is_solid_region(x)) {
         return rho_s;
      } else {
         return rho_g;
      }
   }
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      assert(t < 1e-12);
      v = 0.;
      if (is_projectile_region(x)) { v[0] = v_proj; }

      return;
   }

   double sie0(const Vector &x, const double & t) override
   {
      cout << "sie0\n";
      if (is_solid_region(x))
      {
         // return 0.;
         // double _gam = this->get_gamma();
         return (p0(x,t) + p_inf * this->get_gamma(50)) / this->rho0(x, t) / (this->get_gamma(50) - 1.0);
         // return this->rho0(x, t) * this->sie0(x,t) / (_gam - 1.0) - p_inf * _gam;
      }
      else
      {
         // Value corresponds to 1 atm = 101.325 kPa
         // return 2.5331125E5; // J / kg
         const double _gam = _gamma_g;
         return p0(x,t) * (_gam - 1.) / rho0(x,t);
      }
   }
}; // End class

} // ns hydroLO
} // ns mfem
