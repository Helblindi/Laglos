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
class ElasticIsentropicVortex: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   const double _gamma_g = 1.4, _gamma_s = 3.4;
   const int cell_attr_g = 1, cell_attr_s = 50;
   const double solid_multiplier = 1.E3;
   const double p_inf_s = 2.15E10;
   const double _mu = 2.6E10;
   /* options */
   bool _distort_mesh = false;
   bool _known_exact_solution = false;
   bool _thbcs = true;
   bool _mvbcs = true;
   string _indicator = "ElasticIsentropicVortex"; // Possible: saltzmann

   // Free steam conditions
   const double rho_inf = 1., p_inf = 1., T_inf = 1.;
   const double beta = 5.;
   double xc_0 = 0., xc_1 = 0.;
   double vc_0 = 0., vc_1 = 0.; // u_inf

   /* helper function to determine the region based on x,y coords */
   bool is_solid_region(const Vector &x)
   {
      if (abs(x[0]) <= 2.5 && abs(x[1]) <= 2.5) { return true; }
      return false;
   }

public:
ElasticIsentropicVortex()
   {
      this->set_indicator(_indicator);
      this->set_pinf(p_inf);
      this->set_shear_modulus(_mu);
      this->set_thbcs_indicator(_thbcs);
      this->set_mvbcs_indicator(_mvbcs);
      this->set_distort_mesh(_distort_mesh);
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
    * Problem Description functions
    *********************************************************/
   double pressure(const double &rho, const double &sie, const int &cell_attr=0) override
   {
      if (cell_attr == cell_attr_s) {
         return rho*sie*(this->get_gamma(cell_attr) - 1.) - p_inf_s * this->get_gamma(cell_attr);
      } else 
      {
         assert(cell_attr == cell_attr_g);;
         return pow(rho, this->get_gamma(cell_attr));
      }  
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double p0(const Vector &x, const double &t) override
   {
      cout << "p0\n";
      double _gam = _gamma_g;
      double rho = rho0(x,t);
      if (is_solid_region(x)) { 
         cout << "solid\n";
         _gam = _gamma_s; 
         // double sie = sie0(x,t);
         // return rho*sie*(_gam - 1.) - p_inf_s *  _gam;
      } 
      // else {
         // cout << "rho: " << rho << ", gamma: " << _gam << endl;
      // }
      return pow(rho, _gam);
   }
   double rho0(const Vector &x, const double & t) override
   {
      cout << "rho0\n";
      double _gam = _gamma_g, _solid_mult = 1.;
      if (is_solid_region(x))
      {
         _gam = _gamma_s;
         _solid_mult = solid_multiplier;
      }

      Vector center(2);
      center[0] = xc_0 + vc_0 * t;
      center[1] = xc_1 + vc_1 * t;

      Vector x_bar = x;
      x_bar -= center;
      const double r = x_bar.Norml2();

      const double dT = -1. * (exp(1 - pow(r,2)) * pow(beta, 2) * (_gam - 1.)) / (8. * _gam * pow(M_PI, 2));
      double rho = pow(T_inf + dT, 1./(_gam - 1.));
      rho *= _solid_mult;

      return rho;
   }
   void v0(const Vector &x, const double & t, Vector &v) override
   {
      cout << "v0\n";
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
   double sie0(const Vector &x, const double & t) override
   {
      cout << "sie0\n";
      if (is_solid_region(x)) {
         const double _gam = _gamma_s;
         return (p0(x,t) + p_inf * _gam) / this->rho0(x, t) / (_gam - 1.0);
      } else
      {
         const double _gam = _gamma_g;
         return p0(x,t) / this->rho0(x, t) / (_gam - 1.0);
      }
   }
}; // End class

} // ns hydroLO
} // ns mfem