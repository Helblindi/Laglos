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
#include "shear_closure.hpp"
#include <cmath>

using namespace mfem;
using namespace std;

namespace mfem
{
namespace hydroLO
{

class ElasticHumanLegImpact: public ProblemBase
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _gamma = 4.2;
   bool _known_exact_solution = false;
   bool _thbcs = false; // Indicator for thermal boundary conditions
   bool _mvbcs = false; // Indicator for mv boundary conditions
   string _indicator = "ElasticHumanLegImpact";

   //https://www.sciencedirect.com/science/article/pii/S0021999109002654?fr=RR-2&ref=pdf_download&rr=928faaf93aca69c5
   // 5.2 elastic projectile plate section 7.1
   double rho_bone = 1986., rho_flesh = 1000., rho_bar = 2400.; // kg/m^3
   double v_proj = -13.4112; // m/s, equivalent to 30 mph impact speed
   const double p_inf = 2.39E10; // Pa
   const double _mu = 1.3E10; // Pa
   /* Different shear moduli */
   // const double p_inf = 3.42E8; // Pa
   // const double _mu = 9.2E8; // Pa
   /* Bone and flesh params */
   double E_bone = 1.89E10, EA_bone = 1.029E10;
   double GA_bone = 5.63E9, nu_bone = 0.312;
   double E_flesh = 1.729E6, EA_flesh = 1.E6;
   double GA_flesh = 4.59E5, nu_flesh = 0.4;
   double _theta = M_PI/2.;
   ShearClosure *shear_closure_bone = nullptr;
   ShearClosure *shear_closure_flesh = nullptr;
   ShearClosure *shear_closure_bar = nullptr;


   /* helper function to determine the region based on x,y coords */
   bool is_bone_region(const Vector &x) const
   {
      double _x = x[0];
      double _y = x[1];
      if (_x >= 0. && _x <= 0.03 && _y >= 0.0 && _y <= 0.24) { return true; }
      else { return false; }
   }
   bool is_flesh_region(const Vector &x) const
   {
      double _x = x[0];
      double _y = x[1];
      if (_x >= 0.03 && _x <= 0.09 && _y >= 0.0 && _y <= 0.24) { return true; }
      else { return false; }
   }
   bool is_bar_region(const Vector &x) const
   {
      double _x = x[0];
      double _y = x[1];
      if (_x >= 0.09 && _x <= 0.19 && _y >= 0.07 && _y <= 0.17) { return true; }
      else { return false; }
   }

public:
   ElasticHumanLegImpact(const int &_dim) : ProblemBase(_dim)
   {
      if (dim != 2)
      {
         MFEM_ABORT("Dimension must be 2 for the elastic human leg impact problem.\n");
      }
      this->set_indicator(_indicator);
      this->set_pinf(p_inf);
      this->set_shear_modulus(_mu);
      this->set_thbcs_indicator(_thbcs);
      this->set_mvbcs_indicator(_mvbcs);
      this->set_exact_solution(_known_exact_solution);
      this->set_gamma(_gamma);

      // Set Equation of state
      this->eos = std::unique_ptr<EquationOfState>(new NobleAbelStiffenedGasEOS(this->b,this->q,this->p_inf));

      // Overwrite mi vector for anisotropic models
      theta = _theta;
      mi_vec.SetSize(3);
      mi_vec(0) = cos(theta);
      mi_vec(1) = sin(theta);
      mi_vec(2) = 0.;
      // Set elastic models
      shear_closure_bone = new ShearClosureTransverselyIsotropic(mu, mi_vec, E_bone, EA_bone, GA_bone, nu_bone);
      shear_closure_flesh = new ShearClosureTransverselyIsotropic(mu, mi_vec, E_flesh, EA_flesh, GA_flesh, nu_flesh);
      // shear_closure_bone = new ShearClosureNeoHookean(mu);
      // shear_closure_flesh = new ShearClosureNeoHookean(mu);
      shear_closure_bar = new ShearClosureNeoHookean(mu);
   }
   
   /* Override specific update functions */
   void lm_update(const double b_covolume) override 
   {
      this->set_b(b_covolume);
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   double p0(const Vector &x, const double & t) const override
   {
      if (t < 1e-12)
      {
         return 0.;
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
      if (is_bone_region(x)) {
         return rho_bone;
      } else if (is_flesh_region(x)) {
         return rho_flesh;
      } else if (is_bar_region(x)) {
         return rho_bar;
      } else {
         MFEM_ABORT("Point not in any defined region.\n");
         return -1.;
      }
   }
   void v0(const Vector &x, const double & t, Vector &v) const override
   {
      assert(t < 1e-12);
      v = 0.;
      if (is_bar_region(x)) { v[0] = v_proj; }

      return;
   }

   double sie0(const Vector &x, const double & t) const override
   {
      double _rho = rho0(x,t);
      double _p = p0(x,t);
      return this->eos->energy(_p, _rho, _gamma);
   }

   /*********************************************************
   * Elastic functions
   *********************************************************/
   double e_shear(const int &e, const int &cell_attr) const override
   {
      const double rho0 = rho0_v(e);
      DenseMatrix F(dim);
      elastic->ComputeAvgF(e,F);
      switch (cell_attr)
      {
         case 51: // bone
            return shear_closure_bone->ComputeShearEnergy(F, rho0);
         case 52: // flesh
            return shear_closure_flesh->ComputeShearEnergy(F, rho0);
         case 53: // bar
            return shear_closure_bar->ComputeShearEnergy(F, rho0);
         default:
         {
            cout << "cell_attr = " << cell_attr << "\n";
            MFEM_ABORT("Invalid cell attribute in ElasticHumanLegImpact::e_shear");
         }
      }
      return 0.;
   }

   void ComputeS(const int &e, const double &rho, DenseMatrix &S, const int &cell_attr) const override
   {
      DenseMatrix F(3);
      elastic->ComputeAvgF(e, F);
      const double rho0 = rho0_v(e);
      switch (cell_attr)
      {
         case 51: // bone
            shear_closure_bone->ComputeCauchyStress(F, rho, rho0, S);
            return;
         case 52: // flesh
            shear_closure_flesh->ComputeCauchyStress(F, rho, rho0, S);
            return;
         case 53: // bar
            shear_closure_bar->ComputeCauchyStress(F, rho, rho0, S);
            return;
         default:
         {
            cout << "cell_attr = " << cell_attr << "\n";
            MFEM_ABORT("Invalid cell attribute in ElasticHumanLegImpact::ComputeS");
         }
      }
      return;
   }

}; // End class

} // ns hydroLO
} // ns mfem
