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

extern "C" {
   // void __vdw_MOD_initialize_vdw(double *rhop, double *in_state, double *in_data, double *out_state);
   void __vdw_MOD_initialize_vdw(double *rhop, double *in_a, double *in_b, double *in_gamma, 
                                 double *in_rhoL, double *in_rhoR, double *out_vL, 
                                 double *out_vR, double *out_pL, double *out_pR);
   void __vdw_MOD_rho_v_p_vdw(double *x0, int *nmax, double *xx, double *rho, double *v, double *p);
   void __vdw_MOD_f_add(double *x, double *y, double *z);
   void __vdw_MOD_f_test();
}

namespace mfem
{
namespace hydroLO
{

template<int dim>
class VdwTest1: public ProblemBase<dim>
{
private:
   /*********************************************************
    * Problem Specific constants
    *********************************************************/
   double _a = 1., _b = 1., _gamma = 1.02;
   bool _distort_mesh = false;
   bool _known_exact_solution = true; // exact is known, information in supplementary material
   bool _bcs = false; // Indicator for boundary conditions
   string _indicator = "Vdw1"; // Possible: saltzmann

   double rhop = 0.35, x0 = 0.;
   double initial_shock = 0.;
   double rhoL = 0.1, rhoR = 0.39;
   double vL = -0.475504638574729, vR = -0.121375781741349;
   double pL = 0.022084258693080, pR = 0.039073167077590;
   double sieL = 14.337916411885988, sieR = 14.560722040683306;
   double *x_gf_sorted, *rho_d, *v_d, *p_d;
   Vector x_sorted_vec;
   int vec_size = 0;

public:
   VdwTest1()
   {
      this->set_a(_a);
      this->set_b(_b);
      this->set_gamma(_gamma);
      this->set_indicator(_indicator);
      this->set_bcs_indicator(_bcs);
      this->set_distort_mesh(_distort_mesh);
      this->set_exact_solution(_known_exact_solution);
   }

   /* Override getters */
   virtual void update(Vector x_gf, double t = 0.) override {
      // cout << "Vdw1::update()\n";
      // cout << "x_gf in update: " << endl;
      // x_gf.Print(cout);
      // x_gf /= t;
      // cout << "t: " << t << endl;
      // cout << "post division x_gf:\n";
      // x_gf.Print(cout);
      compute_vdw_arrays(x_gf);
   }

   /*********************************************************
    * Problem Description functions
    *********************************************************/
   virtual double pressure(const double &rho, const double &sie, const int &cell_attr=0) override
   {
       // Use van der Waals
      double val = (this->get_gamma() - 1.) * (rho * sie + this->get_a() * pow(rho, 2)) / (1. - this->get_b() * rho) - this->get_a() * pow(rho,2);
      return val;
   }

   /*********************************************************
    * Initial State functions
    *********************************************************/
   double p0(const Vector &x, const double &t) override
   {
      if (t < 1.e-16) {
         if (x[0] <= initial_shock)
         {
            return pL;
         }
         else
         {
            return pR;
         }
      }
      else {
         double rho = rho0(x,t);
         double val = (this->get_gamma() - 1.) * (rho * sie0(x,t) + this->get_a() * pow(rho, 2)) / (1. - this->get_b() * rho) - this->get_a() * pow(rho,2);
         return val;
      }
   }

   virtual double rho0(const Vector &x, const double & t) override
   {
      if (t < 1.e-16) {
         if (x[0] <= initial_shock)
         {
            return rhoL;
         }
         else
         {
            assert(x[0] <= 1.);
            return rhoR;
         }
      }
      else {
         int index;
         getIndex(x[0], index);
         return rho_d[index];
      }
   }
   virtual void v0(const Vector &x, const double & t, Vector &v) override
   {
      if (t < 1.e-16) {
         if (x[0] <= initial_shock)
         {
            v[0] = vL;
            return;
         }
         else
         {
            assert(x[0]<=1.);
            v[0] = vR;
            return;
         }
      }
      else {
         int index;
         getIndex(x[0], index);
         v[0] = v_d[index];
      }
   }
   virtual double sie0(const Vector &x, const double & t) override
   {
      if (t < 1.e-16) {
         if (x[0] <= initial_shock)
         {
            return sieL;
         }
         else
         {
            assert(x[0]<=1.);
            return sieR;
         }
      }
      else {
         return 14.; // TODO: Exact representation of sie0
      }
   }

   void compute_vdw_arrays(Vector x_gf)
   {
      // cout << "Vdw1::compute_vdw_arrays()\n";
      // cout << "x_gf in compute_vdw_arrays: " << endl;
      // x_gf.Print(cout);

      // First sort x_gf
      this->vec_size = x_gf.Size();
      x_gf_sorted = x_gf.GetData();
      x_sorted_vec.SetSize(vec_size);

      std::sort(x_gf_sorted, x_gf_sorted + vec_size);

      // Fill vector as pointer will fall out of scope
      for (int i = 0; i < vec_size; i++)
      {
         x_sorted_vec[i] = x_gf_sorted[i];
         // cout << "x: " << x_gf_sorted[i] << endl;
      }
      
      Vector rho(vec_size), v(vec_size), p(vec_size);
      double out_vL, out_vR, out_pL, out_pR;
      
      // Run Fortran code to compute exact solution
      __vdw_MOD_initialize_vdw(&rhop, &_a, &_b, &_gamma, &rhoL, &rhoR, &out_vL, &out_vR, &out_pL, &out_pR);
      __vdw_MOD_rho_v_p_vdw(&x0, &vec_size, x_gf_sorted, rho.GetData(), v.GetData(), p.GetData());

      // Stuff Fortran results into class to be accessed
      rho_d = rho.GetData();
      v_d = rho.GetData();
      p_d = rho.GetData();

      // cout << "Done computing vdw arrays\n";
   }

   void getIndex(const double val, int & index)
   {
      for (int j = 0; j < vec_size; j++)
      {
         if (x_sorted_vec[j] == val)
         {
            
            index = j;
            return;
         }
      }
      // cout << "Never found index for val: " << val << endl;
   }

// private:
   // double rk_rho(const double & _rho, const double & _S)
   // {
   //    double num = _S * this->get_gamma() * pow(_rho, this->get_gamma() - 1.);
   //    num *= (this->get_gamma() + 1.) * pow(1. - this->get_b() * _rho, -2. - this->get_gamma());
   //    num -= 6. * this->get_a() * _rho;

   //    double denom = pow(1 - this->get_b() * _rho, -1. - this->get_gamma());
   //    denom *= (_S * this->get_gamma() * pow(_rho, this->get_gamma()) - 2 * this->get_a() * pow(_rho,2) * pow(1.-this->get_b() * _rho, this->get_gamma() + 1.)) * _rho;

   //    return 2. * sqrt(denom) / num;
   // }

   // // In terms of vector x, how can change to just in terms of double x
   // void sol_rho(double & _rhoz, 
   //              double & _Sz, 
   //              double (*phi)(const double &, const double &), // &rk_rho
   //              Vector & _xx,
   //              Vector & _rho)
   // {
   //    double dx = _xx[0] - 0.;
   //    double k1 = phi(_rhoz, _Sz);
   //    double k2 = phi(_rhoz+dx*k1/2.,_Sz);
   //    double k3 = phi(_rhoz+dx*k2/2.,_Sz);
   //    double k4 = phi(_rhoz+dx*k3,_Sz);
   //    _rho[0] = _rhoz + (dx/6.)*(k1+2.*k2+2.*k3+k4);

   //    for (int i = 0; i < _xx.Size() - 1; i++)
   //    {
   //       dx = _xx[i+1] - _xx[i];
   //       k1 = phi(_rho[i], _Sz);
   //       k2 = phi(_rho[i]+dx*k1/2.,_Sz);
   //       k3 = phi(_rho[i]+dx*k2/2.,_Sz);
   //       k4 = phi(_rho[i]+dx*k3,_Sz);
   //       _rho[i+1] = _rho[i] + (dx/6.)*(k1+2.*k2+2.*k3+k4);
   //    }

   //    return;
   // }

   // void rho_v_p_vdw(const double x0, Vector & xx, Vector & rho, Vector & v, Vector & p)
   // {  
   //    bool is_nL_found = false, is_n0_found = false, is_nR_found = false;
   //    int nL = 0, n0 = 0, nR = 0;
   //    int nmax = xx.Size();

   //    for (int i = 0; i < nmax; i++)
   //    {
   //       if (xx[i] >= xR && !is_nR_found) { 
   //          nR = i - 1;
   //          is_nR_found = true;
   //          assert(i-1 >= 0);
   //          assert(xx[i-1] > xR); 
   //       }
   //       if (xx[i] > xL && !is_nL_found) { 
   //          nL = i;
   //          is_nL_found = true;
   //       } 
   //       if (xx[i] > x0 && !is_n0_found) { 
   //          n0 = i;
   //          is_nL_found = true;
   //       } 
   //    }

   //    rho[0:nL-1] = rhoL;
   //    v[0:nL-1]   = vL;
   //    p[0:nL-1]   = pL;

   //    if (n0-nL>0) {
   //       sol_rho(rho_minus,SL,&rk_rho,xx[n0-1:nL:-1],rho[n0-1:nL:-1]);
   //       DO n = nL, n0-1
   //          v(n) = xx(n) - c(rho(n),SL)
   //          p(n) = SL*(rho(n)/(1-bvdw*rho(n)))**gamma_vdw-avdw*rho(n)**2
   //       END DO
   //    }

   //    rho[nR+1:] = rhoR;
   //    v[nR+1:]   = vR;
   //    p[nR+1:]   = pR;
      
   //    if (nR-n0+1>0) {
   //       sol_rho(rho_plus,SR,&rk_rho,xx[n0:nR],rho[n0:nR]);
   //       for (int i = n0; i < nR; i++)
   //       {
   //          v(i) = xx(i) - c(rho(i),SR);
   //          p(i) = SR * pow(rho[i]/(1. - this->get_b() *rho[i]), this->get_gamma()) - this->get_a() * pow(rho[i],2);
   //       }
   //       DO n = n0, nR
            
   //       END DO
   //    }
   // }

}; // End class

} // ns hydroLO
} // ns mfem

/***
 * Exact Solution process outlined in supplementary.pdf
 * 
 *    1) set rhop
 *    2) set rhom
 *    3) ps
 *    4) es pm
 *    5) cs pm
 *    6) cLR
 *    7) vs
 *    8) sLR
 *    9) rhoz exact - requires functional implementation
*/