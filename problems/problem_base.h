#ifndef PROBLEM_BASE
#define PROBLEM_BASE

#include "mfem.hpp"
#include <cmath>
#include <string>

using namespace mfem;
using namespace std;

// Fortran subroutine from Lagrangian code
extern "C" {
   void __arbitrary_eos_lagrangian_lambda_module_MOD_lambda_arbitrary_eos(
      double *in_rhol, double *in_ul, double *in_el, double *in_pl,
      double *in_rhor, double *in_ur, double *in_er, double *in_pr,
      double *in_tol, bool *want_iter,double *lambda_maxl_out,
      double *lambda_maxr_out, double *pstar, int *k, double *b_covolume);
}

// Fortran subroutine from Eulerian code
extern "C" {
   void __arbitrary_eos_lambda_module_MOD_lambda_arbitrary_eos(
      double *in_rhol, double *in_ul, double *in_el, const double *in_pl,
      double *in_rhor, double *in_ur, double *in_er, const double *in_pr,
      double *in_tol, bool *no_iter,double *lambda_maxl_out,
      double *lambda_maxr_out, double *pstar, int *k, const double *b_covolume);
}

namespace mfem
{
namespace hydrodynamics
{

template<int dim>
class ProblemBase
{
private:
   double a = 0., b = 0., gamma = 0.;
   bool distort_mesh = false;
   bool known_exact_solution = false;
   bool bcs = false; // Indicator for boundary conditions
   string indicator = ""; // Possible: saltzmann

   // CFL change
   bool _change_cfl = false;
   constexpr static double CFL_first = 0.5;
   constexpr static double CFL_second = 0.5;
   constexpr static double CFL_time_change = 0.01; // From Boscheri's paper

public:
   virtual bool change_cfl() { return _change_cfl; }
   virtual double get_cfl_first() { return CFL_first; }
   virtual double get_cfl_second() { return CFL_second; }
   virtual double get_cfl_time_change() { return CFL_time_change; }
   
   // Setters
   void set_a(const double &_a) { a = _a; }
   void set_b(const double &_b) { b = _b; }
   void set_gamma(const double &_gamma) { gamma = _gamma; }
   void set_indicator(const string &_ind) { this->indicator = _ind; }


   // Getters
   double get_a() { return a; }
   double get_b() { return b; }
   string get_indicator() { return indicator; }
   virtual double get_gamma(const int &cell_attr = 0) { return gamma; }
   virtual bool get_distort_mesh() { return distort_mesh; }
   virtual bool has_exact_solution() { return known_exact_solution; }
   virtual bool has_boundary_conditions() { return bcs; }
   
   /* Functions that update the class, can be overridden */
   virtual void lm_update(const double b_covolume) {
   }
   virtual void update(Vector x_gf, double t = 0.) {
   }

   /* ProblemDescription */
   static double internal_energy(const Vector &U)
   {
      const double &rho = 1./U[0];
      const double &e = specific_internal_energy(U);
      return rho * e;
   }

   static double specific_internal_energy(const Vector &U)
   {
      const Vector v = velocity(U);
      const double E = U[dim + 1]; // specific total energy
      return E - 0.5 * pow(v.Norml2(), 2);
   }

   static inline Vector velocity(const Vector & U)
   {
      Vector v;
      v.SetSize(dim);
      Array<int> dofs;
      for (int i = 0; i < dim; i++)
      {
         v[i] = U[i+1];
      }

      return v;
   }

   inline double compute_lambda_max(const Vector & U_i,
                                    const Vector & U_j,
                                    const Vector & n_ij,
                                    double in_pl,
                                    double in_pr,
                                    double b_covolume=-1.,
                                    const string flag="NA")
   {
      double in_taul, in_ul, in_el, in_taur, in_ur, in_er, in_rhol, in_rhor;
      if (flag == "testing")
      {
         in_taul = U_i[0];
         in_ul = U_i[1];
         in_el = U_i[dim+1];

         in_taur = U_j[0]; 
         in_ur = U_j[1];
         in_er = U_j[dim+1];
      }
      else 
      {
         assert(flag == "NA");

         in_taul = U_i[0];
         in_ul = velocity(U_i) * n_ij; 
         in_el = specific_internal_energy(U_i);

         in_taur = U_j[0]; 
         in_ur = velocity(U_j) * n_ij; 
         in_er = specific_internal_energy(U_j);
      }

      in_rhol = 1. / in_taul;
      in_rhor = 1. / in_taur;

      double in_tol = 1.e-16,
             _b=0.;
      double lambda_maxl_out = 0.,
             lambda_maxr_out = 0.,
             pstar = 0.,
             vstar = 0.;
      int k = 0; // Tells you how many iterations were needed for convergence
      bool want_iter = false;

      // Handle b_covolume parameter
      if (b_covolume == -1.)
      {
         // if b_covolume is not specified, or is -1., use the problem specific value.
         _b = this->get_b();
      }
      else
      {
         // if b_covolume is specified, use this value
         _b = b_covolume;
      }

      __arbitrary_eos_lagrangian_lambda_module_MOD_lambda_arbitrary_eos(
         &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,
         &want_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&k, &_b);

      // bool no_iter = false; 
      // __arbitrary_eos_lambda_module_MOD_lambda_arbitrary_eos(
      //    &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,
      //    &no_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&k, &b_covolume);

      double lambda_max = std::max(std::abs(lambda_maxl_out), std::abs(lambda_maxr_out));

      if (isnan(lambda_max))
      {
         cout << "nij:\n";
         n_ij.Print(cout);
         cout << "Ui:\n";
         U_i.Print(cout);
         cout << "Uj:\n";
         U_j.Print(cout);
         cout << "in_taul: " << in_taul << ", ul: " << in_ul << ", el: " << in_el << ", pl: " << in_pl << endl;
         cout << "in_taur: " << in_taur << ", ur: " << in_ur << ", er: " << in_er << ", pr: " << in_pr << endl;
         MFEM_ABORT("NaN values returned by lambda max computation!\n");
      }

      // cout << "nij:\n";
      // n_ij.Print(cout);
      // cout << "b: " << b_covolume << endl;
      // cout << "UL. Density: " << 1./U_i[0] << ", vel: " << U_i[1] << ", ste: " << U_i[dim+1] << ", p: " << in_pl << endl;
      // cout << "UR. Density: " << 1./U_j[0] << ", vel: " << U_j[1] << ", ste: " << U_j[dim+1] << ", p: " << in_pr << endl;
      // cout << "lamba L: " << std::abs(lambda_maxl_out) << ", lambda_R: " <<  std::abs(lambda_maxr_out) << endl;

      return lambda_max;
      // return 0.5;
   }

   inline DenseMatrix flux(const Vector &U)
   {
      DenseMatrix result(dim+2, dim);

      const Vector v = velocity(U);
      const double p = pressure(U);

      // * is not overridden for Vector class, but *= is
      Vector v_neg = v, vp = v;
      v_neg *= -1.;
      vp *= p;

      // Set f(U) according to (2.1c)
      result.SetRow(0,v_neg);

      for (int i = 0; i < dim; i++)
      {
         result(i+1, i) = p;
      }

      result.SetRow(dim+1, vp);

      return result;
   }

   inline double sound_speed(const Vector &U)
   {
      double _pressure = this->pressure(U);
      double density = 1. / U[0];

      double val = this->get_gamma() * (_pressure + this->get_a() * pow(density,2)) / (density * (1. - this->get_b() * density));
      val -= 2. * this->get_a() * density;
      val = pow(val, 0.5);
      return val * density;
   }

   /*********************************************
    * Functions describing the initial state
    ********************************************/
   double sv0(const Vector &x, const double & t)
   {
      double val = this->rho0(x,t);
      assert(val != 0.);
      return 1./val;
   }

   double ste0(const Vector &x, const double & t)
   {
      Vector v(dim);
      this->v0(x,t,v);
      return this->sie0(x, t) + 0.5 * pow(v.Norml2(), 2);
   }

   /*********************************************
    * Functions to be overridden
    ********************************************/
   virtual double pressure(const Vector &U, const int &cell_attr=0) {
      MFEM_ABORT("Must override pressure in ProblemBase class.\n");
      return 1.;
   } // virtual function, must be overridden
   virtual double rho0(const Vector &x, const double & t) {
      MFEM_ABORT("Must override rho0 in ProblemBase class.\n");
      return 1.;
   } // virtual function, must be overridden
   virtual void v0(const Vector &x, const double & t, Vector &v) {
      MFEM_ABORT("Must override v0 in ProblemBase class.\n");
      return;
   } // virtual function, must be overridden
   virtual double sie0(const Vector &x, const double & t) {
      MFEM_ABORT("Must override sie0 in ProblemBase class.\n");
      return 1.;
   } // virtual function, must be overridden
   
}; // End ProblemBase


} // ns hydrodynamics
} // ns mfem

#endif // PROBLEM_BASE