#ifndef PROBLEM_BASE
#define PROBLEM_BASE

#include "mfem.hpp"
#include "geometry.hpp"
#include "eos.h"
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

// Fortran subroutine from Lagrangian code (GREEDY)
extern "C" {
   void __arbitrary_eos_lagrangian_greedy_lambda_module_MOD_greedy_lambda_arbitrary_eos(
      double *in_rhol, double *in_ul, double *in_el, double *in_pl,
      double *in_rhor, double *in_ur, double *in_er, double *in_pr,
      double *in_tol, bool *no_iter,double *lambda_max, double *pstar, int *k); ///TODO: , double *b_covolume);
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
namespace hydroLO
{

template<int dim>
class ProblemBase
{
private:
   bool distort_mesh = false;
   bool known_exact_solution = false;
   bool bcs = false; // Indicator for boundary conditions
   // Indicator if mesh velocity boundary conditions need to be updated at each iteration
   // So far the only problem that needs this is Saltzman
   bool mv_bcs_need_updating = false; 
   string indicator = ""; // Possible: saltzmann

   // CFL change
   bool change_cfl = false;
   double cfl_first = 0.5;
   double cfl_second = 0.5;
   double cfl_time_change = 0.;

protected:
   std::unique_ptr<EquationOfState> eos = NULL;  // Pointer to the EOS object
   double a = 0., b = 0., gamma = 0.;

public:
   // Setters
   void set_a(const double &_a) { a = _a; }
   void set_b(const double &_b) { b = _b; }
   void set_gamma(const double &_gamma) { gamma = _gamma; }
   void set_indicator(const string &_ind) { this->indicator = _ind; }
   void set_bcs_indicator(const bool &tvalue) { this->bcs = tvalue; }
   void set_mv_bcs_need_updating_indicator(const bool &tvalue) { this->mv_bcs_need_updating = tvalue; }
   void set_distort_mesh(const bool &_distort_mesh) { distort_mesh = _distort_mesh; }
   void set_exact_solution(const bool &_known_exact_solution) { known_exact_solution = _known_exact_solution; }
   // CFL change
   void set_cfl_change(const bool &_change_cfl) { change_cfl = _change_cfl; }
   void set_cfl_first(const double &_cfl_first) { cfl_first = _cfl_first; }
   void set_cfl_second(const double &_cfl_second) { cfl_second = _cfl_second; }
   void set_cfl_time_change(const double &_cfl_time_change) { cfl_time_change = _cfl_time_change; }

   // Getters
   double get_a() const { return a; }
   double get_b() const { return b; }
   string get_indicator() const { return indicator; }
   bool has_boundary_conditions() const { return bcs; }
   bool get_mv_bcs_need_updating() const { return mv_bcs_need_updating; }
   bool get_distort_mesh() const { return distort_mesh; }
   bool has_exact_solution() const { return known_exact_solution; }
   // CFL change
   bool get_cfl_change() const { return change_cfl; }
   double get_cfl_first() const { return cfl_first; }
   double get_cfl_second() const { return cfl_second; }
   double get_cfl_time_change() const { return cfl_time_change; }

   /* Optionally overridden */
   virtual double get_gamma(const int &cell_attr = 0) const { return gamma; }
   virtual void lm_update(const double b_covolume) {}
   virtual void update(Vector vec, double t = 0.) {}
   virtual void get_additional_BCs(const FiniteElementSpace &fes, Array<int> ess_bdr, Array<int> &add_ess_tdofs, Array<double> &add_bdr_vals, const Geometric<dim> &geom=NULL) { MFEM_ABORT("Function get_additional_BCs must be overridden.\n"); }
   virtual void update_additional_BCs(const double &t, const double timestep_first, Array<double> &add_bdr_vals, const Geometric<dim> &geom=NULL, const ParGridFunction &x_gf=NULL) { MFEM_ABORT("Function get_additional_BCs must be overridden.\n"); }
   virtual void GetBoundaryState(const Vector &x, const int &bdr_attr, Vector &state, const double &t=0.) { MFEM_ABORT("Function GetBoundaryState must be overridden.\n"); }

   /* ProblemDescription */
   static double internal_energy(const Vector &U)
   {
      const double &rho = 1./U[0];
      const double &e = specific_internal_energy(U);
      return rho * e;
   }

   static double specific_internal_energy(const Vector &U)
   {
      Vector v;
      velocity(U, v);
      const double E = U[dim + 1]; // specific total energy
      // return E - 0.5 * pow(v.Norml2(), 2); // Gives issues, see 'quest' branch

      if (E < 1.e-20)
      {
         return 0.;
      }
      else
      {
         double _val = E - 0.5 * pow(v.Norml2(), 2);
         if (_val <= 0.)
         {
            MFEM_ABORT("Negative internal energy computed.\n");
         }
         return _val;
      }

   }

   static inline void velocity(const Vector & U, Vector &vel)
   {
      vel.SetSize(dim);
      for (int i = 0; i < dim; i++)
      {
         vel[i] = U[i+1];
      }
   }

   inline double compute_lambda_max(const Vector & U_i,
                                    const Vector & U_j,
                                    const Vector & n_ij,
                                    double in_pl,
                                    double in_pr,
                                    const bool &use_greedy_viscosity,
                                    double b_covolume=-1.,
                                    const string flag="NA") const
   {
      // cout << "compute_lambda_max\n";
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
         Vector vi, vj;

         in_taul = U_i[0];
         velocity(U_i, vi);
         in_ul = vi * n_ij; 
         in_el = specific_internal_energy(U_i);

         in_taur = U_j[0]; 
         velocity(U_j, vj);
         in_ur = vj * n_ij; 
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

      double lambda_max = 1.;
      if (use_greedy_viscosity)
      {
         want_iter = true; // No iter
         // cout << "inul: " << in_ul << ", inur: " << in_ur << endl;
         __arbitrary_eos_lagrangian_greedy_lambda_module_MOD_greedy_lambda_arbitrary_eos(
            &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,
            &want_iter,&lambda_max, &pstar,&k);
      }
      else 
      {
         __arbitrary_eos_lagrangian_lambda_module_MOD_lambda_arbitrary_eos(
            &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,
            &want_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&k, &_b);

         // bool no_iter = false; 
         // __arbitrary_eos_lambda_module_MOD_lambda_arbitrary_eos(
         //    &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,
         //    &no_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&k, &b_covolume);

         lambda_max = std::max(std::abs(lambda_maxl_out), std::abs(lambda_maxr_out));
      }

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

      // cout << "lambda_max: " << lambda_max << endl;
      return lambda_max;
      // return 0.5;
   }

   inline DenseMatrix flux(const Vector &U, const int &cell_attr=0) const
   {
      DenseMatrix result(dim+2, dim);

      Vector v; velocity(U, v);
      const double rho = 1. / U[0];
      const double sie = specific_internal_energy(U);
      const double p = pressure(rho, sie, cell_attr);

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

   inline double sound_speed(const double &rho, const double &press, const int &cell_attr=0) const
   {
      double val = this->get_gamma(cell_attr) * (press + this->get_a() * pow(rho,2)) / (rho * (1. - this->get_b() * rho));
      val -= 2. * this->get_a() * rho;
      val = pow(val, 0.5);
      return val;
   }

   /*********************************************
    * Functions describing the initial state
    ********************************************/
   double sv0(const Vector &x, const double & t) const
   {
      double val = this->rho0(x,t);
      assert(val != 0.);
      return 1./val;
   }

   double ste0(const Vector &x, const double & t) const
   {
      Vector v(dim);
      this->v0(x,t,v);
      return this->sie0(x, t) + 0.5 * pow(v.Norml2(), 2);
   }

   /* To compute pressure, utilize our equation of state */
   double pressure(const double &rho, const double &sie, const double &gamma) const {
      assert(eos != NULL);
      return this->eos->pressure(rho, sie, gamma);
   }

   double pressure(const double &rho, const double &sie, const int &cell_attr=0) const {
      return pressure(rho, sie, this->get_gamma(cell_attr));
   }

   /*********************************************
    * Functions to be overridden
    ********************************************/
   virtual double rho0(const Vector &x, const double & t) const {
      MFEM_ABORT("Must override rho0 in ProblemBase class.\n");
   } // virtual function, must be overridden
   virtual void v0(const Vector &x, const double & t, Vector &v) const {
      MFEM_ABORT("Must override v0 in ProblemBase class.\n");
   } // virtual function, must be overridden
   virtual double sie0(const Vector &x, const double & t) const {
      MFEM_ABORT("Must override sie0 in ProblemBase class.\n");
   } // virtual function, must be overridden
   virtual double p0(const Vector &x, const double & t) const {
      MFEM_ABORT("Must override p0 in the ProblemBase class.\n");
   } // virtual function, must be overridden
   virtual double gamma_func(const Vector &x = Vector(), const double &t = 0) const {
      return gamma;
   } // optionally overridden
   
}; // End ProblemBase


} // ns hydroLO
} // ns mfem

#endif // PROBLEM_BASE