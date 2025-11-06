#ifndef PROBLEM_BASE
#define PROBLEM_BASE

#include "mfem.hpp"
#include "geometry.hpp"
#include "eos.h"
#include "elastic.hpp"
#include "shear_closure.hpp"
#include "laglos_tools.hpp"
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
      double *lambda_maxr_out, double *pstar, int *k, 
      double *b_covolume, double *p_inf);
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

enum ShearEOS {
   NEO_HOOKEAN,
   MOONEY_RIVLIN,
   AORTIC,         // anisotropic model
   TRANSVERSELY_ISOTROPIC
};

class ProblemBase
{
private:
   bool distort_mesh = false;
   bool known_exact_solution = false;
   bool th_bcs = false; // Indicator for thermo boundary conditions to be imposed on the thermo solution
   bool mv_bcs = false; // Indicator for mesh velocity boundary conditions to be imposed on the mesh velocity solution  
   // Indicator if mesh velocity boundary conditions need to be updated at each iteration
   // So far the only problem that needs this is Saltzman
   bool mv_bcs_need_updating = false; 
   string indicator = ""; // Possible: saltzmann

   // CFL change
   bool change_cfl = false;
   double cfl_first = 0.5;
   double cfl_second = 0.5;
   double cfl_time_change = 0.;
   int elastic_eos_type = -1; // Indicator for type of elastic eos model
   ShearEOS shear_eos;
   Elastic *elastic = nullptr;
   Vector rho0_v;
   /* 
   If using the aortic eos, a fiber direction must be defined. 
   This direction is used to compute the invariants:
      m = (cos(theta), sin(theta), 0) 
   */
   double mu = -1.;
   double theta = M_PI/4.; // angle of fiber direction
   Vector mi_vec; // Fiber direction vector for anisotropic models
   DenseMatrix Gi;
   ShearClosure *shear_closure_model = nullptr;

   // double A1 = 21.5802, B1 = 9.9007, D1 = 0.8849, w1 = 0.4189; /* PP too squished */
   // double A1 = 771.8033, B1 = 21.2093, D1 = 3.8068, w1 = 0.4971;
   // double A2 = 1., B2 = 1.;
   /* These params were closer to the Neo Hook results when w1 = 0 */
   double stiffness = 9.63E5;
   double A1 = 0.5 * stiffness, B1 = 0.5 * stiffness, A2 = 1., B2 = 1., D1 = 0.5 * (1.5*stiffness), w1 = 0.49;

protected:
   const int dim;
   double a = 0., b = 0., gamma = 0.;
   double q = 0., p_inf = 0.;
   std::unique_ptr<EquationOfState> eos = NULL;  // Pointer to the EOS object

   /* helpers */
   static std::string toFixed(double x, int prec = 6) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(prec) << x;
        return ss.str();
    }

    static std::string shearEOSToString(ShearEOS e) {
        switch (e) {
            case ShearEOS::NEO_HOOKEAN:            return "NEO_HOOKEAN";
            case ShearEOS::MOONEY_RIVLIN:          return "MOONEY_RIVLIN";
            case ShearEOS::AORTIC:                 return "AORTIC";
            case ShearEOS::TRANSVERSELY_ISOTROPIC: return "TRANSVERSELY_ISOTROPIC";
        }
        return "Unknown";
    }

    static std::string vecSizeString(const Vector& v) {
      return std::to_string(v.Size());
    }
    static std::string vecHead(const Vector& v, int k = 3) {
      // implement: "[v0, v1, v2]" or "(empty)"1
      if (k == 0) return "(empty)";
      else {
         for (int i = 0; i < std::min(k, v.Size()); i++) {
            std::cout << v[i];
            if (i < std::min(k, v.Size()) - 1) std::cout << ", ";
         }
         if (v.Size() > k) std::cout << ", ...";
         return "";
      }
    }

public:
   ProblemBase(const int &_dim) : dim(_dim) {}
   // Setters
   void set_a(const double &_a) { a = _a; }
   void set_b(const double &_b) { b = _b; }
   void set_gamma(const double &_gamma) { gamma = _gamma; }
   void set_pinf(const double &_p_inf) { p_inf = _p_inf; }
   void set_shear_modulus(const double &_mu) { mu = _mu; }
   void set_indicator(const string &_ind) { this->indicator = _ind; }
   void set_thbcs_indicator(const bool &tvalue) { this->th_bcs = tvalue; }
   void set_mvbcs_indicator(const bool &tvalue) { this->mv_bcs = tvalue; }
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
   bool has_th_boundary_conditions() const { return th_bcs; }
   bool has_mv_boundary_conditions() const { return mv_bcs; }
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
   virtual double get_pinf(const int &cell_attr = 0) const { return p_inf; }
   virtual double get_shear_modulus(const int &cell_attr = 0) const { return mu; }
   virtual void lm_update(const double b_covolume) {}
   virtual void update(Vector vec, double t = 0.) {}
   virtual void get_additional_BCs(const FiniteElementSpace &fes, Array<int> ess_bdr, Array<int> &add_ess_tdofs, Array<double> &add_bdr_vals, const Geometric *geom=NULL) { MFEM_ABORT("Function get_additional_BCs must be overridden.\n"); }
   virtual void update_additional_BCs(const double &t, const double timestep_first, Array<double> &add_bdr_vals, const Geometric *geom=NULL, const ParGridFunction *x_gf=NULL) { MFEM_ABORT("Function get_additional_BCs must be overridden.\n"); }
   virtual void GetBoundaryState(const Vector &x, const int &bdr_attr, Vector &state, const double &t=0.) { MFEM_ABORT("Function GetBoundaryState must be overridden.\n"); }

   void PrintInfo() const
   {
      const int W = 24; // label width for alignment
      auto line = [&] { std::cout << "@------------------------------------------@\n"; };
      auto head = [&](const std::string& title) {
         std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
         std::cout << "@ " << std::setw(41) << std::left << title << "@\n";
         line();
      };
      auto field = [&](const std::string& key, const std::string& val) {
         std::cout << "@ " << std::setw(W) << std::left << (key + " : ") 
                     << std::setw(17) << std::left << val << "@\n";
      };
      auto b2s = [](bool b){ return b ? "true" : "false"; };

      head("ProblemBase Configuration");

      // Booleans
      field("distort_mesh",         b2s(distort_mesh));
      field("known_exact_solution", b2s(known_exact_solution));
      field("th_bcs",               b2s(th_bcs));
      field("mv_bcs",               b2s(mv_bcs));
      field("mv_bcs_need_updating", b2s(mv_bcs_need_updating));

      // Indicator string
      field("indicator", indicator.empty() ? "(none)" : indicator);

      line();

      // CFL
      field("change_cfl",       b2s(change_cfl));
      field("cfl_first",        toFixed(cfl_first,2));
      field("cfl_second",       toFixed(cfl_second,2));
      field("cfl_time_change",  toFixed(cfl_time_change,2));

      line();

      // Pointers (presence + optional nested printing)
      field("elastic", elastic ? "set" : "null");
      if (elastic) {
         line();
         // Nested elastic configuration output
         field("shear_eos", shearEOSToString(shear_eos));
         field("shear_closure_model", shear_closure_model ? "set" : "null");

         line();

         // Material / anisotropy
         field("mu",          toFixed(mu,1));
         field("theta (rad)", toFixed(theta));
         field("theta (deg)", toFixed(theta * 180.0 / M_PI));

         // Vectors / matrices: sizes + preview
         field("rho0_v size", vecSizeString(rho0_v));
         field("mi_vec size", vecSizeString(mi_vec));
         field("mi_vec head", vecHead(mi_vec));
         field("Gi shape",    "3x3");
         line();
      }

      std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
   }
   /* ProblemDescription */
   double internal_energy(const Vector &U, const double e_shear)
   {
      const double &rho = 1./U[0];
      const double &e = specific_internal_energy(U, e_shear);
      return rho * e;
   }

   double specific_internal_energy(const Vector &U, const double e_shear)
   {
      // Verify e_shear > 0.
      // TODO: This is a temporary fix. The correct fix is to ensure that e_shear is always positive.
      // assert(e_shear >=  -std::numeric_limits<double>::epsilon());
      /* Subtract out kinetic and shear energy */
      Vector v;
      velocity(U, v);
      double E = U[dim + 1]; // specific total energy
      double _ske = 0.5 * pow(v.Norml2(), 2); // specific kinetic energy
      if (std::isnan(e_shear) || std::isnan(E) || std::isnan(_ske))
      {
         cout << "E: " << E << ", v: " << v.Norml2() << ", e_shear: " << e_shear << ", _ske: " << _ske << endl;
         MFEM_ABORT("NaN values in specific internal energy computation.\n");
      }
      double val = E - _ske - e_shear;

      // Verify sie > 0.
      if (val <  -std::numeric_limits<double>::epsilon())
      {
         cout << "sie computation. E: " << E << ", ke: " << 0.5 * pow(v.Norml2(), 2) << ", esheer: " << e_shear << ", sie: " << val << endl;
         MFEM_ABORT("Specific internal energy is negative or zero.\n");
      }

      return val;
   }

   inline void velocity(const Vector & U, Vector &vel)
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
                                    double esl, // e_sheer_left
                                    double esr, // e_sheer_right
                                    double in_pl,
                                    double in_pr,
                                    const bool &use_greedy_viscosity,
                                    const string flag="NA")
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
         in_el = specific_internal_energy(U_i, esl);

         in_taur = U_j[0]; 
         velocity(U_j, vj);
         in_ur = vj * n_ij; 
         in_er = specific_internal_energy(U_j, esr);
      }

      in_rhol = 1. / in_taul;
      in_rhor = 1. / in_taur;

      double in_tol = 1.e-16;
      double lambda_maxl_out = 0.,
             lambda_maxr_out = 0.,
             pstar = 0.,
             vstar = 0.;
      int k = 0; // Tells you how many iterations were needed for convergence
      bool want_iter = false;

      double _b = this->get_b();
      double _p_inf = this->get_pinf();

      double lambda_max = 1.;
      if (use_greedy_viscosity)
      {
         want_iter = true; // No iter
         // cout << "inul: " << in_ul << ", inur: " << in_ur << endl;
         // MFEM_ABORT("Need to incorporate b_covolume and p_inf");
         __arbitrary_eos_lagrangian_greedy_lambda_module_MOD_greedy_lambda_arbitrary_eos(
            &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,
            &want_iter,&lambda_max, &pstar,&k);
      }
      else 
      {
         __arbitrary_eos_lagrangian_lambda_module_MOD_lambda_arbitrary_eos(
            &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,
            &want_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&k, &_b, &_p_inf);

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

   /**
    * @brief Computes the flux matrix for a given state vector U.
    *
    * This function calculates the flux matrix based on the input state vector U
    * and an optional cell attribute. The flux matrix is used in the context of
    * solving partial differential equations, particularly in fluid dynamics.
    *
    * @param U The state vector for which the flux is to be computed.
    * @param cell_attr An optional integer representing cell attributes, default is 0.
    * @return A DenseMatrix representing the flux of the state vector U.
    *
    * Note: This function assumes elasticity is not being used.
    */
   inline DenseMatrix flux(const Vector &U, const int &cell_attr=0)
   {
      if (cell_attr == 50)
      {
         MFEM_ABORT("Cell attr 50 is reserved for elastic.\n");
      }
      DenseMatrix result(dim+2, dim);

      Vector v; velocity(U, v);
      const double rho = 1. / U[0];
      const double sie = specific_internal_energy(U, 0.);
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

   inline DenseMatrix flux_ex_p(const Vector &U, const Vector &x, const double &t, const int &cell_attr=0)
   {
      MFEM_ABORT("This function is deprecated. Use flux instead.\n");
      if (cell_attr == 50)
      {
         MFEM_ABORT("Cell attr 50 is reserved for elastic.\n");
      }
      DenseMatrix result(dim+2, dim);

      Vector v(dim); 
      // velocity(U, v);
      this->v0(x, t, v); // Use initial velocity function
      // const double sie = specific_internal_energy(U, 0.);
      const double p = p0(x,t);

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

   //NF//MS add a compute sigma routine here you have to compute the F tensor in each points and its average
   // on the cell. then compute C then sigma, you will need the pressure in that subroutine
   // this as to be compute 
   /**
    * @brief Computes the stress tensor sigma by adding the deviatoric stress tensor and the pressure term.
    *
    * This function calculates the stress tensor `sigma` by first extracting the dim x dim deviatoric stress tensor 
    * from `sig_dev` and then adding the pressure term to it. The pressure is computed based on the 
    * state vector `U` and an optional cell attribute `cell_attr`.
    *
    * @param[in] e The element index.
    * @param[in] pressure The pressure value.
    * @param[out] sigma The resulting stress tensor.
    * @param[in] cell_attr An optional cell attribute used in the pressure calculation (default is 0).
    *
    * NOTE: @param sig_dev is not assumed to be dim x dim.
    */
   inline void ComputeSigma(const int &e, const double &pressure, DenseMatrix &sigma, const int&cell_attr=0)
   {
      // cout << "ProblemBase::ComputeSigma\n";
      DenseMatrix sig_dev(3);
      ComputeS(e, sig_dev);

      DenseMatrix I(dim);
      I = 0.;
      Array<int> idx(dim);
      for (int i = 0; i < dim; i++) {
         idx[i] = i;
         I(i,i) = 1.;
      }

      sig_dev.GetSubMatrix(idx, idx, sigma);
      assert(sigma.NumRows() == dim && sigma.NumCols() == dim);
      mfem::Add(sigma, I, -1. * pressure, sigma);
   }

   /**
    * @brief Computes the elastic flux for a given deviatoric stress tensor and displacement vector.
    *
    * This function calculates the elastic flux based on the provided deviatoric stress tensor (`sig_dev`),
    * state vector (`U`), and an optional cell attribute (`cell_attr`). It first computes a flux matrix
    * using the state vector and cell attribute, then computes the stress tensor `sigma` and sets it
    * in the appropriate submatrix of the result. The function currently contains an assertion that always fails.
    *
    * @param sig_dev The deviatoric stress tensor as a DenseMatrix.
    * @param U The state vector.
    * @param cell_attr An optional integer representing the cell attribute (default is 0).
    * @return A DenseMatrix representing the computed elastic flux.
    *
    * NOTE: @param sig_dev is not assumed to be dim x dim.
    */
   inline DenseMatrix ElasticFlux(const int &e, const Vector &U, const int &cell_attr=0)
   {
      if (cell_attr != 50)
      {
         MFEM_WARNING("Elastic cell must have attribute 50.\n");
         return flux(U, cell_attr);
      }
      // cout << "ProblemBase::ElasticFlux\n";
      DenseMatrix result(dim+2, dim);
      // Indices from stress tensor in result
      Array<int> idx(dim), idy(dim);
      for (int i = 0; i < dim; i++) { 
         idx[i] = 1 + i;
         idy[i] = i;
      }

      /* Compute dim x dim stress tensor */
      DenseMatrix sigma(dim);
      const double rho = 1./ U[0];
      const double es = e_shear(e);
      const double sie = specific_internal_energy(U, es);
      const double _pressure = pressure(rho, sie, cell_attr);

      ComputeSigma(e, _pressure, sigma, cell_attr);

      // * is not overridden for Vector class, but *= is
      Vector v; velocity(U, v);
      Vector v_neg = v, sigmav(dim);
      v_neg *= -1.;
      sigma *= -1.;
      sigma.Mult(v, sigmav);

      // Set entries in flux
      result.SetRow(0,v_neg);
      result.SetSubMatrix(idx, idy, sigma);
      result.SetRow(dim+1, sigmav);

      return result;
   }

   void SetElastic(Elastic &_elastic) { elastic = &_elastic; }

   bool HasElastic() const { return (elastic != nullptr); }

   void SetElasticEOSType(const int &_elastic_eos_type) { 
      elastic_eos_type = _elastic_eos_type; 
      switch (elastic_eos_type)
      {
         case 1: // Neo-Hookean
            shear_eos = ShearEOS::NEO_HOOKEAN;
            shear_closure_model = new ShearClosureNeoHookean(mu);
            break;
         case 2: // Mooney-Rivlin
            shear_eos = ShearEOS::MOONEY_RIVLIN;
            shear_closure_model = new ShearClosureMooneyRivlin(mu);
            break;
         case 3: // Aortic
         {
            /* Invariants */
            mi_vec.SetSize(3);
            Gi.SetSize(3);
            /* Set fiber direction */
            mi_vec(0) = cos(theta);
            mi_vec(1) = sin(theta);
            mi_vec(2) = 0.;
            tensor(mi_vec, mi_vec, Gi);
            shear_eos = ShearEOS::AORTIC;
            shear_closure_model = new ShearClosureAortic(mu, mi_vec, w1, D1, A1, B1);
            break;
         }
         case 4: // Transversely Isotropic
         {
            shear_eos = ShearEOS::TRANSVERSELY_ISOTROPIC;
            double E = 1.729E6, EA = 1.E6, GA = 4.59E5, nu = 0.4;
            mi_vec.SetSize(3);
            /* Set fiber direction */
            mi_vec(0) = cos(theta);
            mi_vec(1) = sin(theta);
            mi_vec(2) = 0.;
            shear_closure_model = new ShearClosureTransverselyIsotropic(mu, mi_vec, E, EA, GA, nu);
            
            break;
         }
         default:
            MFEM_ABORT("Invalid value for shear_eos.");
      }
      
   }

   void set_rho0_vec(const Vector &_rho0_vec) { rho0_v = _rho0_vec; }

   virtual double e_shear(const int &e) const
   {
      const double rho0 = rho0_v(e);

      DenseMatrix F(dim);
      elastic->ComputeAvgF(e, F);
      double _e_shear = shear_closure_model->ComputeShearEnergy(F, rho0);
      return _e_shear;
   }

   void ComputeS(const int &e, DenseMatrix &S) const
   {
      /* Objects needed in function */
      DenseMatrix F(3);

      /* Compute F */
      elastic->ComputeAvgF(e, F);

      shear_closure_model->ComputeCauchyStress(F, S);
      return;
   }

   void ComputeF(const int &e, DenseMatrix &F) const
   {
      /* Compute F */
      elastic->ComputeAvgF(e, F);
      return;
   }

   double sound_speed(const double &rho, const double &press, const double &gamma) const
   {
      assert(eos != NULL);
      return this->eos->sound_speed(rho, press, gamma);
   }

   double sound_speed(const double &rho, const double &press, const int &cell_attr=0) const {
      return sound_speed(rho, press, this->get_gamma(cell_attr));
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

   /**
    * @brief Computes the total specific energy at a given point and time.
    *
    * This function calculates the total specific energy by adding the specific internal energy
    * and the contribution of the kinetic energy at a given point and time.
    *
    * @param x The spatial coordinates as a Vector.
    * @param t The time as a double.
    * @return The total specific energy as a double.
    *
    * NOTE: There is no shear energy at initial time.
    */
   double ste0(const Vector &x, const double & t)
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
      const double _g = this->get_gamma(cell_attr);
      return pressure(rho, sie, _g);
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

   /* Vectorized density and energy functions for extrapolated test cases */
   void rho0_vec(const Vector &x, const double & t, Vector &rho_v) const {
      rho_v.SetSize(1);
      rho_v[0] = this->rho0(x,t);
   }
   void sie0_vec(const Vector &x, const double & t, Vector &sie_v) const {
      sie_v.SetSize(1);
      sie_v[0] = this->sie0(x,t);
   }
   
}; // End ProblemBase


} // ns hydroLO
} // ns mfem

#endif // PROBLEM_BASE