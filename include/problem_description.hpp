#include "mfem.hpp"
#include <cmath>

using namespace mfem;

// Fortran subroutine
extern "C" {
   void __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      double *in_rhol, double *in_ul, double *in_el, double *in_pl,
      double *in_rhor, double *in_ur, double *in_er, double *in_pr,
      double *in_tol, bool *no_iter,double *lambda_maxl_out,
      double *lambda_maxr_out, double *pstar, int *k);
}

template <int dim>
class ProblemDescription
{
public:
   static constexpr unsigned int problem_dimension = 2 + dim;

   using state_type = Vector;
   using flux_type = DenseMatrix;

   static constexpr double gamma = 7. / 5.;

   static inline double internal_energy(const state_type &U);

   static inline double pressure(const state_type &U);

   static inline double compute_lambda_max(const state_type & U_i,
                                           const state_type & U_j,
                                           const Vector & n_ij);
   
   static inline Vector velocity(const state_type & U);

   static inline flux_type flux(const state_type &U);
};

/* Function definitions */
template <int dim>
inline Vector
ProblemDescription<dim>::velocity(const state_type & U)
{
   Vector v;
   const double &rho = 1./U[0];
   v.SetSize(dim);
   Array<int> dofs;
   for (int i = 0; i < dim; i++)
   {
      dofs.Append(i + 1);
   }
   U.GetSubVector(dofs, v);

   return v;
}

template <int dim>
inline double
ProblemDescription<dim>::internal_energy(const Vector & U)
{
   const double &rho = 1./U[0];
   const Vector  v   = velocity(U);
   const double &E   = U[dim + 1];
   return rho * (E - 0.5 * pow(v.Norml2(), 2));
}

template <int dim>
inline double
ProblemDescription<dim>::pressure(const Vector & U)
{
   return (gamma - 1.) * internal_energy(U);
}

template<int dim>
inline double 
ProblemDescription<dim>::compute_lambda_max(const Vector & U_i,
                                            const Vector & U_j,
                                            const Vector & n_ij)
{
   std::cout << "Computing Viscosity\n";
   double in_rhol = 1. / U_i[0];
   double in_ul = U_i * n_ij; // TODO: Validate.
   double in_el = internal_energy(U_i);
   double in_pl = pressure(U_i);

   double in_rhor = 1. / U_j[0]; 
   double in_ur = U_j * n_ij; // TODO: Check that this shouldn't be n_ji.
   double in_er = internal_energy(U_j);
   double in_pr = pressure(U_j);

   double in_tol,lambda_maxl_out,lambda_maxr_out,pstar;
   bool no_iter = true;
   int k = 10;
   __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,
      &in_tol,&no_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&k);

   return 1.0;
}

/* Explicit Instantiation */
// template class ProblemDescription<1>;
// template class ProblemDescription<2>;
// template class ProblemDescription<3>;