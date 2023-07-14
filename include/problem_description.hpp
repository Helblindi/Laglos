#include "mfem.hpp"
#include <cmath>

using namespace mfem;
using namespace std;

// Fortran subroutine
extern "C" {
   void __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      double *in_rhol, double *in_ul, double *in_el, double *in_pl,
      double *in_rhor, double *in_ur, double *in_er, double *in_pr,
      double *in_tol, bool *no_iter,double *lambda_maxl_out,
      double *lambda_maxr_out, double *pstar, double *vstar, int *k, double *b_covolume);
}

namespace mfem
{
namespace hydrodynamics
{

template <int dim, int problem>
class ProblemDescription
{
public:
   static constexpr unsigned int problem_dimension = 2 + dim;

   static constexpr double gamma = 7. / 5.;

   static inline double internal_energy(const Vector &U);

   static inline double specific_internal_energy(const Vector &U);

   static inline double pressure(const Vector &U);

   static inline double gamma_func(const int shocktube = 0);

   static inline double compute_lambda_max(const Vector & U_i,
                                           const Vector & U_j,
                                           const Vector & n_ij,
                                           const string flag="NA");
   
   static inline Vector velocity(const Vector & U);

   static inline DenseMatrix flux(const Vector &U);

   static inline void get_shocktube_vals(const int shocktube, 
                                         double & _rhoL, 
                                         double & _rhoR, 
                                         double & _rhoLstar, 
                                         double & _rhoRstar, 
                                         double & _vL, 
                                         double & _vR, 
                                         double & _vstar, 
                                         double & _pL, 
                                         double & _pR, 
                                         double & _pstar, 
                                         double & _lambda1m, 
                                         double & _lambda1p, 
                                         double & _lambda3);
};

/* Function definitions */
template <int dim, int problem>
inline Vector
ProblemDescription<dim, problem>::velocity(const Vector & U)
{
   Vector v;
   v.SetSize(dim);
   Array<int> dofs;
   for (int i = 0; i < dim; i++)
   {
      dofs.Append(i + 1);
   }
   U.GetSubVector(dofs, v);

   return v;
}

template<int dim, int problem> inline
double ProblemDescription<dim, problem>::internal_energy(const Vector & U)
{
   const double &rho = 1./U[0];
   const double &e = specific_internal_energy(U);
   return rho * e;
}

template<int dim, int problem> inline
double ProblemDescription<dim, problem>::specific_internal_energy(const Vector & U)
{
   const Vector  v   = velocity(U);
   const double &E   = U[dim + 1]; // specific total energy
   return E - 0.5 * pow(v.Norml2(), 2);
}

template<int dim, int problem> inline
double ProblemDescription<dim, problem>::pressure(const Vector & U)
{
   switch (problem)
   {
      case 0:
      {
         assert(dim==1);
         return 1.;
      }
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      case 7:
      default:
      {
         double gamma = gamma_func();
         return (gamma - 1.) * internal_energy(U);
      }
   }
}

template<int dim, int problem> inline
DenseMatrix ProblemDescription<dim, problem>::flux(const Vector &U)
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

template<int dim, int problem> inline
double ProblemDescription<dim, problem>::gamma_func(const int shocktube)
{
   switch (problem)
   {
      case 0:
      case 1:
      case 2:
      case 3:
      case 5:
      case 8:
      {
         return 7./5.;
      }
      case 6:
      {
         switch (shocktube)
         {
            case 0:
            case 1:                    // Sod
            case 2: { return 7./5.; }  // Lax
            case 3: { return 5./3.; }  // Leblanc
            default:
            {
               MFEM_ABORT("Bad number given for shocktube case number.");
               return 0.0;
            }
         }
      }
      case 7:
      case 4:
      {
         return 5./3.;
      }
      default:
      {
         MFEM_ABORT("Bad number given for problem id!"); 
         return 0.0;
      }
   }
}

template<int dim, int problem> inline 
double ProblemDescription<dim, problem>::compute_lambda_max(
   const Vector & U_i,
   const Vector & U_j,
   const Vector & n_ij,
   const string flag) // default is 'NA'
{
   double in_rhol, in_ul, in_el, in_pl, in_rhor, in_ur, in_er, in_pr;
   if (flag == "testing")
   {
      in_rhol = U_i[0];
      in_ul = U_i[1];
      in_pl = U_i[2];
      in_el = U_i[3];

      in_rhor = U_j[0]; 
      in_ur = U_j[1];
      in_pr = U_j[2];
      in_er = U_j[3];
   }
   else 
   {
      assert(flag == "NA");

      in_rhol = 1. / U_i[0];
      in_ul = velocity(U_i) * n_ij; 
      in_el = specific_internal_energy(U_i);
      in_pl = pressure(U_i);

      in_rhor = 1. / U_j[0]; 
      in_ur = velocity(U_j) * n_ij; 
      in_er = specific_internal_energy(U_j);
      in_pr = pressure(U_j);
   }

   // double b_covolume = 0.1/max(in_rhol,in_rhor);
   double b_covolume = 0.;

   double in_tol = 0.0000000000001,
          lambda_maxl_out = 0.,
          lambda_maxr_out = 0.,
          pstar = 0.,
          vstar = 0.;
   bool no_iter = true; // Had to change for test case to run
   int k = 0; // Tells you how many iterations were needed for convergence

   // cout << "CLM pre fortran function.\n";
   __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,
      &no_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&vstar,&k,&b_covolume);

   double d = std::max(std::abs(lambda_maxl_out), std::abs(lambda_maxr_out));

   if (isnan(d))
   {
      cout << "nij:\n";
      n_ij.Print(cout);
      cout << "Ui:\n";
      U_i.Print(cout);
      cout << "Uj:\n";
      U_j.Print(cout);
      cout << "in_rhol: " << in_rhol << ", ul: " << in_ul << ", el: " << in_el << ", pl: " << in_pl << endl;
      cout << "in_rhor: " << in_rhor << ", ur: " << in_ur << ", er: " << in_er << ", pr: " << in_pr << endl;
      MFEM_ABORT("NaN values returned by lambda max computation!\n");
   }

   return d;
   // return 1.;
}

template<int dim, int problem>
inline void ProblemDescription<dim, problem>::get_shocktube_vals(
   const int shocktube, double & _rhoL, double & _rhoR, double & _rhoLstar, 
   double & _rhoRstar, double & _vL, double & _vR, double & _vstar, double & _pL, 
   double & _pR, double & _pstar, double & _lambda1m, double & _lambda1p, 
   double & _lambda3)
{
   switch (problem)
   {
      case 6:
      {
         switch (shocktube)
         {
            case 1: // Sod
            {  
               _rhoL = 1.0;
               _rhoR = 0.125;
               _rhoLstar = 0.4263194281784952;
               _rhoRstar = 0.26557371170530708;
               _vL = 0.;
               _vR = 0.;
               _vstar = 0.92745262004894991;
               _pL = 1.0;
               _pR= 0.1;
               _pstar = 0.3031301780506468;
               _lambda1m = -1.183215956619923;
               _lambda1p = -0.07027281256118334;
               _lambda3 = 1.7521557320301779;
               return;
            }
            case 2: // Lax
            {
               _rhoL = 0.445;
               _rhoR = 0.5;
               _rhoLstar = 0.34456847418960945;
               _rhoRstar = 1.3040845320261998;
               _vL = 0.698;
               _vR = 0.;
               _vstar = 1.5287230266328840;
               _pL = 3.528;
               _pR= 0.571;
               _pstar = 2.4660979192073564;
               _lambda1m = -2.6335650740600323;
               _lambda1p = -1.6366974421005713;
               _lambda3 = 2.4793214809898405;
               return;
            }
            case 3: // Leblanc
            {
               _rhoL = 1.0;
               _rhoR = 0.001;
               _rhoLstar = 5.4079335349316249 * pow(10., -2);
               _rhoRstar = 3.9999980604299963 * pow(10., -3);
               _vL = 0.;
               _vR = 0.;
               _vstar = 0.62183867139173454;
               _pL = 2./30.;
               _pR= (2./3.) * pow(10., -10);
               _pstar = 5.1557792765096996 * pow(10., -4);
               _lambda1m = -1./3.;
               _lambda1p = 0.49578489518897934;
               _lambda3 = 0.82911836253346982;
               return;
            }
            default:
            {
               MFEM_ABORT("Bad number given for shocktube case number.");
            }
         }
      }
      case 7: // Saltzman
      {
         double _gamma = ProblemDescription<dim,problem>::gamma_func();
         _rhoL = 1.0;
         _rhoR = 1.0;
         _rhoLstar = 3.9992502342988532;
         _rhoRstar = 3.9992502342988532;
         _vL = 2.;
         _vR = 0.;
         _vstar = 1;
         _pL = (_gamma - 1.) * pow(10., -4); // change needed?
         _pR= (_gamma - 1.) * pow(10., -4); // change needed?
         _pstar = 1.3334833281256511;
         _lambda1m = 0.66658333854101581;
         _lambda1p = 0.66658333854101581; // No provided val?
         _lambda3 = 1.3334166614589844;
         return;
      }
      default:
      {
         MFEM_ABORT("Bad number given for problem case number.");
         return;
      }

   }
}

/* Explicit Instantiation */
template class ProblemDescription<1, 2>;
template class ProblemDescription<2, 2>;
template class ProblemDescription<3, 2>;

template class ProblemDescription<1, 1>;
template class ProblemDescription<2, 1>;
template class ProblemDescription<3, 1>;

template class ProblemDescription<1, 0>;
template class ProblemDescription<2, 0>;
template class ProblemDescription<3, 0>;

template class ProblemDescription<1, 3>;
template class ProblemDescription<2, 3>;
template class ProblemDescription<3, 3>;

template class ProblemDescription<1, 4>;
template class ProblemDescription<2, 4>;
template class ProblemDescription<3, 4>;

template class ProblemDescription<1, 5>;
template class ProblemDescription<2, 5>;
template class ProblemDescription<3, 5>;

template class ProblemDescription<1, 6>;
template class ProblemDescription<2, 6>;
template class ProblemDescription<3, 6>;

template class ProblemDescription<1, 7>;
template class ProblemDescription<2, 7>;
template class ProblemDescription<3, 7>;

template class ProblemDescription<1, 8>;
template class ProblemDescription<2, 8>;
template class ProblemDescription<3, 8>;

} // end ns hydrodynamics
} // end ns mfem