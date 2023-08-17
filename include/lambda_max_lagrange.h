#ifndef LAMBDA_MAX_LAGRANGE
#define LAMBDA_MAX_LAGRANGE

#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;

namespace mfem
{
namespace hydrodymanics
{

double fz(const double p, const double pz, const double A, 
          const double B, const double C,  const double expoz)
{
   assert(p >= 0.);
   if (p < pz)
   {
      return C * (pow(p / pz, expoz) - 1.);
   }
   else
   {
      return (p - pz) * sqrt(A / (p + B));
   }
}

double phi(const double p, const double pL, const double uL, const double AL, 
           const double BL, const double CL,  const double expoL,
           const double pR, const double uR, const double AR, 
           const double BR, const double CR,  const double expoR)
{
   return fz(p, pL, AL, BL, CL, expoL) + fz(p, pR, AR, BR, CR, expoR) + uR - uL;
}

void init(const double tau, const double e, const double p, const double b_covolume, double & gamma, 
          double &a, double &A, double &B, double &C, double &expo, double &alpha)
{
   double x = tau - b_covolume;
   // Compute local gamma, gamma_z
   gamma = 1. + (p / e) * x;
   assert(gamma > 1.);
   // Compute local sound speed a_z
   a = tau * sqrt(gamma * p / x);
   // Other quantities
   A = 2. * x / (gamma + 1.);
   B = ((gamma - 1.) / (gamma + 1.)) * p;
   C = (2. * a * x) / (tau * (gamma - 1.)); // 4.1
   expo = 0.5 * (gamma - 1.) / gamma;
   // alpha = 1; // TODO
}

void lagrange_lambda_max(const double in_taul, const double in_ul, const double in_el, const double in_pl, 
                         const double in_taur, const double in_ur, const double in_er, const double in_pr, 
                         const double in_tol, const  double b_covolume, const bool WANT_ITERATION, 
                         double & lambda_maxl_out, double & lambda_maxr_out, double & pstar, int &k)
{
   double gammaL, aL, AL, BL, CL, expoL, alphaL;
   double gammaR, aR, AR, BR, CR, expoR, alphaR;
   init(in_taul, in_el, in_pl, b_covolume, gammaL, aL, AL, BL, CL, expoL, alphaL);
   init(in_taur, in_er, in_pr, b_covolume, gammaR, aR, AR, BR, CR, expoR, alphaR);

   // Compute phi(pmin), phi(pmax) from (4.2)  
   double p_min, tau_min, gamma_min, alpha_min, A_min, B_min, C_min;
   double p_max, tau_max, gamma_max, alpha_max, A_max, B_max, C_max;
   double gamma_lm, gamma_uM;
   double phi_pmin, phi_pmax;

   if (in_pl <= in_pr)
   {
      p_min = in_pl;
      tau_min = in_taul;
      gamma_min = gammaL;
      A_min = AL;
      B_min = BL;
      p_max = in_pr;
      tau_max = in_taur;
      gamma_max = gammaR;
      C_max = CR;
   }
   else
   {
      p_min = in_pr;
      tau_min = in_taur;
      gamma_min = gammaR;
      A_min = AR;
      B_min = BR;
      p_max = in_pl;
      tau_max = in_taul;
      gamma_max = gammaL;
      C_max = CL;
   }

   if (gammaL <= gammaR)
   {
      gamma_lm = gammaL;
      gamma_uM = gammaR;
   }
   else
   {
      gamma_lm = gammaR;
      gamma_uM = gammaL;
   }

   phi_pmin = phi(p_min, in_pl, in_ul, AL, BL, CL, expoL,
                  in_pr, in_ur, AR, BR, CR,  expoR);
   phi_pmax = phi(p_max, in_pl, in_ul, AL, BL, CL, expoL,
                  in_pr, in_ur, AR, BR, CR,  expoR);

   double phat_star = 0.;
   if (phi_pmin >= 0.)
   {
      // compute ptildestar from 5.5
      double ptilde_star = 0.;
      // set phatstar
      phat_star = max(p_min, ptilde_star);
   }
   else if (phi_pmax >= 0.)
   {

   }

}

}
}

#endif