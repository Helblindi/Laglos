#include <iostream>

extern "C" {
   void __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      double *in_rhol, double *in_ul, double *in_el, double *in_pl,
      double *in_rhor, double *in_ur, double *in_er, double *in_pr,
      double *in_tol, bool *no_iter,double *lambda_maxl_out,
      double *lambda_maxr_out, double *pstar, int *k, double *b_covolume);
}


double gamma_law_internal(double rho, double p, double gamma, double b_covolume)
{
   double vv = p * (1 - b_covolume * rho) / ((gamma - 1) * rho);
   return vv;
}

int main()
{
   double in_rhol = 1., in_ul = 0., in_el, in_pl = 0.06666666666666666;
   double in_rhor = 0.001, in_ur = 0., in_er, in_pr = 0.00000000006666666666666666;
   double in_tol = 0.000000000000001;
   bool no_iter = false;
   double lambda_maxl_out, lambda_maxr_out, pstar;
   int k = 10;
   double b_covolume = 0.1 / std::max(in_rhol, in_rhor);
   in_el = gamma_law_internal(in_rhol, in_pl, 7./5., b_covolume);
   in_el = gamma_law_internal(in_rhor, in_pr, 7./5., b_covolume);
   

   __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,
      &no_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&k,&b_covolume);
   
   std::cout << "Lambda max: " << std::max(std::abs(lambda_maxl_out), std::abs(lambda_maxr_out)) << std::endl;

   return 0;
}