#include "mfem.hpp"
#include "var-config.h"
#include "stdbool.h"
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>


using namespace mfem;
using namespace std;

/* Global variables */
double _gamma = 1.02, a = 0., b = 0.;

// Fortran subroutine from Lagrangian code
extern "C" {
   //  ,-- 2 Leading underscores to start
   //   | ,-- then the module name
   //   | |                                       ,-- then _MOD_
   //   | |                                       |    ,-- then the subroutine name
   //   V V                                       V    V
   void __arbitrary_eos_lagrangian_lambda_module_MOD_lambda_arbitrary_eos(
      double *in_rhol, double *in_ul, double *in_el, double *in_pl,
      double *in_rhor, double *in_ur, double *in_er, double *in_pr,
      double *in_tol, bool *want_iter,double *lambda_maxl_out,
      double *lambda_maxr_out, double *pstar, int *k, double* b_covolume, double* p_inf);

   double* __arbirary_eos_lagrangian_lambda_module_MOD_phi_(double *p);
}

double gamma_law_internal(double rho, double p);

int main(int argc, char **argv)
{
   std::ofstream out(argv[1]);
   
   std::string test_data = "/Fortran/lagrangian-mws/data_c";
   std::string data = std::string(LAGLOS_DIR) + test_data;
   std::ifstream infile(data);

   std::string line;

   std::string ind;
   int num_cases = 0, case_num = 0, k = 0;
   double rhoL, rhoR, uL, uR, pL, pR, tol;
   double tauL, tauR;
   double eL, eR;
   double a_vdw = 1.;
   double b_vdw = 1.;
   double p_inf = 0.;
   double gamma_vdw = 1.02;
   double gamma_ideal = 1.4;
   double lambdaL, lambdaR, pstar, vstar;

   bool next_num_cases = false;
   bool next_initial_data = false;
   bool next_tol = false;
   bool want_iter = true;

   while (std::getline(infile, line))
   {
      std::istringstream iss(line);
      iss >> ind;

      /* Check indicators */
      if (next_num_cases)
      {
         iss >> num_cases;
         next_num_cases = false;
      }
      else if (next_initial_data)
      {
         // Set initial data
         // rhol, rhor, ul, ur, pl, pr
         rhoL = stod(ind);
         iss >> rhoR
             >> uL
             >> uR
             >> pL
             >> pR;            

         tauL = 1./rhoL;
         tauR = 1./rhoR;

         next_initial_data = false;
         next_tol = true;
      }
      else if (next_tol)
      {
         // Set tol
         tol = stod(ind);

         if (case_num < 15)
         {
            // Run case
            _gamma = gamma_ideal;
            a = 0.;
            b = .1/(max(rhoL, rhoR));
         }
         else
         {
            _gamma = gamma_vdw;
            a = a_vdw;
            b = b_vdw;
         }

         eL = gamma_law_internal(rhoL, pL);
         eR = gamma_law_internal(rhoR, pR);
         __arbitrary_eos_lagrangian_lambda_module_MOD_lambda_arbitrary_eos(
            &rhoL, &uL, &eL, &pL, &rhoR, &uR, &eR, &pR, &tol, 
            &want_iter, &lambdaL, &lambdaR, &pstar, &k, &b, &p_inf);
         next_tol = false;

         // Output
         out << "===Case " << case_num << endl;
         out << "lambda_max=" << max(abs(lambdaL), abs(lambdaR)) << ", pstar=" << pstar << " k=" << k << endl;
         // out << "relative residual   " << residual << endl;
      }

      /* Set indicators for next line */
      if (ind == "===Number")
      {
         next_num_cases = true;
      }
      else if (ind == "===Case")
      {
         next_initial_data = true;

         iss >> case_num;
      }
   }

   return 0;
}


double gamma_law_internal(double rho, double p)
{
   cout << "gamma: " << _gamma << ", a: " << a << ", b: " << b << endl;
   return ((p + a * pow(rho,2)) * (1. - b*rho) / (rho * (_gamma - 1.))) - a * rho;
}
