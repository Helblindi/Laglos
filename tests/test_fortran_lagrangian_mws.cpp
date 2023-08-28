#include "mfem.hpp"
#include "var-config.h"
#include "stdbool.h"
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>


using namespace mfem;
using namespace std;

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
      double *in_tol, bool *no_iter,double *lambda_maxl_out,
      double *lambda_maxr_out, double *pstar, int *k, double* b_covolume,
      double *a_vdw, double *b_vdw, double *gamma_vdw);

   double* __arbirary_eos_lagrangian_lambda_module_MOD_phi_(double *p);
}

double gamma_law_internal(double b_covolume, double rho, double p, double gamma);
double van_der_waals_internal(const double & a_vdw, const double & b_vdw, const double & gamma_vdw, 
                              const double & rho, const double & p);

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
   double b_covolume, gamma, eL, eR;
   double a_vdw = 1.;
   double b_vdw = 1.;
   double gamma_vdw = 1.02;
   double lambdaL, lambdaR, pstar, vstar;

   bool next_num_cases = false;
   bool next_initial_data = false;
   bool next_tol = false;
   bool no_iter = true;

   while (std::getline(infile, line))
   {
      std::istringstream iss(line);
      iss >> ind;
      // cout << "ind: " << ind << endl;

      /* Check indicators */
      if (next_num_cases)
      {
         // cout << "Grabbing num_cases\n";
         iss >> num_cases;
         next_num_cases = false;
      }
      else if (next_initial_data)
      {
         // cout << "Grabbing initial data for Case: " << case_num << endl;
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
         // cout << "Grabbing tol\n";
         // Set tol
         tol = stod(ind);

         if (case_num < 15)
         {
            // Run case
            b_covolume = .1/(max(rhoL, rhoR));
            gamma = 1.4;

            eL = gamma_law_internal(b_covolume, rhoL, pL, gamma);
            eR = gamma_law_internal(b_covolume, rhoR, pR, gamma);
         }
         else
         {
            eL = van_der_waals_internal(a_vdw, b_vdw, gamma_vdw, rhoL, pL);
            eR = van_der_waals_internal(a_vdw, b_vdw, gamma_vdw, rhoR, pR);
         }

         // cout << "rhoL: " << rhoL << ", rhoR: " << rhoR << ", uL: " << uL << ", uR: " << uR << ", pL: " << pL << ", pR: " << pR << endl;


         __arbitrary_eos_lagrangian_lambda_module_MOD_lambda_arbitrary_eos(
            &rhoL, &uL, &eL, &pL, &rhoR, &uR, &eR, &pR, &tol, 
            &no_iter, &lambdaL, &lambdaR, &pstar, &k, &b_covolume,
            &a_vdw, &b_vdw, &gamma_vdw);
         
         next_tol = false;

         /* Compute relative residual */
         // double phi_star = *__arbirary_eos_lagrangian_lambda_module_MOD_phi_(&pstar);
         // double phi_pl = __arbirary_eos_lagrangian_lambda_module_MOD_phi_(pL);
         // double phi_pr = __arbirary_eos_lagrangian_lambda_module_MOD_phi_(pR);
         // double residual = phi_pl / max( abs(phi_pl), abs(phi_pr));

         // Output
         out << "===Case " << case_num << endl;
         out << "lambda_max=" << max(abs(lambdaL), abs(lambdaR)) << ", pstar=" << pstar << " k=" << k << endl;
         // out << "relative residual   " << residual << endl;
      }

      /* Set indicators for next line */
      if (ind == "===Number")
      {
         next_num_cases = true;
         // cout << "we have the total num case:\n";
      }
      else if (ind == "===Case")
      {
         next_initial_data = true;
         // cout << "we have a case\n";

         iss >> case_num;
      }
   }

   return 0;
}


double gamma_law_internal(double b_covolume, double rho, double p, double gamma)
{
   return p*(1.-b_covolume*rho)/((gamma-1.)*rho);
}

double van_der_waals_internal(const double & a_vdw,
                              const double & b_vdw,
                              const double & gamma_vdw,
                              const double & rho, 
                              const double & p)
{
   return (p + a_vdw*rho*rho)*(1. - b_vdw*rho)/((gamma_vdw - 1.)*rho) - a_vdw*rho;
}