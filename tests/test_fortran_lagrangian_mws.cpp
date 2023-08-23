#include "mfem.hpp"
#include "var-config.h"
#include "stdbool.h"
#include <cmath>
#include <sstream>
#include <fstream>


using namespace mfem;
using namespace std;

// Fortran subroutine from Lagrangian code
extern "C" {
   void __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      double *in_taul, double *in_ul, double *in_el, double *in_pl,
      double *in_taur, double *in_ur, double *in_er, double *in_pr,
      double *in_tol, bool *no_iter,double *lambda_maxl_out,
      double *lambda_maxr_out, double *pstar, int *k, double* b_covolume);
}

double gamma_law_internal(double b_covolume, double rho, double p, double gamma);

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
   double lambdaL, lambdaR, pstar, vstar;

   bool next_num_cases = false;
   bool next_initial_data = false;
   bool next_tol = false;
   bool no_iter = false;

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

         // Run case
         b_covolume = .1/(max(rhoL, rhoR));
         // b_covolume = 0.;
         gamma = 1.4;
         eL = gamma_law_internal(b_covolume, rhoL, pL, gamma);
         eR = gamma_law_internal(b_covolume, rhoR, pR, gamma);
         __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
            &tauL, &uL, &eL, &pL, &tauR, &uR, &eR, &pR, &tol, 
            &no_iter, &lambdaL, &lambdaR, &pstar, &k, &b_covolume);
         
         next_tol = false;

         // Output
         out << "===Case " << case_num << endl;
         out << "rhoL: " << rhoL << ", rhoR: " << rhoR << ", uL: " << uL << ", uR: " << uR << ", pL: " << pL << ", pR: " << pR << endl;
         out << "lambda_max=" << max(abs(lambdaL), abs(lambdaR)) << ", pstar=" << pstar << " k=" << k << endl;
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