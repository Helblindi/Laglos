#include "mfem.hpp"
#include "var-config.h"
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
      double *lambda_maxr_out, double *pstar, int *k);
}

// Fortran subroutine from Eulerian code
extern "C" {
   void __arbitrary_eos_lambda_module_MOD_eulerian_lambda_arbitrary_eos(
      double *in_rhol, double *in_ul, double *in_el, double *in_pl,
      double *in_rhor, double *in_ur, double *in_er, double *in_pr,
      double *in_tol, bool *no_iter,double *lambda_maxl_out,
      double *lambda_maxr_out, double *pstar, double *vstar, int *k, double *b_covolume);
}

int test_lagrangian_lambda_max();
int test_eulerian_lambda_max();

int main()
{
   int d = test_eulerian_lambda_max();
   // d += test_lagrangian_lambda_max();

   return d;
}

int test_eulerian_lambda_max()
{
   std::string test_data = "/Fortran/eulerian-mws/data_c";
   std::string data = std::string(LAGLOS_DIR) + test_data;
   std::ifstream infile(data);

   std::string line;

   double a;
   std::string ind;
   while (std::getline(infile, line))
   {
      std::istringstream iss(line);
      iss >> ind;
      cout << "ind: " << ind << endl;

      if (iss >> a)
      {
         cout << "wrote to double\n";
      }

      // int a, b;
      // if (!(iss >> a >> b)) { break; } // error
      cout << "string: " << iss.str() << endl;
      // process pair (a,b)
   }

   return 0;
}

int test_lagrangian_lambda_max()
{
   std::string test_data = "/Fortran/lagrangian-mws/data_c";
   std::string data = std::string(LAGLOS_DIR) + test_data;
   std::ifstream infile(data);

   std::string line;
   double num_cases = 0;
   while (std::getline(infile, line))
   {
      std::istringstream iss(line);
      // int a, b;
      // if (!(iss >> a >> b)) { break; } // error
      cout << "string: " << iss.str() << endl;
      // process pair (a,b)
   }

   return 0;
}