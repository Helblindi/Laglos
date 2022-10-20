#include "mfem.hpp"
#include "laglos_solver.hpp"
#include <cassert>
#include <fstream>
#include <sstream>

using namespace std;

const int dim = 1;
static double b_covolume = 0;

void test_lambda_max_computation();
double gamma_law_internal(double rho, double p, double b_covolume);

int main()
{
   test_lambda_max_computation();
   return 0;
}

void test_lambda_max_computation()
{
   ifstream infile("../Fortran/data_c"); // Import file
   string line;
   string rhol_str, rhor_str, ul_str, ur_str, pl_str, pr_str;
   double rhol, rhor, ul, ur, pl, pr;
   double el, er;
   double tol;
   double lambda_max;
   int num_cases = 0;
   string val, next_step;
   mfem::Vector n(dim), Ul_v(dim+3), Ur_v(dim+3);
   n  = 1.;
   int _case = 0;

   while(getline(infile, line))
   {
      istringstream iss(line);
      // cout << "line: " << line << endl;

      // Check conditions for current line
      if (next_step == "num_cases") { 
         iss >> val;
         num_cases = stoi(val);
         next_step = "NA";
      } else if (next_step == "read_vals") {
         // cout << "Setting doubles\n";
         iss >> rhol_str >> rhor_str >> ul_str >> ur_str >> pl_str >> pr_str;
         rhol = stod(rhol_str);
         rhor = stod(rhor_str);
         ul = stod(ul_str);
         ur = stod(ur_str);
         pl = stod(pl_str);
         pr = stod(pr_str);

         double b_covolume = 0.1/max(rhol,rhor); // Needed for gamma_law_internal

         el = gamma_law_internal(rhol,pl, b_covolume);
         er = gamma_law_internal(rhor,pr, b_covolume);

         Ul_v[0] = rhol;
         Ul_v[1] = ul;
         Ul_v[2] = pl;
         Ul_v[3] = el;

         Ur_v[0] = rhor;
         Ur_v[1] = ur;
         Ur_v[2] = pr;
         Ur_v[3] = er;

         // cout << "rhol: " << rhol
         //      << ", ul: " << ul
         //      << ", pl: " << pl
         //      << ", el: " << el << endl;
         // cout << "rhor: " << rhor
         //      << ", ur: " << ur
         //      << ", pr: " << pr
         //      << ", er: " << er << endl;

         
         // Call case
         lambda_max = mfem::hydrodynamics::LagrangianLOOperator<dim>::compute_lambda_max(Ul_v, Ur_v, n, "testing");
         cout << "case: " << _case << ", lambda max: " << lambda_max << endl;

         next_step = "read_tol";

      } else if (next_step == "read_tol") {
         iss >> val;
         tol = stod(val);
         next_step = "NA";
      }
      else {
         // cout << "else statement\n";
         iss >> val;
      }

      // cout << endl << "val: " << val << endl << endl;;

      if (val == "===Number") { next_step = "num_cases"; }
      else if (val == "===Case") { 
         iss >> val;
         _case = stoi(val);
         next_step = "read_vals"; 
      }
      
   }
   return;
}

double gamma_law_internal(double rho, double p, double b_covolume)
{
   double vv = (p*(1-b_covolume*rho)) / ((mfem::hydrodynamics::LagrangianLOOperator<dim>::gamma - 1.) * rho); 

   // cout << "b_covolume: "
   //      << b_covolume
   //      << ", rho: " 
   //      << rho 
   //      << ", p: " 
   //      << p 
   //      << ", gamma: "
   //      << mfem::hydrodynamics::LagrangianLOOperator<dim>::gamma
   //      << ", vv: " 
   //      << vv
   //      << endl;

   return vv;
}
   //  vv = p*(1-b_covolume*rho)/((gamma-1)*rho)
