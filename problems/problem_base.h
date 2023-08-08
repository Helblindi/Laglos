#include "mfem.hpp"
#include <cmath>

using namespace mfem;
using namespace std;

namespace mfem
{
namespace hydrodynamics
{

template<int dim>
class ProblemBase
{
public:
   double a = 0., b = 0., gamma = 0.;
   bool distort_mesh = false;
   bool known_exact_solution = false;
   // TODO: BC specifier

   // ProblemDescription
   double internal_energy(const Vector &U)
   {
      const double &rho = 1./U[0];
      const double &e = specific_internal_energy(U);
      return rho * e;
   }

   virtual double specific_internal_energy(const Vector &U)
   {
      const Vector  v   = velocity(U);
      const double &E   = U[dim + 1]; // specific total energy
      return E - 0.5 * pow(v.Norml2(), 2);
   }

   double gamma_func() { return gamma; }

   virtual double pressure(const Vector &U) = 0; // virtual function, must be overridden

   // InitialState
   virtual void sv0(const Vector &x, const double & t)
   {
      double val = rho0(x,t);
      assert(val != 0.);
      return 1./val;
   }

   virtual void ste0(const Vector &x, const double & t)
   {
      Vector v(dim);
      v0(x,t,v);
      return sie0(x, t) + 0.5 * pow(v.Norml2(), 2);
   }

   virtual double rho0(const Vector &x, const double & t) = 0; // virtual function, must be overridden
   virtual void v0(const Vector &x, const double & t, Vector &v) = 0; // virtual function, must be overridden
   virtual void sie0(const Vector &x, const double & t) = 0; // virtual function, must be overridden
   
}; // End ProblemBase

} // ns hydrodynamics
} // ns mfem