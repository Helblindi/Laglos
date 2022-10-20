/*
* The purpose of this class is to codify initial conditions corresponding
* to various test cases:
*
*     Case 0: Constant density and specific internal energy. 0 velocity.
*             In this case, we shouldn't see any values changing from
*             the prescribed initial conditions.
*     Case 1: The above but with normal vector with equal components for
*             the cell velocity.
*/

#include "mfem.hpp"

namespace mfem
{
namespace hydrodynamics
{
template<int problem, int dim>
class InitialValues
{
public:
   static inline double rho0(const Vector &x);        // density
   static inline double sv0(const Vector &x);         // specific volume = 1 / rho0
   static inline double sie0(const Vector &x);        // Specific internal energy
   static inline double ste0(const Vector &x);        // specific total energy
   static inline void v0(const Vector &x, Vector &v); // velocity
   static inline double rad(double x, double y);
   static inline double gamma_func(const Vector &x);
};

template<int problem, int dim>
inline double InitialValues<problem, dim>::sv0(const Vector &x)
{
   double val = rho0(x);
   assert(val != 0.);
   return 1./val;
}

template<int problem, int dim>
inline double InitialValues<problem, dim>::rho0(const Vector &x)
{
   switch (problem)
   {
      case 0: 
      case 1:
      {
         return 1.0;
      }
      default: 
      {
         MFEM_ABORT("Bad number given for problem id!"); 
         return 0.0;
      }
   }
}

// template<int problem, int dim>
// inline double InitialValues<problem, dim>::gamma_func(const Vector &x)
// {
//    switch (problem)
//    {
//       case 0: return 5.0 / 3.0;
//       case 1: return 1.4;
//       case 2: return 1.4;
//       case 3: return (x(0) > 1.0 && x(1) <= 1.5) ? 1.4 : 1.5;
//       case 4: return 5.0 / 3.0;
//       case 5: return 1.4;
//       case 6: return 1.4;
//       case 7: return 5.0 / 3.0;
//       default: MFEM_ABORT("Bad number given for problem id!"); return 0.0;
//    }
// }

// template<int problem, int dim>
// inline double InitialValues<problem, dim>::rad(double x, double y) 
// { 
//    return sqrt(x*x + y*y); 
// }

template<int problem, int dim>
inline double InitialValues<problem, dim>::sie0(const Vector &x)
{
   switch (problem)
   {
      case 0:
      case 1:
      {
         return 1.;
      }
      default: 
      {
         MFEM_ABORT("Bad number given for problem id!"); 
         return 0.0;
      }
   }
}

template<int problem, int dim>
inline double InitialValues<problem, dim>::ste0(const Vector &x)
{
   Vector v;
   v0(x,v);
   return sie0(x) + 0.5 * pow(v.Norml2(), 2);
}

template<int problem, int dim>
inline void InitialValues<problem, dim>::v0(const Vector &x, Vector &v)
{
   switch (problem)
   {
      case 0:
      {
         v = 0.;
         return;
      }
      case 1:
      {
         v = 1;
         v /= v.Norml2();
         return;
      }
      default: 
      {
         MFEM_ABORT("Bad number given for problem id!");
         return;
      }
   }
}

} // ns hydrodynamics

} // ns mfem
