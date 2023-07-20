/*
* The purpose of this class is to codify initial conditions corresponding
* to various test cases:
*
*     Case 0: Constant density and specific internal energy. 0 velocity.
*             In this case, we shouldn't see any values changing from
*             the prescribed initial conditions.
*     Case 1: Sod Shocktube 1D
*     Case 2: Lax Shocktube 1D ---TODO---
*     Case 3: Leblanc Shocktube 1D ---TODO---
*     Case 4: Noh problem as described in (6.3) ---TODO---
*     Case 5: Isentropic vortex as described in (6.1) with a stationary 
*             center at the origin. ---TODO---
*     Case 6: Isentropic vortex (see above) with moving center. ---TODO---
*     Case 7: Saltzman problem. See https://people.tamu.edu/~guermond/PUBLICATIONS/guermond_popov_Saavedra_JCP_2020.pdf.
*             Requires Neumann BC on right face and Dirichlet elsewhere. ---TODO---
*/

#include "mfem.hpp"
#include "compile_time_vals.h"
#include "riemann1D.hpp"
#include <cmath>


namespace mfem
{
namespace hydrodynamics
{
template<int dim, int problem>
class InitialValues
{
public:
   const static int shocktube = CompileTimeVals::shocktube;
   const static int shocktube_rotation = CompileTimeVals::rotation_angle;
   static inline double rho0(const Vector &x, const double & t);        // density
   static inline void v0(const Vector &x, const double & t, Vector &v); // velocity
   static inline double sv0(const Vector &x, const double & t);         // specific volume = 1 / rho0
   static inline double sie0(const Vector &x, const double & t);        // Specific internal energy
   static inline double ste0(const Vector &x, const double & t);        // specific total energy
   static inline double IV_pressure(const Vector &x, const double & t); // Pressure (for Sod tube problem)
};

template<int dim, int problem>
inline double InitialValues<dim, problem>::rho0(const Vector &x, const double & t)
{
   switch (problem)
   {
      case 0: 
      {
         assert(dim == 1);
         double x0 = 0.1, x1 = 0.3;
         if (x[0] - t < x0)
         {
            return 1.;
         }
         else if (x[0] - t < x1)
         {
            double _rho = 1 + pow(2.,6) * (1. / pow(x1-x0,6)) * pow(x[0]-t-x0,3) * pow(x1-x[0]+t,3);
            return _rho;
         }
         else 
         {
            return 1.;
         }
      }
      case 1: return (x(0) < 0.5) ? 1.0 : 0.125; // Sod shocktube
      case 2:
      case 3:
      case 4: // NOH
      {
         // const double norm_x = x.Norml2();
         // if (norm_x <= t / 3)
         // {
         //    return 16.;
         // }
         // else
         // {
         //    return 1. + (t / norm_x);
         // }
         MFEM_ABORT("Not Implemented.\n"); 
      }
      
      case 5: // Isentropic vortex
      {
         const double beta = 5.; // vortex strength
         const double gamma = ProblemDescription<dim,problem>::gamma_func(); 
         Vector center(2);
         center = 0.;
         Vector x_bar = x;
         x_bar -= center;
         const double r = x_bar.Norml2();
         const double dT = -1. * (exp(1 - pow(r,2)) * pow(beta, 2) * (gamma - 1)) / (8 * gamma * pow(M_PI, 2));
         const double rho = pow(1. + dT, 1./(gamma - 1.));
         return rho;
      }
      case 6: // Isentropic vortex with moving center
      {
         const double beta = 5.; // vortex strength
         const double gamma = ProblemDescription<dim,problem>::gamma_func(); 
         Vector center(2);
         
         center[0] = 2. * t;
         center[1] = 0.;

         Vector x_bar = x;
         x_bar -= center;
         const double r = x_bar.Norml2();
         const double dT = -1. * (exp(1 - pow(r,2)) * pow(beta, 2) * (gamma - 1)) / (8 * gamma * pow(M_PI, 2));
         const double rho = pow(1. + dT, 1/(gamma - 1));
         return rho;
      }
      case 7:{
         if (t == 0) { return 1.; }
         else 
         {
            Vector _n(dim);
            _n[0] = cos(CompileTimeVals::rotation_angle);
            _n[1] = sin(CompileTimeVals::rotation_angle);
            double _x_tilde = x * _n;

            double _rhoL, _rhoR, _rhoLstar, _rhoRstar, _vL, _vR, _vstar, 
                  _pL, _pR, _pstar, _lambda1m, _lambda1p, _lambda3;
            double _x0 = 0.0000000001; // Initial shock position
            double _gamma = ProblemDescription<dim,problem>::gamma_func(CompileTimeVals::shocktube);
            ProblemDescription<dim,problem>::get_shocktube_vals(CompileTimeVals::shocktube, _rhoL, _rhoR, _rhoLstar, _rhoRstar, 
                              _vL, _vR, _vstar, _pL, _pR, _pstar, _lambda1m, 
                              _lambda1p, _lambda3);
            double _cL = sqrt(_gamma * _pL / _rhoL);
            double _xi = (_x_tilde - _x0) / t;

            // Saltzman TODO
            if (_x_tilde <= _lambda1p * t)
            {
               return _rhoL;
            }
            else if (_x_tilde <= _lambda3 * t)
            {
               return _rhoRstar;
            }
            else
            {
               return _rhoR;
            }
         }
         return 0.;
      }
      default: 
      {
         MFEM_ABORT("Bad number given for problem id!"); 
         return 0.0;
      }
   }
}

template<int dim, int problem>
inline void InitialValues<dim, problem>::v0(const Vector &x, const double & t, Vector &v)
{
   switch (problem)
   {
      case 0:
      {
         assert(dim==1);
         v = 1.;
         return;
      }
      case 1: v = 0.0; break; // Sod
      case 2:
      case 3:
      case 4: // Noh Problem
      {
         // const double norm_x = x.Norml2();
         // if (norm_x <= t / 3)
         // {
         //    v = 0.;
         //    return;
         // }
         // else
         // {
         //    v = x;
         //    v /= norm_x;
         //    v *= -16.;
         //    return;
         // }
         MFEM_ABORT("Not Implemented.\n"); 
      }
      case 5: // Isentropic vortex with 0 velocity
      {
         const double beta = 5.; // vortex strength
         const double gamma = ProblemDescription<dim,problem>::gamma_func(); 
         Vector center(2);

         center = 0.;
         Vector x_bar = x;
         x_bar -= center;

         const double r = x_bar.Norml2();
         const double coeff = (exp((1. - pow(r,2)) / 2.) * beta) / (2. * M_PI);

         v[0] = coeff * x_bar[1] * -1.;
         v[1] = coeff * x_bar[0];
         return;
      }
      case 6: // isentropic vortex with moving center
      {
         const double beta = 5.; // vortex strength
         const double gamma = ProblemDescription<dim,problem>::gamma_func(); 
         Vector center(2);

         center[0] = 2. * t;
         center[1] = 0.;

         Vector x_bar = x;
         x_bar -= center;

         const double r = x_bar.Norml2();
         const double coeff = (exp((1. - pow(r,2)) / 2.) * beta) / (2. * M_PI);

         v[0] = 2. + coeff * x_bar[1] * -1.;
         v[1] = coeff * x_bar[0];
         return;
      }
      case 7: // Saltzman
      {
         if (t == 0) { v = 0.; }
         else 
         {
            Vector _n(2);
            _n[0] = cos(CompileTimeVals::rotation_angle);
            _n[1] = sin(CompileTimeVals::rotation_angle);
            double _x_tilde = x * _n;

            double _rhoL, _rhoR, _rhoLstar, _rhoRstar, _vL, _vR, _vstar, 
                  _pL, _pR, _pstar, _lambda1m, _lambda1p, _lambda3;
            double _x0 = 0.0000000001; // Initial shock position
            double _gamma = ProblemDescription<dim,problem>::gamma_func(CompileTimeVals::shocktube);
            ProblemDescription<dim,problem>::get_shocktube_vals(CompileTimeVals::shocktube, _rhoL, _rhoR, _rhoLstar, _rhoRstar, 
                              _vL, _vR, _vstar, _pL, _pR, _pstar, _lambda1m, 
                              _lambda1p, _lambda3);
            double _cL = sqrt(_gamma * _pL / _rhoL);
            double _xi = (_x_tilde - _x0) / t;

            // our velocity will be some scaling of the normal vector
            v = _n;

            // Saltzman TODO
            if (_x_tilde <= _lambda1p * t)
            {
               v *= _vL;
            }
            else if (_x_tilde <= _lambda3 * t)
            {
               v *= _vstar;
            }
            else
            {
               v *= _vR;
            }
         }
         return;
      }
      default: 
      {
         MFEM_ABORT("Bad number given for problem id!");
         return;
      }
   }
}

template<int dim, int problem>
inline double InitialValues<dim, problem>::sv0(const Vector &x, const double & t)
{
   double val = rho0(x,t);
   assert(val != 0.);
   return 1./val;
}

template<int dim, int problem>
inline double InitialValues<dim, problem>::sie0(const Vector &x, const double & t)
{
   switch (problem)
   {
      case 1: return (x(0) < 0.5) ? 1.0 / rho0(x, t) / (ProblemDescription<dim,problem>::gamma_func() - 1.0) // Sod
                        : 0.1 / rho0(x, t) / (ProblemDescription<dim,problem>::gamma_func() - 1.0);
      case 4: // Noh
      {
         const double norm_x = x.Norml2();
         if (norm_x <= t / 3)
         {
            return 8.;
         }
         else
         {
            // return 0.5 * (1 + (t / norm_x));
            double p = 0.00001, gamma = ProblemDescription<dim,problem>::gamma_func();
            return p / ((gamma - 1) * rho0(x,t));
         }
      }
      case 7:
      {
         if (t == 0) { return pow(10, -4); }
         else {
            double _p = IV_pressure(x,t);
            double _gamma = ProblemDescription<dim,problem>::gamma_func();
            double _rho = rho0(x,t);
            return _p / (_rho * (_gamma - 1.));
         }
      }
      /* In all other cased use pressure law to get specific internal energy */
      case 0:
      case 2:
      case 3: 
      case 5: 
      case 6:
      {
         const double gamma = ProblemDescription<dim,problem>::gamma_func();
         double rho = rho0(x, t);
         double p = pow(rho,gamma);
         const double val = p / ((gamma - 1) * rho);
         return val;
      }
      default: 
      {
         MFEM_ABORT("Bad number given for problem id!"); 
         return 0.0;
      }
   }
}

template<int dim, int problem>
inline double InitialValues<dim, problem>::ste0(const Vector &x, const double & t)
{
   switch (problem)
   {
      default:
      {
         Vector v(dim);
         v0(x,t,v);
         return sie0(x, t) + 0.5 * pow(v.Norml2(), 2);
      }
   }
   
}


// TODO: Get rid of this function
template<int dim, int problem>
inline double InitialValues<dim, problem>::IV_pressure(const Vector & x, const double & t)
{
   switch (problem)
   {
      case 6:
      {
         double _rhoL, _rhoR, _rhoLstar, _rhoRstar, _vL, _vR, _vstar, 
                _pL, _pR, _pstar, _lambda1m, _lambda1p, _lambda3;
         double _x0 = 0.5; // Initial shock position
         double _gamma = ProblemDescription<dim,problem>::gamma_func(CompileTimeVals::shocktube);
         ProblemDescription<dim,problem>::get_shocktube_vals(CompileTimeVals::shocktube, _rhoL, _rhoR, _rhoLstar, _rhoRstar, 
                            _vL, _vR, _vstar, _pL, _pR, _pstar, _lambda1m, 
                            _lambda1p, _lambda3);
         double _cL = sqrt(_gamma * _pL / _rhoL);
         
         if (t == 0)
         {
            if (x[0] <= _x0) { return _pL; } else { return _pR; }
         }
         else 
         {
            double _xi = (x[0] - _x0) / t;
            if (_xi <= _lambda1m)
            {
               return _pL;
            }
            else if (_xi <= _lambda1p)
            {
               double val = _pL * pow( 2 / (_gamma + 1) + ((_gamma - 1)/((_gamma + 1) * _cL)) * (_vL - _xi), 2 * _gamma / (_gamma - 1) );
               return val;
            }
            else if (_xi <= _vstar)
            {
               return _pstar;
            }
            else if (_xi <= _lambda3)
            {
               return _pstar;
            }
            else { return _pR; }
         }
      }
      case 7:
      {
         if (t == 0) { return 1.; } // TODO: Change pressure
         else {
            Vector _n(dim);
            _n[0] = cos(CompileTimeVals::rotation_angle);
            if (dim > 1)
            {
               _n[1] = sin(CompileTimeVals::rotation_angle);
            }
            double _x_tilde = x * _n;

            double _rhoL, _rhoR, _rhoLstar, _rhoRstar, _vL, _vR, _vstar, 
                  _pL, _pR, _pstar, _lambda1m, _lambda1p, _lambda3;
            double _x0 = 0.0000000001; // Initial shock position
            double _gamma = ProblemDescription<dim,problem>::gamma_func(CompileTimeVals::shocktube);
            ProblemDescription<dim,problem>::get_shocktube_vals(CompileTimeVals::shocktube, _rhoL, _rhoR, _rhoLstar, _rhoRstar, 
                              _vL, _vR, _vstar, _pL, _pR, _pstar, _lambda1m, 
                              _lambda1p, _lambda3);
            double _cL = sqrt(_gamma * _pL / _rhoL);
            double _xi = (_x_tilde - _x0) / t;

            // Saltzman TODO
            if (_x_tilde <= _lambda1p * t)
            {
               return _pL;
            }
            else if (_x_tilde <= _lambda3 * t)
            {
               return _pstar;
            }
            else
            {
               return _pR;
            }
         }
         
      }
      default:
      {
         MFEM_ABORT("Bad value given. IV_pressure should only be used for Sod tube!");
         return 0.0;
      }
   }
}

} // ns hydrodynamics

} // ns mfem
