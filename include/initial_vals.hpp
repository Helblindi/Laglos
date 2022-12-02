/*
* The purpose of this class is to codify initial conditions corresponding
* to various test cases:
*
*     Case 0: Constant density and specific internal energy. 0 velocity.
*             In this case, we shouldn't see any values changing from
*             the prescribed initial conditions.
*     Case 1: The above but with normal vector with equal components for
*             the cell velocity.
*     Case 2: Isentropic vortex as described in (6.1) with a stationary 
*             center at the origin.
*     Case 3: Isentropic vortex (see above) with moving center.
*     Case 4: Noh problem as described in (6.3)
*     Case 5: 1D Horizontal Movement on 2D mesh
*     Case 6: Shocktubes 1 - Sod, 2 - Lax, 3 - Leblanc
*     Case 7: Saltzman problem. See https://people.tamu.edu/~guermond/PUBLICATIONS/guermond_popov_Saavedra_JCP_2020.pdf.
*             Requires Neumann BC on right face and Dirichlet elsewhere.
*/

#include "mfem.hpp"
#include "compile_time_vals.h"
#include <cmath>


namespace mfem
{
namespace hydrodynamics
{
template<int problem, int dim>
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
   static inline double gamma_func(const int shocktube = 0);                                   // Gamma for problem
   static inline double IV_pressure(const Vector &x, const double & t); // Pressure (for Sod tube problem)
   static inline void get_shocktube_vals(const int shocktube, 
                                         double & _rhoL, 
                                         double & _rhoR, 
                                         double & _rhoLstar, 
                                         double & _rhoRstar, 
                                         double & _vL, 
                                         double & _vR, 
                                         double & _vstar, 
                                         double & _pL, 
                                         double & _pR, 
                                         double & _pstar, 
                                         double & _lambda1m, 
                                         double & _lambda1p, 
                                         double & _lambda3);
};

template<int problem, int dim>
inline double InitialValues<problem, dim>::rho0(const Vector &x, const double & t)
{
   switch (problem)
   {
      case 0: 
      case 1:
      {
         return 1.0;
      }
      case 2: // Isentropic vortex
      {
         const double beta = 5.; // vortex strength
         const double gamma = gamma_func(); 
         Vector center(2);
         center = 0.;
         Vector x_bar = x;
         x_bar -= center;
         const double r = x_bar.Norml2();
         const double dT = -1. * (exp(1 - pow(r,2)) * pow(beta, 2) * (gamma - 1)) / (8 * gamma * pow(M_PI, 2));
         const double rho = pow(1. + dT, 1/(gamma - 1));
         return rho;
      }
      case 3: // Isentropic vortex with moving center
      {
         const double beta = 5.; // vortex strength
         const double gamma = gamma_func(); 
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
      case 4: // Noh problem
      {
         const double norm_x = x.Norml2();
         if (norm_x <= t / 3)
         {
            return 16.;
         }
         else
         {
            return 1. + (t / norm_x);
         }
      }
      case 5: // 1D Horizontal motion
      {
         return x[0];
      }
      case 6: // Shock tube
      {
         Vector _n(2);
         _n[0] = cos(CompileTimeVals::rotation_angle);
         _n[1] = sin(CompileTimeVals::rotation_angle);
         double _x_tilde = x * _n;

         double _rhoL, _rhoR, _rhoLstar, _rhoRstar, _vL, _vR, _vstar, 
                _pL, _pR, _pstar, _lambda1m, _lambda1p, _lambda3;
         double _x0 = 0.000001; // Initial shock position
         double _gamma = gamma_func(CompileTimeVals::shocktube);
         get_shocktube_vals(CompileTimeVals::shocktube, _rhoL, _rhoR, _rhoLstar, _rhoRstar, 
                            _vL, _vR, _vstar, _pL, _pR, _pstar, _lambda1m, 
                            _lambda1p, _lambda3);
         double _cL = sqrt(_gamma * _pL / _rhoL);
         double _xi = (_x_tilde - _x0) / t;

         if (t == 0)
         {
            if (_x_tilde <= _x0) { return _rhoL; } else { return _rhoR; }
         }
         else 
         {
            if (_x_tilde <= _lambda1m * t)
            {
               return _rhoL;
            }
            else if (_x_tilde <= _lambda1p * t)
            {
               double val = _rhoL * pow(2/(_gamma + 1) + ((_gamma - 1)/((_gamma + 1) * _cL))*(_vL - _xi), 2/(_gamma - 1));
               return val;
            }
            else if (_x_tilde <= _vstar * t)
            {
               return _rhoLstar;
            }
            else if (_x_tilde <= _lambda3 * t)
            {
               return _rhoRstar;
            }
            else { return _rhoR; }
         }
      }
      case 7:{ // Saltzman TODO???
         if (t == 0) { return 1.; }
         else 
         {
            Vector _n(2);
            _n[0] = cos(CompileTimeVals::rotation_angle);
            _n[1] = sin(CompileTimeVals::rotation_angle);
            double _x_tilde = x * _n;

            double _rhoL, _rhoR, _rhoLstar, _rhoRstar, _vL, _vR, _vstar, 
                  _pL, _pR, _pstar, _lambda1m, _lambda1p, _lambda3;
            double _x0 = 0.000001; // Initial shock position
            double _gamma = gamma_func(CompileTimeVals::shocktube);
            get_shocktube_vals(CompileTimeVals::shocktube, _rhoL, _rhoR, _rhoLstar, _rhoRstar, 
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

template<int problem, int dim>
inline void InitialValues<problem, dim>::v0(const Vector &x, const double & t, Vector &v)
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
      case 2: // Isentropic vortex with 0 velocity
      {
         const double beta = 5.; // vortex strength
         const double gamma = gamma_func(); 
         Vector center(2);

         center = 0.;
         Vector x_bar = x;
         x_bar -= center;

         const double r = x_bar.Norml2();
         const double coeff = (exp((1 - pow(r,2)) / 2) * beta) / (2 * M_PI);

         v[0] = coeff * x_bar[1] * -1.;
         v[1] = coeff * x_bar[0];
         return;
      }
      case 3: // isentropic vortex with moving center
      {
         const double beta = 5.; // vortex strength
         const double gamma = gamma_func(); 
         Vector center(2);

         center[0] = 2. * t;
         center[1] = 0.;

         Vector x_bar = x;
         x_bar -= center;

         const double r = x_bar.Norml2();
         const double coeff = (exp((1 - pow(r,2)) / 2) * beta) / (2 * M_PI);

         v[0] = 2. + coeff * x_bar[1] * -1.;
         v[1] = coeff * x_bar[0];
         return;
      }
      case 4: // Noh Problem
      {
         const double norm_x = x.Norml2();
         if (norm_x <= t / 3)
         {
            v = 0.;
            return;
         }
         else
         {
            v = x;
            v /= norm_x;
            v *= -16.;
            return;
         }
      }
      case 5:
      {
         // v[0] = atan(x[0]);
         v[0] = x[0];
         v[1] = 0.;
         return;
      }
      case 6: // Shocktubes
      {
         Vector _n(2);
         _n[0] = cos(CompileTimeVals::rotation_angle);
         _n[1] = sin(CompileTimeVals::rotation_angle);
         double _x_tilde = x * _n;

         double _rhoL, _rhoR, _rhoLstar, _rhoRstar, _vL, _vR, _vstar, 
                _pL, _pR, _pstar, _lambda1m, _lambda1p, _lambda3;
         double _x0 = 0.000001; // Initial shock position
         double _gamma = gamma_func(CompileTimeVals::shocktube);
         get_shocktube_vals(CompileTimeVals::shocktube, _rhoL, _rhoR, _rhoLstar, _rhoRstar, 
                            _vL, _vR, _vstar, _pL, _pR, _pstar, _lambda1m, 
                            _lambda1p, _lambda3);
         double _cL = sqrt(_gamma * _pL / _rhoL);
         double _xi = (_x_tilde - _x0) / t;

         // our velocity will be some scaling of the normal vector
         v = _n;

         if (t == 0)
         {
            if (_x_tilde <= _x0) { v *= _vL; } else { v *= _vR; }
         }
         else 
         {
            if (_x_tilde <= _lambda1m * t)
            {
               v *= _vL;
            }
            else if (_x_tilde <= _lambda1p * t)
            {
               double val = (2 / (_gamma + 1)) * (_cL + ((_gamma -1)/2.)*_vL + _xi);
               v *= val;
            }
            else if (_x_tilde <= _vstar * t)
            {
               v *= _vstar;
            }
            else if (_x_tilde <= _lambda3 * t)
            {
               v *= _vstar;
            }
            else { v *= _vR; }
         }
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
            double _x0 = 0.000001; // Initial shock position
            double _gamma = gamma_func(CompileTimeVals::shocktube);
            get_shocktube_vals(CompileTimeVals::shocktube, _rhoL, _rhoR, _rhoLstar, _rhoRstar, 
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

template<int problem, int dim>
inline double InitialValues<problem, dim>::sv0(const Vector &x, const double & t)
{
   double val = rho0(x,t);
   assert(val != 0.);
   return 1./val;
}

template<int problem, int dim>
inline double InitialValues<problem, dim>::sie0(const Vector &x, const double & t)
{
   switch (problem)
   {
      case 2: // Isentropic vortex
      case 3: 
      {
         // Use pressure law to get specific internal energy
         const double gamma = gamma_func();
         double rho = rho0(x, t);
         double p = pow(rho,gamma);
         const double val = p / ((gamma - 1) * rho);
         return val;
      }
      case 4:
      {
         const double norm_x = x.Norml2();
         if (norm_x <= t / 3)
         {
            return 8.;
         }
         else
         {
            // return 0.5 * (1 + (t / norm_x));
            double p = 0.00001, gamma = gamma_func();
            return p / ((gamma - 1) * rho0(x,t));
         }
      }
      case 5: // 1D horizontal motion
      {
         return 0.;
      }
      case 6: // Sod tube, pressure is given
      {
         double _p = IV_pressure(x,t);
         double _gamma = gamma_func();
         double _rho = rho0(x,t);
         return _p / (_rho * (_gamma - 1));
      }
      case 7:
      {
         if (t == 0) { return pow(10, -4); }
         else {
            double _p = IV_pressure(x,t);
            double _gamma = gamma_func();
            double _rho = rho0(x,t);
            return _p / (_rho * (_gamma - 1));
         }
      }
      default: 
      {
         MFEM_ABORT("Bad number given for problem id!"); 
         return 0.0;
      }
   }
}

template<int problem, int dim>
inline double InitialValues<problem, dim>::ste0(const Vector &x, const double & t)
{
   switch (problem)
   {
      case 0:
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      case 7:
      default:
      {
         Vector v(2);
         v0(x,t,v);
         double sie = sie0(x, t);
         return sie + 0.5 * pow(v.Norml2(), 2);
      }
   }
   
}

template<int problem, int dim>
inline double InitialValues<problem, dim>::gamma_func(const int shocktube)
{
   switch (problem)
   {
      case 0:
      case 1:
      case 2:
      case 3:
      case 5:
      {
         return 7./5.;
      }
      case 6:
      {
         switch (shocktube)
         {
            case 0:
            case 1:                    // Sod
            case 2: { return 7./5.; }  // Lax
            case 3: { return 5./3.; }  // Leblanc
            default:
            {
               MFEM_ABORT("Bad number given for shocktube case number.");
               return 0.0;
            }
         }
      }
      case 7:
      case 4:
      {
         return 5./3.;
      }
      default:
      {
         MFEM_ABORT("Bad number given for problem id!"); 
         return 0.0;
      }

   }
}

// template<int problem, int dim>
// inline double InitialValues<problem, dim>::get_shocktube_vals()
// {

// }

template<int problem, int dim>
inline double InitialValues<problem, dim>::IV_pressure(const Vector & x, const double & t)
{
   switch (problem)
   {
      case 6:
      {
         Vector _n(2);
         _n[0] = cos(CompileTimeVals::rotation_angle);
         _n[1] = sin(CompileTimeVals::rotation_angle);
         double _x_tilde = x * _n;

         double _rhoL, _rhoR, _rhoLstar, _rhoRstar, _vL, _vR, _vstar, 
                _pL, _pR, _pstar, _lambda1m, _lambda1p, _lambda3;
         double _x0 = 0.000001; // Initial shock position
         double _gamma = gamma_func(CompileTimeVals::shocktube);
         get_shocktube_vals(CompileTimeVals::shocktube, _rhoL, _rhoR, _rhoLstar, _rhoRstar, 
                            _vL, _vR, _vstar, _pL, _pR, _pstar, _lambda1m, 
                            _lambda1p, _lambda3);
         double _cL = sqrt(_gamma * _pL / _rhoL);
         double _xi = (_x_tilde - _x0) / t;
         if (t == 0)
         {
            if (_x_tilde <= _x0) { return _pL; } else { return _pR; }
         }
         else 
         {
            if (x[0] <= _lambda1m * t)
            {
               return _pL;
            }
            else if (x[0] <= _lambda1p * t)
            {
               double val = _pL * pow( 2 / (_gamma + 1) + ((_gamma - 1)/((_gamma + 1) * _cL)) * (_vL - _xi), 2 * _gamma / (_gamma - 1) );
               return val;
            }
            else if (x[0] <= _vstar * t)
            {
               return _pstar;
            }
            else if (x[0] <= _lambda3 * t)
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
            Vector _n(2);
            _n[0] = cos(CompileTimeVals::rotation_angle);
            _n[1] = sin(CompileTimeVals::rotation_angle);
            double _x_tilde = x * _n;

            double _rhoL, _rhoR, _rhoLstar, _rhoRstar, _vL, _vR, _vstar, 
                  _pL, _pR, _pstar, _lambda1m, _lambda1p, _lambda3;
            double _x0 = 0.000001; // Initial shock position
            double _gamma = gamma_func(CompileTimeVals::shocktube);
            get_shocktube_vals(CompileTimeVals::shocktube, _rhoL, _rhoR, _rhoLstar, _rhoRstar, 
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

template<int problem, int dim>
inline void InitialValues<problem, dim>::get_shocktube_vals(
   const int shocktube, double & _rhoL, double & _rhoR, double & _rhoLstar, 
   double & _rhoRstar, double & _vL, double & _vR, double & _vstar, double & _pL, 
   double & _pR, double & _pstar, double & _lambda1m, double & _lambda1p, 
   double & _lambda3)
{
   switch (problem)
   {
      case 6:
      {
         switch (shocktube)
         {
            case 1: // Sod
            {  
               _rhoL = 1.0;
               _rhoR = 0.125;
               _rhoLstar = 0.4263194281784952;
               _rhoRstar = 0.26557371170530708;
               _vL = 0.;
               _vR = 0.;
               _vstar = 0.92745262004894991;
               _pL = 1.0;
               _pR= 0.1;
               _pstar = 0.3031301780506468;
               _lambda1m = -1.183215956619923;
               _lambda1p = -0.07027281256118334;
               _lambda3 = 1.7521557320301779;
               return;
            }
            case 2: // Lax
            {
               _rhoL = 0.445;
               _rhoR = 0.5;
               _rhoLstar = 0.34456847418960945;
               _rhoRstar = 1.3040845320261998;
               _vL = 0.698;
               _vR = 0.;
               _vstar = 1.5287230266328840;
               _pL = 3.528;
               _pR= 0.571;
               _pstar = 2.4660979192073564;
               _lambda1m = -2.6335650740600323;
               _lambda1p = -1.6366974421005713;
               _lambda3 = 2.4793214809898405;
               return;
            }
            case 3: // Leblanc
            {
               _rhoL = 1.0;
               _rhoR = 0.001;
               _rhoLstar = 5.4079335349316249 * pow(10., -2);
               _rhoRstar = 3.9999980604299963 * pow(10., -3);
               _vL = 0.;
               _vR = 0.;
               _vstar = 0.62183867139173454;
               _pL = 2./30.;
               _pR= (2./3.) * pow(10., -10);
               _pstar = 5.1557792765096996 * pow(10., -4);
               _lambda1m = -1./3.;
               _lambda1p = 0.49578489518897934;
               _lambda3 = 0.82911836253346982;
               return;
            }
            default:
            {
               MFEM_ABORT("Bad number given for shocktube case number.");
            }
         }
      }
      case 7: // Saltzman
      {
         double _gamma = gamma_func();
         _rhoL = 1.0;
         _rhoR = 1.0;
         _rhoLstar = 3.9992502342988532;
         _rhoRstar = 3.9992502342988532;
         _vL = 2.;
         _vR = 0.;
         _vstar = 1;
         _pL = (_gamma - 1.) * pow(10., -4); // change needed?
         _pR= (_gamma - 1.) * pow(10., -4); // change needed?
         _pstar = 1.3334833281256511;
         _lambda1m = 0.66658333854101581;
         _lambda1p = 0.66658333854101581; // No provided val?
         _lambda3 = 1.3334166614589844;
         return;
      }
      default:
      {
         MFEM_ABORT("Bad number given for problem case number.");
         return;
      }

   }
}

} // ns hydrodynamics

} // ns mfem
