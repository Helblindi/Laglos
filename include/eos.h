// EquationOfState.hpp
#pragma once
#include <cmath>
#include <stdexcept>

class EquationOfState {
public:
   virtual ~EquationOfState() = default;

   virtual double pressure(const double & rho, const double & e, const double & gamma) const = 0;
   virtual double sound_speed(const double & rho, const double & p, const double & gamma) const = 0;
   virtual double energy(const double & pressure, const double & density, const double & gamma) const = 0; // specific internal energy
};

class IdealGasEOS : public EquationOfState {
public:
   double pressure(const double &rho, const double &e, const double &gamma) const override {
      return (gamma - 1.0) * rho * e;
   }

   double sound_speed(const double &rho, const double &p, const double &gamma) const override {
      return sqrt(gamma * p / rho);
   }

   double energy(const double &pressure, const double &density, const double &gamma) const override {
      return pressure / ((gamma - 1.0) * density);
  }
};

class VanDerWaalsEOS : public EquationOfState {
private:
   double a, b;
public:
   VanDerWaalsEOS(double _a, double _b) : a(_a), b(_b) {}

   double pressure(const double &rho, const double &e, const double &gamma) const override {
      double denom = 1.0 - b * rho;
      if (denom <= 0.0) {
         std::cout << "rho: " << rho << " b*rho: " << b * rho << std::endl;
         throw std::runtime_error("Van der Waals: unphysical density (b * rho >= 1)");
      }
      return ((gamma - 1.0) * (rho * e + a * rho * rho)) / denom - a * rho * rho;
   }

   double sound_speed(const double &rho, const double &p, const double &gamma) const override {
      double denom = rho * (1.0 - b * rho);
      if (denom <= 0.0) {
         throw std::runtime_error("Van der Waals: unphysical density in sound speed");
      }
      // dp/drho derived assuming e is constant
      double num = gamma * (p + a * rho * rho);
      double _val = num / denom - 2. * a * rho;
      if (_val < 0.0) {
         throw std::runtime_error("Van der Waals: complex sound speed");
      }
      return sqrt(_val);
   }

   double energy(const double &pressure, const double &density, const double &gamma) const override {
      double denom = (gamma - 1.) * density;
      return (1. - b * density) * (pressure + a * density * density) / denom - a * density;
  }
};

class NobleAbelStiffenedGasEOS : public EquationOfState {
private:
   double b, q, p_inf;
public:
   NobleAbelStiffenedGasEOS(const double &_b, const double &_q, const double &_p_inf)
      : b(_b), q(_q),  p_inf(_p_inf) {}

   double pressure(const double &rho, const double &e, const double &gamma) const override {
      if (rho <= 0.0) {
         throw std::runtime_error("NASG: unphysical density (rho <= 1)");
      }
      double denom = 1. / rho - b;
      if (denom <= 0.0) {
         throw std::runtime_error("NASG: unphysical density (b * rho >= 1)");
      }
      return (gamma - 1.0) * (e - q) / denom - gamma * p_inf;
   }

   double sound_speed(const double &rho, const double &p, const double &gamma) const override {
      double denom = rho * (1.0 - b * rho);
      if (denom <= 0.0) {
         throw std::runtime_error("NASG: unphysical density in sound speed");
      }
      // dp/drho derived assuming e is constant
      double num = gamma * (p + p_inf);
      return sqrt(num / denom);
   }

   double energy(const double &pressure, const double &density, const double &gamma) const override {
      double denom = (gamma - 1.);
      if (denom <= 0.0) {
         throw std::runtime_error("NASG: unphysical gamma (gamma <= 1)");
      }
      return q + (pressure + gamma * p_inf) * (1. / density - b) / denom;
   }
};

class PolytropicEOS : public EquationOfState {
private:
   double K;
public:
   PolytropicEOS(const double &_K) : K(_K) {}

   double pressure(const double &rho, const double &e, const double &gamma) const override {
      return K * std::pow(rho, gamma);
   }

   double sound_speed(const double &rho, const double &p, const double &gamma) const override {
      return sqrt(gamma * p / rho);
   }

   double energy(const double &pressure, const double &density, const double &gamma) const override {
      return pressure / ((gamma - 1.0) * density);
   }
};
   