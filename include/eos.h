// EquationOfState.hpp
#pragma once
#include <cmath>
#include <stdexcept>

class EquationOfState {
public:
   virtual ~EquationOfState() = default;

   virtual double pressure(double rho, double e, double gamma) const = 0;
   virtual double sound_speed(double rho, double e, double gamma) const = 0;
   virtual double energy(double pressure, double density, double gamma) const = 0; // specific internal energy
};

class IdealGasEOS : public EquationOfState {
public:
   double pressure(double rho, double e, double gamma) const override {
      return (gamma - 1.0) * rho * e;
   }

   double sound_speed(double rho, double e, double gamma) const override {
      // double p = pressure(rho, e, gamma);
      // return std::sqrt(gamma * p / rho);
      throw std::runtime_error("Sound speed not implemented");
      return 0.;
   }

   double energy(double pressure, double density, double gamma) const override {
      return pressure / ((gamma - 1.0) * density);
  }
};

class VanDerWaalsEOS : public EquationOfState {
private:
   double a, b;
public:
   VanDerWaalsEOS(double _a, double _b) : a(_a), b(_b) {}

   double pressure(double rho, double e, double gamma) const override {
      double denom = 1.0 - b * rho;
      if (denom <= 0.0) {
         throw std::runtime_error("Van der Waals: unphysical density (b * rho >= 1)");
      }
      return ((gamma - 1.0) * (rho * e + a * rho * rho)) / denom - a * rho * rho;
   }

   double sound_speed(double rho, double e, double gamma) const override {
      // double p = pressure(rho, e, gamma);
      // double denom = (1.0 - b * rho);
      // if (denom <= 0.0) {
      //    throw std::runtime_error("Van der Waals: unphysical density in sound speed");
      // }
      // // dp/drho derived assuming e is constant
      // double dp_drho = (gamma - 1.0) * e * (1.0 + b * rho / denom) / denom / denom - 2.0 * a * rho;
      // return std::sqrt(dp_drho);
      throw std::runtime_error("Sound speed not implemented");
      return 0.;
   }

   double energy(double pressure, double density, double gamma) const override {
      double denom = (gamma - 1.) * density;
      return (1. - b * density) * (pressure + a * density * density) / denom - a * density;
  }
};

class NobleAbelStiffenedGasEOS : public EquationOfState {
private:
   double q, b, p_inf;
public:
   NobleAbelStiffenedGasEOS(const double &_q,const double &_b, const double &_p_inf)
      : q(_q), b(_b), p_inf(_p_inf) {}

   double pressure(double rho, double e, double gamma) const override {
      if (rho <= 0.0) {
         throw std::runtime_error("NASG: unphysical density (rho <= 1)");
      }
      double denom = 1. / rho - b;
      if (denom <= 0.0) {
         throw std::runtime_error("NASG: unphysical density (b * rho >= 1)");
      }
      return (gamma - 1.0) * (e - q) / denom - gamma * p_inf;
   }

   double sound_speed(double rho, double e, double gamma) const override {
      // double p = pressure(rho, e, gamma);
      // double term1 = gamma * (p - p_inf + (a * rho * rho) / (1.0 - b * rho)) / rho;
      // double term2 = (a * rho) / ((1.0 - b * rho) * (1.0 - b * rho));
      // return std::sqrt(term1 + term2);
      throw std::runtime_error("Sound speed not implemented");
      return 0.;
   }

   double energy(double pressure, double density, double gamma) const override {
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

   double pressure(double rho, double e, double gamma) const override {
      return K * std::pow(rho, gamma);
   }

   double sound_speed(double rho, double e, double gamma) const override {
      // double p = pressure(rho, e, gamma);
      // return std::sqrt(gamma * p / rho);
      throw std::runtime_error("Sound speed not implemented");
      return 0.;
   }

   double energy(double pressure, double density, double gamma) const override {
      return pressure / ((gamma - 1.0) * density);
   }
};
   