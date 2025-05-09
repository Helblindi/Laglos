// tests/test_eos.cpp
#include <iostream>
#include <cmath>
#include <memory>
#include <cassert>
#include "eos.h" // Adjust path to your eos.h as needed

bool approximatelyEqual(double a, double b, double tol = 1e-10) {
    return std::abs(a - b) < tol;
}

void test_ideal_gas() {
    double rho = 1.2;
    double e = 2.5;
    double gamma = 1.4;

    IdealGasEOS eos;
    double p = eos.pressure(rho, e, gamma);
    double e_back = eos.energy(p, rho, gamma);
    assert(approximatelyEqual(e, e_back));
}

void test_van_der_waals() {
    double rho = 1.2;
    double e = 2.5;
    double gamma = 1.4;
    double a = 0.1;
    double b = 0.05;

    VanDerWaalsEOS eos(a, b);
    double p = eos.pressure(rho, e, gamma);
    double e_back = eos.energy(p, rho, gamma);
    assert(approximatelyEqual(e, e_back, 1e-8)); // looser tolerance due to complexity
}

void test_noble_abel() {
    double rho = 1.5;
    double e = 3.0;
    double gamma = 1.4;
    double b = 0.1;
    double q = 0.5;
    double p_inf = 1.0;

    NobleAbelStiffenedGasEOS eos(b, q, p_inf);
    double p = eos.pressure(rho, e, gamma);
    double e_back = eos.energy(p, rho, gamma);
    assert(approximatelyEqual(e, e_back, 1e-8));
}

void test_polytropic() {
    double rho = 1.3;
    double gamma = 1.4;
    double K = 0.5;

    PolytropicEOS eos(K);
    double p = eos.pressure(rho, 0.0, gamma); // e not used
    double e = eos.energy(p, rho, gamma);
    double p_back = eos.pressure(rho, e, gamma);
    assert(approximatelyEqual(p, p_back));
}

int main() {
    test_ideal_gas();
    test_van_der_waals();
    test_noble_abel();
    test_polytropic();

    std::cout << "All EOS tests passed!" << std::endl;
    return 0;
}
