#ifndef SHEAR_CLOSURE
#define SHEAR_CLOSURE

#include "mfem.hpp"

namespace mfem {
namespace hydroLO {

// Shear closure model base class
class ShearClosure
{
public:
   virtual ~ShearClosure() {}
   virtual double ComputeShearEnergy(const DenseMatrix &F, const double &rho0=1.) const {
      MFEM_ABORT("Must override.\n");
   } // virtual function, must be overridden;
   virtual void ComputeCauchyStress(const DenseMatrix &F, DenseMatrix &sigma) const {
      MFEM_ABORT("Must override.\n");
   } // virtual function, must be overridden;
   void ComputeIsotropicInvariants(const DenseMatrix &C, double &J1, double &J2, double &J3) const;
   void ComputeIsotropicInvariantsReduced(const DenseMatrix &C, double &j1, double &j2) const;

}; // end of ShearClosure class

class ShearClosureAnisotropic : public ShearClosure
{
protected:
   Vector mi; // fiber direction
public:
   ShearClosureAnisotropic(const Vector &_mi) : mi(_mi) {}
   virtual ~ShearClosureAnisotropic() {}
   void ComputeAnisotropicInvariants(const Vector &mi, const DenseMatrix &C, double &j4, double &j5) const;

}; // end of ShearClosureAnisotropic class

class ShearClosureNeoHookean : public ShearClosure
{
private:
   double mu; // shear modulus
public:
   ShearClosureNeoHookean(const double &_mu) : mu(_mu) {}
   virtual ~ShearClosureNeoHookean() {}
   double ComputeShearEnergy(const DenseMatrix &F, const double &rho0=1.) const override;
   void ComputeCauchyStress(const DenseMatrix &F, DenseMatrix &sigma) const override;
}; // end of ShearClosureNeoHookean class

class ShearClosureMooneyRivlin : public ShearClosure
{
private:
   double mu; // shear modulus
public:
   ShearClosureMooneyRivlin(const double &_mu) : mu(_mu) {}
   virtual ~ShearClosureMooneyRivlin() {}
   double ComputeShearEnergy(const DenseMatrix &F, const double &rho0=1.) const override;
   void ComputeCauchyStress(const DenseMatrix &F, DenseMatrix &sigma) const override;
}; // end of ShearClosureMooneyRivlin class

class ShearClosureAortic : public ShearClosureAnisotropic
{
private:
   double w1, D1, A1, B1;
   double A2 = 1., B2 = 1.;
   DenseMatrix Gi;
   void ComputeAnisotropicInvariantsAortic(const DenseMatrix &C, double &i4, double &i5) const;
public:
   ShearClosureAortic(const Vector &_mi, const double &_w1, const double &_D1,
                      const double &_A1, const double &_B1);
   virtual ~ShearClosureAortic() {}
   double ComputeShearEnergy(const DenseMatrix &F, const double &rho0=1.) const override;
   void ComputeCauchyStress(const DenseMatrix &F, DenseMatrix &sigma) const override;
}; // end of ShearClosureAortic class

class ShearClosureTransverselyIsotropic : public ShearClosureAnisotropic
{
private:
   double E, EA, GA, nu;
   double lambda, mu, alpha, beta, gamma, m, n;
   void ComputeMaterialParameters();
public:
   ShearClosureTransverselyIsotropic(const Vector &_mi, const double &_E, const double &_EA,
                                    const double &_GA, const double &_nu)
      : ShearClosureAnisotropic(_mi), E(_E), EA(_EA), GA(_GA), nu(_nu) {}
   virtual ~ShearClosureTransverselyIsotropic() {}
   double ComputeShearEnergy(const DenseMatrix &F, const double &rho0=1.) const override;
   void ComputeCauchyStress(const DenseMatrix &F, DenseMatrix &sigma) const override;
}; // end of ShearClosureTransverselyIsotropic class

} // namespace hydroLO
} // namespace mfem

#endif // SHEAR_CLOSURE