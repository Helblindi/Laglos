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
   virtual double ComputeShearEnergy(const DenseMatrix &F) const {
      MFEM_ABORT("Must override.\n");
   } // virtual function, must be overridden;
   virtual void ComputeCauchyStress(const DenseMatrix &F, DenseMatrix &sigma) const {
      MFEM_ABORT("Must override.\n");
   } // virtual function, must be overridden;
   void ComputeIsotropicInvariants(const DenseMatrix &C, double &I1, double &I2, double &I3) const;

}; // end of ShearClosure class

class ShearClosureAnisotropic : public ShearClosure
{
protected:
   Vector mi; // fiber direction
public:
   ShearClosureAnisotropic(const Vector &_mi) : mi(_mi) {}
   virtual ~ShearClosureAnisotropic() {}
   void ComputeAnisotropicInvariants(const Vector &mi, const DenseMatrix &C, double &I4, double &I5) const;

}; // end of ShearClosureAnisotropic class

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
   double ComputeShearEnergy(const DenseMatrix &F) const override;
   void ComputeCauchyStress(const DenseMatrix &F, DenseMatrix &sigma) const override;
}; // end of ShearClosureTransverselyIsotropic class

} // namespace hydroLO
} // namespace mfem

#endif // SHEAR_CLOSURE