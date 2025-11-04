#include "shear_closure.hpp"
using namespace std;

namespace mfem {
namespace hydroLO {

void ShearClosure::ComputeIsotropicInvariants(const DenseMatrix &C, double &I1, double &I2, double &I3) const
{
   I1 = C.Trace();
   DenseMatrix C2(3);
   mfem::Mult(C, C, C2);
   // MFEM_WARNING("SHEAR ENERGY CHECK");
   I2 = 0.5 * (I1 * I1 - C2.Trace());
   I3 = C.Det();
}

void ShearClosureAnisotropic::ComputeAnisotropicInvariants(const Vector &mi, const DenseMatrix &C, double &I4, double &I5) const
{
   // MFEM_WARNING("SHEAR ENERGY CHECK");
   Vector Cmi(3), C2mi(3);
   // cout << "C:\n";
   // C.Print(cout);
   // cout << "mi:\n";
   // mi.Print(cout);
   C.Mult(mi, Cmi);
   I4 = mfem::InnerProduct(mi, Cmi);

   DenseMatrix C2(3);
   mfem::Mult(C, C, C2);
   C2.Mult(mi, C2mi);
   I5 = mfem::InnerProduct(mi, C2mi);
}

void ShearClosureTransverselyIsotropic::ComputeMaterialParameters()
{
   assert(1 + nu != 0.);
   m = 1. - nu - 2. * n * nu * nu;
   assert(m != 0.);
   n = EA / E;
   lambda = E * nu * (1. + n * nu) / (m * (1. + nu));
   mu = 0.5 * E / (1. + nu);
   alpha = mu - GA;
   beta = E * nu * nu * (1 - n) / (4 * m * (1 + nu));
   gamma = EA * (1. - nu) / (8 * m) - (lambda + 2. * mu) / 8. + 0.5 * alpha - beta;
}

double ShearClosureTransverselyIsotropic::ComputeShearEnergy(const DenseMatrix &F) const
{
   double e_shear = 0.0;
   /* NH part (isotropic) */
   double I1, I2, I3, J;
   DenseMatrix FT(3), C(3);
   FT.Transpose(F);
   mfem::Mult(FT, F, C);
   ComputeIsotropicInvariants(C, I1, I2, I3);
   J = F.Det();
   e_shear = 0.5 * mu * (I1 - 3.) - mu * log(J) + 0.5 * lambda * pow(J - 1., 2);

   /* Transversely isotropic part */
   double I4, I5;
   ComputeAnisotropicInvariants(mi, C, I4, I5);
   e_shear += (alpha + beta * log(J) + gamma * (I4 - 1.) ) * (I4 - 1.) - 0.5 * alpha * (I5 - 1.);

   return e_shear;
}

void ShearClosureTransverselyIsotropic::ComputeCauchyStress(const DenseMatrix &F, DenseMatrix &sigma) const
{
   // Implement Cauchy stress computation for transversely isotropic material
   /* NH part (isotropic) */

   /* Transversely isotropic part */
   MFEM_ABORT("NEED OVERRIDE");
}

} // namespace hydroLO
} // namespace mfem
