#include "shear_closure.hpp"
using namespace std;

namespace mfem {
namespace hydroLO {

void ShearClosure::ComputeIsotropicInvariants(const DenseMatrix &C, double &J1, double &J2, double &J3) const
{
   J1 = C.Trace();
   DenseMatrix C2(3);
   mfem::Mult(C, C, C2);
   // MFEM_WARNING("SHEAR ENERGY CHECK");
   J2 = 0.5 * (J1 * J1 - C2.Trace());
   J3 = C.Det();
}
void ShearClosure::ComputeIsotropicInvariantsReduced(const DenseMatrix &C, double &j1, double &j2) const
{
   DenseMatrix c(3);
   c = C;
   c *= std::pow(C.Det(), -1./3.);
   j1 = c.Trace();
   DenseMatrix c2(3);
   mfem::Mult(c, c, c2);
   j2 = c2.Trace();
}

void ShearClosure::ComputeAnisotropicInvariants(const Vector &mi, const DenseMatrix &C, double &j4, double &j5) const
{
   // MFEM_WARNING("SHEAR ENERGY CHECK");
   Vector Cmi(3), C2mi(3);
   C.Mult(mi, Cmi);
   j4 = mfem::InnerProduct(mi, Cmi);

   DenseMatrix C2(3);
   mfem::Mult(C, C, C2);
   C2.Mult(mi, C2mi);
   j5 = mfem::InnerProduct(mi, C2mi);
}

void ShearClosure::ComputeCauchyStress(const DenseMatrix &F, DenseMatrix &sigmaD) const
{
   DenseMatrix FT(3), C(3), c(3), c2(3), B(3), b(3),b2(3), I(3);
   for (int i = 0; i < 3; i++) { I(i,i) = 1.; }
   sigmaD.SetSize(3);
   sigmaD = 0.;

   FT.Transpose(F);
   mfem::Mult(FT, F, C);
   c = C;
   c *= std::pow(C.Det(), -1./3.);
   mfem::Mult(c, c, c2);
   mfem::Mult(F, FT, B);
   b = B;
   b *= std::pow(B.Det(), -1./3.);
   mfem::Mult(b, b, b2);
   
   /* Compute deviatoric part of stess tensor */
   DenseMatrix b_tf(3), b2_tf(3); // 'trace-free objects'
   // remove multiple of trace from diagonal
   mfem::Add(b, I, -1./3. * b.Trace(), b_tf);
   mfem::Add(b2, I, -1./3. * b2.Trace(), b2_tf);

   double des_dj1, des_dj2;
   this->ComputeShearEnergyIsotropicDerivatives(C, des_dj1, des_dj2);
   double c_coeff = des_dj1;
   double c2_coeff = 2.*des_dj2;
   mfem::Add(c_coeff, b_tf, c2_coeff, b2_tf, sigmaD);
   if (this->is_anisotropic)
   {
      DenseMatrix sigmaD_an(3);
      this->ComputeCauchyStressAnisotropicComponent(F, sigmaD_an);
      sigmaD.Add(1., sigmaD_an);
   }

   sigmaD *= 2. / F.Det();
}




/**************************************************************************
* Neo Hookean Shear Closure Class
***************************************************************************/
double ShearClosureNeoHookean::ComputeShearEnergy(const DenseMatrix &F, const double &rho0) const
{
   // cout << "\n\nShearClosureNeoHookean::ComputeShearEnergy\n\n";
   double j1, j2;
   DenseMatrix FT(3), C(3);
   FT.Transpose(F);
   mfem::Mult(FT, F, C);
   ComputeIsotropicInvariantsReduced(C, j1, j2);
   return mu / 2. / rho0 * (j1 - 3.);
}

void ShearClosureNeoHookean::ComputeShearEnergyIsotropicDerivatives(const DenseMatrix &C, double &des_dj1, double &des_dj2) const
{
   des_dj1 = mu / 2.;
   des_dj2 = 0.;
}




/**************************************************************************
* Mooney Rivlin Shear Closure Class
***************************************************************************/
double ShearClosureMooneyRivlin::ComputeShearEnergy(const DenseMatrix &F, const double &rho0) const
{
   // cout << "\n\nShearClosureMooneyRivlin::ComputeShearEnergy\n\n";
   double j1, j2;
   DenseMatrix FT(3), C(3);
   FT.Transpose(F);
   mfem::Mult(FT, F, C);
   ComputeIsotropicInvariantsReduced(C, j1, j2);
   return mu / 32. / rho0 * (pow(j1,4) - 2*j2*j1*j1 + j2*j2 - 8.*j1 - 12.);
}

void ShearClosureMooneyRivlin::ComputeShearEnergyIsotropicDerivatives(const DenseMatrix &C, double &des_dj1, double &des_dj2) const
{
   double j1, j2;
   ComputeIsotropicInvariantsReduced(C, j1, j2);
   des_dj1 = mu / 8. * (pow(j1,3) - j1*j2 - 2.);
   des_dj2 = mu / 16. * (-j1*j1 + j2);
}


/**************************************************************************
* Aortic Shear Closure Class
***************************************************************************/
ShearClosureAortic::ShearClosureAortic(const double &_mu, const Vector &_mi, const double &_w1, const double &_D1,
                     const double &_A1, const double &_B1)
      : ShearClosure(_mu, _mi), Gi(3), w1(_w1), D1(_D1), A1(_A1), B1(_B1) 
{
   /* Form Gi */
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         Gi.Elem(i,j) = mi[i] * mi[j];
      }
   }
}
void ShearClosureAortic::ComputeAnisotropicInvariantsAortic(const DenseMatrix &C, double &j4, double &j5) const
{
   /* compute j4 and j5 */
   DenseMatrix c(3);
   c = C;
   c *= std::pow(C.Det(), -1./3.);
   DenseMatrix _res(3), _res2(3);
   Mult(c, Gi, _res);
   j4 = _res.Trace();
   Mult(c, c, _res);
   Mult(_res, Gi, _res2);
   j5 = _res2.Trace();
}
double ShearClosureAortic::ComputeShearEnergy(const DenseMatrix &F, const double &rho0) const
{
   // cout << "\n\nShearClosureAortic::ComputeShearEnergy\n\n";
   double j1, j2;
   double j4, j5;
   DenseMatrix FT(3), C(3);
   FT.Transpose(F);
   mfem::Mult(FT, F, C);
   ComputeIsotropicInvariantsReduced(C, j1, j2);
   ComputeAnisotropicInvariantsAortic(C, j4, j5);
   
   /* Compute shear energy */
   double val = (1. - 2 * w1) * ((A1 / 2.) * (j1 - 3.) + (B1 / 4.) * (j1*j1 - j2 - 6.));

   double _tt1 = (1. / A2) * (exp(A2 * (j4 - 1.)) - 1.);
   double _tt2 = (1. / B2) * (exp(B2 * (j5 - j1 * j4 + (j1*j1 - j2) / 2. - 1.)) - 1.);
   double _anisotropic_val = 2. * w1 * D1 * (_tt1 + _tt2);
   // if (val > 1.E-10 || _anisotropic_val > 1.E-10)
   // {
   //    cout << "isotropic: " << val << ", anisotropic: " << _anisotropic_val << endl;
   // }

   
   if (!std::isnan(_anisotropic_val))
   {
      val += _anisotropic_val;
   }
   else
   {
      // cout << "_anisotropic_val is nan\n";
      // cout << "w1: " << w1 << ", D1: " << D1 << ", A1: " << A1 << ", B1: " << B1 << endl;
      // cout << "trc: " << trc << ", trc2: " << trc2 << endl;
      // cout << "i4: " << i4 << ", i5: " << i5 << endl;
      // cout << "tt1: " << _tt1 << ", tt2: " << _tt2 << endl;
      MFEM_WARNING("Aortic EOS anisotropic is NaN.\n");
   }

   // cout << "anisotropic ratio: " << _anisotropic_val / val << endl;

   return val;
}

void ShearClosureAortic::ComputeShearEnergyIsotropicDerivatives(const DenseMatrix &C, double &des_dj1, double &des_dj2) const
{
   /* Compute trace values */
   double j1, j2, j4, j5;
   ComputeIsotropicInvariantsReduced(C, j1, j2);
   ComputeAnisotropicInvariantsAortic(C, j4, j5);

   des_dj1 = (1. - 2. * w1) * (A1 / 2. + (B1 / 2.) * j1);
   des_dj1 += 2. * w1 * D1 * exp(B2 * (j5 - j1 * j4 + (j1*j1 - j2) / 2. - 1.)) * (-j4 + j1);
   
   des_dj2 = (2. * w1 - 1.) * (B1 / 4.);
   des_dj2 -= w1 * D1 * B2 / B1 * exp(B2 * (j5 - j1 * j4 + (j1*j1 - j2) / 2. - 1.)) * (-j4 + j1);
}

void ShearClosureAortic::ComputeShearEnergyAnisotropicDerivatives(const DenseMatrix &C, double &des_dj4, double &des_dj5) const
{
   /* Compute trace values */
   double j1, j2, j4, j5;
   ComputeIsotropicInvariantsReduced(C, j1, j2);
   ComputeAnisotropicInvariantsAortic(C, j4, j5);

   des_dj4 = 2. * w1 * D1 * (exp(A2 * (j4 - 1.)) - j1 * exp(B2 * (j5 - j1 * j4 + (j1*j1 - j2) / 2. - 1.)));
   des_dj5 = 2. * w1 * D1 * B2 * exp(B2 * (j5 - j1 * j4 + (j1*j1 - j2) / 2. - 1.));
}

void ShearClosureAortic::ComputeCauchyStressAnisotropicComponent(const DenseMatrix &F, DenseMatrix &sigmaD_an) const
{
   DenseMatrix FT(3), C(3);
   FT.Transpose(F);
   mfem::Mult(FT, F, C);

   double des_dj4, des_dj5;
   ComputeAnisotropicInvariantsAortic(C, des_dj4, des_dj5);

   DenseMatrix c(3), c2(3), I(3);
   c = C;
   c *= std::pow(C.Det(), -1./3.);
   mfem::Mult(c, c, c2);

   /* des_dj4 portion */
   double c4_coeff = des_dj4;
   DenseMatrix c4_tf(3), _tmp(3);
   mfem::Mult(F, Gi, _tmp);
   mfem::Mult(_tmp, FT, c4_tf);
   c4_tf *= std::pow(C.Det(), -1./3.);
   mfem::Mult(c, Gi, _tmp);
   mfem::Add(c4_tf, I, -1./3. * _tmp.Trace(), c4_tf);
   mfem::Add(sigmaD_an, c4_tf, c4_coeff, sigmaD_an);

   /* des_dj5 portion */
   double c5_coeff = 2. * des_dj5;
   DenseMatrix c5_tf(3), _tmp2(3);
   mfem::Mult(F, C, _tmp);
   mfem::Mult(_tmp, Gi, _tmp2);
   mfem::Mult(_tmp2, F, c5_tf);
   c5_tf *= std::pow(C.Det(), -2./3.);
   mfem::Mult(c2, Gi, _tmp);
   mfem::Add(c5_tf, I, -1./3. * _tmp.Trace(), c5_tf);
   mfem::Add(sigmaD_an, c5_tf, c5_coeff, sigmaD_an);
}


/**************************************************************************
* Transversely Isotropic Shear Closure Class
***************************************************************************/
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

double ShearClosureTransverselyIsotropic::ComputeShearEnergy(const DenseMatrix &F, const double &rho0) const
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

void ShearClosureTransverselyIsotropic::ComputeShearEnergyIsotropicDerivatives(const DenseMatrix &C, double &des_dj1, double &des_dj2) const
{
   // Implement isotropic derivative computations
   /* NH part (isotropic) */
   MFEM_ABORT("NEED OVERRIDE");
}

void ShearClosureTransverselyIsotropic::ComputeShearEnergyAnisotropicDerivatives(const DenseMatrix &C, double &des_dj4, double &des_dj5) const
{
   // Implement anisotropic derivative computations
   /* NH part (isotropic) */
   MFEM_ABORT("NEED OVERRIDE");
}

void ShearClosureTransverselyIsotropic::ComputeCauchyStressAnisotropicComponent(const DenseMatrix &F, DenseMatrix &sigmaD_an) const
{
   MFEM_ABORT("NEED OVERRIDE FOR ANISOTROPIC");
}


} // namespace hydroLO
} // namespace mfem
