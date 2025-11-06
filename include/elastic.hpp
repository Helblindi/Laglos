//NF//MS
/*
* This is the elastic class, it is used to compute the elastic energy 
* and the deviatoric part of the stress tensor.
*/
#ifndef ELASTIC
#define ELASTIC

#include "mfem.hpp"
#include "problem_base.h"
#include "shear_closure.hpp"
#include "laglos_assembly.hpp"

using namespace std;

namespace mfem
{

namespace hydroLO
{

enum ShearEOS {
   NEO_HOOKEAN,
   MOONEY_RIVLIN,
   AORTIC,         // anisotropic model
   TRANSVERSELY_ISOTROPIC
};

enum ShearEnergyMethod {
   AVERAGE_F,
   AVERAGE_C,
   AVERAGE_ENERGY
};
 
class Elastic
{
private:
   const int dim;
   ShearEnergyMethod shear_method;
   ShearEOS shear_eos;
   ParFiniteElementSpace &H1;
   const ParGridFunction &rho0_gf;
   const IntegrationRule &ir;
   const int NE, nqp;

   ShearClosure *shear_closure_model = nullptr;

   // Reference to physical Jacobian for the initial mesh
   // These are computed only at time zero and stored here
   const QuadratureData &quad_data;
   double mu = -1.; // Shear modulus

   /* 
   If using the aortic eos, a fiber direction must be defined. 
   This direction is used to compute the invariants:
      m = (cos(theta), sin(theta), 0) 
   */
   double theta = M_PI/4.; // angle of fiber direction
   Vector mi_vec;
   DenseMatrix Gi;
   
   // double A1 = 21.5802, B1 = 9.9007, D1 = 0.8849, w1 = 0.4189; /* PP too squished */
   // double A1 = 771.8033, B1 = 21.2093, D1 = 3.8068, w1 = 0.4971;
   // double A2 = 1., B2 = 1.;

   /* These params were closer to the Neo Hook results when w1 = 0 */
   double stiffness = 9.63E5;
   double A1 = 0.5 * stiffness, B1 = 0.5 * stiffness, A2 = 1., B2 = 1., D1 = 0.5 * (1.5*stiffness), w1 = 0.49;

public:
   Elastic(const int &_dim,
           const int &_elastic_eos,
           QuadratureData &_quad_data,
           ParFiniteElementSpace &h1_fes,
           const ParGridFunction &rho0_gf,
           const IntegrationRule &ir,
           ShearEnergyMethod method = ShearEnergyMethod::AVERAGE_F) : 
      dim(_dim),
      quad_data(_quad_data),
      H1(h1_fes), 
      rho0_gf(rho0_gf),
      ir(ir),
      NE(H1.GetParMesh()->GetNE()),
      nqp(ir.GetNPoints()),
      shear_method(method)
   {
      /* Set equation of state */
      switch(_elastic_eos)
      {
         case 1: // Neo-Hookean
            shear_eos = ShearEOS::NEO_HOOKEAN;
            shear_closure_model = new ShearClosureNeoHookean(mu);
            break;
         case 2: // Mooney-Rivlin
            shear_eos = ShearEOS::MOONEY_RIVLIN;
            shear_closure_model = new ShearClosureMooneyRivlin(mu);
            break;
         case 3: // Aortic
         {
            /* Invariants */
            mi_vec.SetSize(3);
            Gi.SetSize(3);
            /* Set fiber direction */
            mi_vec(0) = cos(theta);
            mi_vec(1) = sin(theta);
            mi_vec(2) = 0.;
            tensor(mi_vec, mi_vec, Gi);
            shear_eos = ShearEOS::AORTIC;
            shear_closure_model = new ShearClosureAortic(mu, mi_vec, w1, D1, A1, B1);
            break;
         }
         case 4: // Transversely Isotropic
         {
            shear_eos = ShearEOS::TRANSVERSELY_ISOTROPIC;
            double E = 1.729E6, EA = 1.E6, GA = 4.59E5, nu = 0.4;
            mi_vec.SetSize(3);
            /* Set fiber direction */
            mi_vec(0) = cos(theta);
            mi_vec(1) = sin(theta);
            mi_vec(2) = 0.;
            shear_closure_model = new ShearClosureTransverselyIsotropic(mu, mi_vec, E, EA, GA, nu);
            
            break;
         }
         default:
            MFEM_ABORT("Invalid value for shear_eos.");
      }

      /*************** CONFIGURATION OUTPUT ***************/
      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
      << "@         Elastic Class Configuration      @\n"
      << "@------------------------------------------@\n"
      << "@ shear_eos         : " << std::setw(20) << std::left;
      switch (shear_eos)
      {
         case NEO_HOOKEAN:
            cout << "NEO_HOOKEAN";
            break;
         case MOONEY_RIVLIN:
            cout << "MOONEY_RIVLIN";
            break;
         case AORTIC:
            cout << "AORTIC";
            break;
         case TRANSVERSELY_ISOTROPIC:
            cout << "TRANSVERSELY_ISOTROPIC";
            break;
         default:
            MFEM_ABORT("Invalid value for shear_eos.");
      }
      cout << " @\n";
      if (shear_eos == ShearEOS::AORTIC)
      {
         cout << "@ -- Aortic parameters -- " << std::setw(18) << std::right << "@" << "\n"
              << "@ \ttheta         : " << std::setw(21) << std::left << theta << "@" << "\n"
              << "@ \tm(0)        : " << std::setw(21) << std::left << mi_vec(0) << "@" << "\n"
              << "@ \tm(1)        : " << std::setw(21) << std::left << mi_vec(1) << "@" << "\n"
              << "@ \tm(2)        : " << std::setw(21) << std::left << mi_vec(2) << "@" << "\n";
      }

      
      cout << "@ ShearEnergyMethod : " << std::setw(20) << std::left;
      switch (shear_method)
      {
         case ShearEnergyMethod::AVERAGE_F:
            cout << "AVERAGE_F";
            break;
         case ShearEnergyMethod::AVERAGE_C:
            cout << "AVERAGE_C";
            break;
         default:
            MFEM_ABORT("Unknown shear energy method");
      }
      cout << " @\n";
      
      cout << "@ nqp               : " << std::setw(20) << std::left << nqp      << " @\n"
           << "@ NE                : " << std::setw(20) << std::left << NE       << " @\n"
           << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n";
      /***************** END CONFIGURATION *****************/
   }
 
   void tensor(const Vector & v1, const Vector & v2, DenseMatrix & dm) const
   {
      const int v1_len = v1.Size(), v2_len = v2.Size();
      for (int i = 0; i < v1_len; i++)
      {
         for (int j = 0; j < v2_len; j++)
         {
            dm.Elem(i,j) = v1[i] * v2[j];
         }
      }
   }
   void set_shear_modulus(const double &_mu) { 
      this->mu = _mu; 
      shear_closure_model->UpdateShearModulus(_mu);
   }

   void setShearEnergyMethod(ShearEnergyMethod method) { 
      shear_method = method; 
      if (shear_method != ShearEnergyMethod::AVERAGE_F && shear_eos == ShearEOS::AORTIC)
      {
         MFEM_ABORT("Aortic EOS only supports shear energy method AVERAGE_F.\n");
      }
   }

   /* Wrapper function that needs to be moved to ElasticProblemBase class */
   double e_sheer(const int &e) const
   {
      const double rho0 = rho0_gf(e);

      DenseMatrix F(dim);
      ComputeAvgF(e, F);
      double _e_shear = shear_closure_model->ComputeShearEnergy(F, rho0);
      return _e_shear;
   }   

   void compute_i4_i5(const DenseMatrix &c, double &i4, double &i5) const
   {
      assert(shear_eos == ShearEOS::AORTIC);

      /* compute i4 and i5 */
      DenseMatrix _res(3), _res2(3);
      Mult(c, Gi, _res);
      i4 = _res.Trace();
      Mult(c, c, _res);
      Mult(_res, Gi, _res2);
      i5 = _res2.Trace();
   }
 
   double des_di4(const DenseMatrix &c) const
   {
      /* Compute trace values */
      const double trc = c.Trace();
      DenseMatrix c2(3);
      mfem::Mult(c,c,c2);
      const double trc2 = c2.Trace();

      assert(shear_eos == ShearEOS::AORTIC);

      if (mu == -1.)
      {
         MFEM_ABORT("Must set shear modulus.\n");
      }

      /* compute i4 and i5 */
      double i4, i5;
      compute_i4_i5(c, i4, i5);

      double val = 2. * w1 * D1 * (exp(A2 * (i4 - 1.)) - trc * exp(B2 * (i5 - trc * i4 + (trc*trc - trc2) / 2. - 1.)));
      if (std::isnan(val))
      {
         val = 0.;
         MFEM_WARNING("Aortic EOS shear energy derivative computation failed in des_di4.\n");
      }
      return val;
   }
 
   double des_di5(const DenseMatrix &c) const
   {
      /* Compute trace values */
      const double trc = c.Trace();
      DenseMatrix c2(3);
      mfem::Mult(c,c,c2);
      const double trc2 = c2.Trace();

      assert(shear_eos == ShearEOS::AORTIC);

      if (mu == -1.)
      {
         MFEM_ABORT("Must set shear modulus.\n");
      }

      /* compute i4 and i5 */
      double i4, i5;
      compute_i4_i5(c, i4, i5);

      double val = 2. * w1 * D1 * exp(B2 * (i5 - trc * i4 + (trc*trc - trc2) / 2. - 1.));
      if (std::isnan(val))
      {
         val = 0.;
         MFEM_WARNING("Aortic EOS shear energy derivative computation failed in des_di5.\n");
      }
      return val;
   }
 
   void ComputeAvgF(const int e, DenseMatrix &F) const
   {
      DenseMatrix F_dim(dim); F_dim = 0.;
      ElementTransformation *T = H1.GetElementTransformation(e);
      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         T->SetIntPoint(&ip);
         DenseMatrix Jr = T->Jacobian();
         /* Instead of a mapping from the reference element,
            we need a mapping from the initial configuration
            of the given element */
         DenseMatrix Ji(dim);
         mfem::Mult(Jr, quad_data.Jac0inv(e*nqp + q), Ji);
         F_dim.Add(ip.weight, Ji);
      }
      /* If dim != 3, will need to augment F with 1s on the diagonal */
      if (dim == 3)
      {
         F = F_dim;
      }
      else {
         /* Set submatrix to F_dim */
         F.SetSize(3);
         F = 0.;
         Array<int> idx(dim);
         for (int i = 0; i < dim; i++)
         {
            idx[i] = i;
         }
         F.SetSubMatrix(idx, idx, F_dim);

         /* Fill remaining entries with a 1 */
         for (int i = dim; i < 3; i++)
         {
            F(i,i) = 1.;
         }
      }
   }
 
   void Compute_cFromAvgF(const int &e, DenseMatrix &c) const
   {
      DenseMatrix F(3), FT(3), C(3);
      c.SetSize(3);
      ComputeAvgF(e,F);
      FT.Transpose(F);
      mfem::Mult(FT, F, C);
      c = C;
      c *= std::pow(C.Det(), -1./3.);
   }
 
   void Compute_cAvg(const int &e, DenseMatrix &c) const
   {
      // cout << "!!!!!!!!!!!!!!!!!Elastic::Compute_cAvg!!!!!!!!!!!!!!!!!\n";
      DenseMatrix F_dim(dim);
      DenseMatrix F(3), FT(3), FTF(3);
      DenseMatrix C(3), cqp(3);
      ElementTransformation *T = H1.GetElementTransformation(e);

      c.SetSize(3);
      c = 0.;
      for (int q = 0; q < nqp; q++)
      {
         F_dim = 0.;
         F = 0.; FT = 0.; FTF = 0.; C = 0.;
         const IntegrationPoint &ip = ir.IntPoint(q);
         T->SetIntPoint(&ip);
         DenseMatrix Jr = T->Jacobian();
         /* Instead of a mapping from the reference element,
            we need a mapping from the initial configuration
            of the given element */
         DenseMatrix Ji(dim);
         mfem::Mult(Jr, quad_data.Jac0inv(e*nqp + q), F_dim);
         if (dim == 3)
         {
            F = F_dim;
         }
         else { /* If dim != 3, will need to augment F with 1s on the diagonal */
            /* Set submatrix to F_dim */
            Array<int> idx(dim);
            for (int i = 0; i < dim; i++)
            {
               idx[i] = i;
            }
            F.SetSubMatrix(idx, idx, F_dim);

            /* Fill remaining entries with a 1 */
            for (int i = dim; i < 3; i++)
            {
               F(i,i) = 1.;
            }
         }

         /* Compute FTF at quadrature point, and normalize */
         FT.Transpose(F);
         mfem::Mult(F, FT, cqp);
         cqp *= std::pow(cqp.Det(), -1./3.);

         /* Add in contribution from quadrature point */
         c.Add(ip.weight, cqp);
      }
   }
 
   void ComputeS(const int &e, DenseMatrix &S) const
   {
      /* Objects needed in function */
      DenseMatrix F(3), FT(3), B(3), b(3), b2(3), C(3), c(3), c2(3), I(3);
      for (int i = 0; i < 3; i++) { I(i,i) = 1.; }
      const double rho0 = rho0_gf(e);

      S.SetSize(3);
      S = 0.;

      /* Compute F, B, b, C, and c */
      assert(shear_method == ShearEnergyMethod::AVERAGE_F);
      ComputeAvgF(e, F);
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
      shear_closure_model->ComputeShearEnergyIsotropicDerivatives(C, des_dj1, des_dj2);
      double c_coeff = des_dj1;
      double c2_coeff = 2.*des_dj2;
      mfem::Add(c_coeff, b_tf, c2_coeff, b2_tf, S);

      /* Handle anistropic contribution */
      if (shear_eos == ShearEOS::AORTIC)
      {
         // cout << "===== Anisotropic contribution =====\n";
         /* des_di4 portion */
         double c4_coeff = des_di4(c);
         DenseMatrix c4_tf(3), c5_tf(3), _tmp(3), _tmp2(3);
         mfem::Mult(F, Gi, _tmp);
         mfem::Mult(_tmp, FT, c4_tf);
         c4_tf *= std::pow(C.Det(), -1./3.);
         mfem::Mult(c, Gi, _tmp);
         mfem::Add(c4_tf, I, -1./3. * _tmp.Trace(), c4_tf);
         mfem::Add(S, c4_tf, c4_coeff, S);

         /* des_di5 portion */
         double c5_coeff = 2. * des_di5(c);
         mfem::Mult(F, C, _tmp);
         mfem::Mult(_tmp, Gi, _tmp2);
         mfem::Mult(_tmp2, F, c5_tf);
         c5_tf *= std::pow(C.Det(), -2./3.);
         mfem::Mult(c2, Gi, _tmp);
         mfem::Add(c5_tf, I, -1./3. * _tmp.Trace(), c5_tf);
         mfem::Add(S, c5_tf, c5_coeff, S);
         // cout << "=====================================\n";
      }

      /* Finally, multiply by 2rho */
      S *= 2. / F.Det();
   }
}; // End Elastic class
} // end ns hydroLO

} // end ns mfem

#endif // ELASTIC
