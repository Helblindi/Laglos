//NF//MS
/*
 * This is the elastic class, it is used to compute the elastic energy 
 * and the deviatoric part of the stress tensor.
 */
#ifndef ELASTIC
#define ELASTIC

#include "mfem.hpp"
#include "problem_base.h"

using namespace std;

namespace mfem
{

namespace hydroLO
{

enum class ShearEnergyMethod {
   AVERAGE_F,
   AVERAGE_C,
   AVERAGE_ENERGY
};

template<int dim>
class Elastic
{
private:
   ShearEnergyMethod shear_method;
   ParFiniteElementSpace &H1, &L2;
   const ParGridFunction &rho0_gf;
   const IntegrationRule &ir;
   const int NE, NDofs_L2, nqp;

   // Reference to physical Jacobian for the initial mesh.
   // These are computed only at time zero and stored here.
   DenseTensor Jac0inv;
   double mu = -1.; // Shear modulus

public:
   Elastic(ParFiniteElementSpace &h1_fes,
           ParFiniteElementSpace &l2_fes,
           const ParGridFunction &rho0_gf,
           const IntegrationRule &ir,
           ShearEnergyMethod method = ShearEnergyMethod::AVERAGE_C) : 
      H1(h1_fes), 
      L2(l2_fes),
      rho0_gf(rho0_gf),
      ir(ir),
      NE(H1.GetParMesh()->GetNE()),
      NDofs_L2(L2.GetNDofs()), 
      nqp(ir.GetNPoints()),
      Jac0inv(dim, dim, NE * nqp),
      shear_method(method)
   {
      cout << "=== Elastic constructor ===\n";

      if (NDofs_L2 != NE)
      {
         MFEM_ABORT("Error: Number of L2 dofs must be equal to the number of elements.\n");
      }

      /* Construct initial inverse jacobian transformations*/
      Jac0inv = 0.;
      for (int e = 0; e < NE; e++)
      {
         ElementTransformation *T = H1.GetElementTransformation(e);
         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            T->SetIntPoint(&ip);
            DenseMatrixInverse Jinv(T->Jacobian());
            Jinv.GetInverseMatrix(Jac0inv(e*nqp + q));
         }
      }
   }

   void set_shear_modulus(const double &_mu) { this->mu = _mu; }

   void setShearEnergyMethod(ShearEnergyMethod method) { shear_method = method; }

   double e_sheer(const int &e) const
   {
      DenseMatrix c(3), c2(3);
      const double rho0 = rho0_gf(e);
      switch (shear_method) {
         case ShearEnergyMethod::AVERAGE_F:
         {
            Compute_cFromAvgF(e, c);
            mfem::Mult(c, c, c2);
            return e_sheer(c.Trace(), c2.Trace(), rho0);
         }
         case ShearEnergyMethod::AVERAGE_C:
         {
            Compute_cAvg(e,c);
            mfem::Mult(c, c, c2);
            return e_sheer(c.Trace(), c2.Trace(), rho0);
         }
         case ShearEnergyMethod::AVERAGE_ENERGY:
            MFEM_ABORT("Not implemented");
         default:
            MFEM_ABORT("Unknown shear energy method");
      }
   }
   
   double e_sheer(const double &trc, const double &trc2, const double &rho0) const
   {
      if (mu == -1.)
      {
         MFEM_ABORT("Must set shear modulus.\n");
      }
      return mu / 2 * (trc - 3.) / rho0; // Neo hookean
   }

   double des_dtrc(const double &trc, const double &trc2) const
   {
      return mu / 2.; // Neo hookean
   }

   double des_dtrc2(const double &trc, const double &trc2) const
   {
      return 0.; // Neo hookean
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
         mfem::Mult(Jr, Jac0inv(e*nqp + q), Ji);
         F_dim.Add(ip.weight, Ji);
         // TODO: Try this averaging on F F^T instead of F.
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
         mfem::Mult(Jr, Jac0inv(e*nqp + q), F_dim);
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
         mfem::Mult(FT, F, cqp);
         cqp *= std::pow(cqp.Det(), -1./3.);

         /* Add in contribution from quadrature point */
         c.Add(ip.weight, cqp);
      }
      // cout << "c: ";
      // c.Print(cout);
   }

   void ComputeS(const int &e, const double &rho, DenseMatrix &S) const
   {
      /* Objects needed in function */
      DenseMatrix c(3), c2(3), I(3);
      for (int i = 0; i < 3; i++) { I(i,i) = 1.; }

      S.SetSize(3);
      S = 0.;

      /* Compute c */
      Compute_cFromAvgF(e, c);

      /* Compute sheer energy, and save in class member */
      mfem::Mult(c, c, c2);
      const double rho0 = rho0_gf(e);
      
      /* Compute deviatoric part of stess tensor */
      DenseMatrix c_tf(3), c2_tf(3); // 'trace-free objects'
      // remove multiple of trace from diagonal
      mfem::Add(c, I, -1./3. * c.Trace(), c_tf);
      mfem::Add(c2, I, -1./3. * c2.Trace(), c2_tf);

      double c_coeff = des_dtrc(c.Trace(), c2.Trace());
      double c2_coeff = 2.*des_dtrc2(c.Trace(), c2.Trace());
      mfem::Add(c_coeff, c_tf, c2_coeff, c2_tf, S);
      S *= -2. * rho / rho0_gf(e);
   }

};
} // end ns hydroLO

} // end ns mfem

#endif // ELASTIC