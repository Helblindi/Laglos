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

template<int dim>
class Elastic
{
private:
   ParFiniteElementSpace &H1, &L2;
   const ParGridFunction &rho0_gf;
   const IntegrationRule &ir;
   const int NE, NDofs_L2, nqp;

   // Reference to physical Jacobian for the initial mesh.
   // These are computed only at time zero and stored here.
   DenseTensor Jac0inv;

public:
   Elastic(ParFiniteElementSpace &h1_fes,
           ParFiniteElementSpace &l2_fes,
           const ParGridFunction &rho0_gf,
           const IntegrationRule &ir) : 
      H1(h1_fes), 
      L2(l2_fes),
      rho0_gf(rho0_gf),
      ir(ir),
      NE(H1.GetParMesh()->GetNE()),
      NDofs_L2(L2.GetNDofs()), 
      nqp(ir.GetNPoints()),
      Jac0inv(dim, dim, NE * nqp)
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

   double e_sheer(const int &e) const
   {
      DenseMatrix c(3), c2(3);
      Compute_c(e, c);

      mfem::Mult(c, c, c2);
      return e_sheer(c.Trace(), c2.Trace());
   }
   
   double e_sheer(const double &trc, const double &trc2) const
   {
      double _mu = 1.;
      return _mu * (trc/3. - 1.); // Neo hookean
   }

   double des_dtrc(const double &trc, const double &trc2) const
   {
      double _mu = 1.;
      return _mu / 3.; // Neo hookean
   }

   double des_dtrc2(const double &trc, const double &trc2) const
   {
      double _mu = 1.;
      return 0.; // Neo hookean
   }

   void ComputeF(const int e, DenseMatrix &F) const
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

   void Compute_c(const int &e, DenseMatrix &c) const
   {
      DenseMatrix F(3), FT(3), C(3);
      c.SetSize(3);
      ComputeF(e,F);
      FT.Transpose(F);
      mfem::Mult(FT, F, C);
      c = C;
      c *= std::pow(C.Det(), -1./3.);
   }

   void ComputeS(const int &e, const double &rho, double &es, DenseMatrix &S) const
   {
      /* Objects needed in function */
      DenseMatrix c(3), c2(3), I(3);
      for (int i = 0; i < 3; i++) { I(i,i) = 1.; }

      S.SetSize(3);
      S = 0.;

      /* Compute c */
      Compute_c(e, c);

      /* Compute sheer energy, and save in class member */
      mfem::Mult(c, c, c2);
      es = e_sheer(c.Trace(), c2.Trace());
      
      /* Compute deviatoric part of stess tensor */
      DenseMatrix c_tf(3), c2_tf(3); // 'trace-free objects'
      // remove multiple of trace from diagonal
      mfem::Add(c, I, -1./3. * c.Trace(), c_tf);
      mfem::Add(c2, I, -1./3. * c2.Trace(), c2_tf);

      double c_coeff = des_dtrc(c.Trace(), c2.Trace());
      double c2_coeff = 2.*des_dtrc2(c.Trace(), c2.Trace());
      mfem::Add(c_coeff, c_tf, c2_coeff, c2_tf, S);
      S *= 2. * rho / rho0_gf(e);
   }

};
} // end ns hydroLO

} // end ns mfem

#endif // ELASTIC