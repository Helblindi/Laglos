//NF//MS
/*
* This is the elastic class, it is used to compute the elastic energy 
* and the deviatoric part of the stress tensor.
*/
#ifndef ELASTIC
#define ELASTIC

#include "mfem.hpp"
#include "laglos_assembly.hpp"

using namespace std;

namespace mfem
{

namespace hydroLO
{

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
   ParFiniteElementSpace &H1;
   const IntegrationRule &ir;
   const int NE, nqp;

   // Reference to physical Jacobian for the initial mesh
   // These are computed only at time zero and stored here
   const QuadratureData &quad_data;   

public:
   Elastic(const int &_dim,
           QuadratureData &_quad_data,
           ParFiniteElementSpace &h1_fes,
           const IntegrationRule &ir,
           ShearEnergyMethod method = ShearEnergyMethod::AVERAGE_F) : 
      dim(_dim),
      quad_data(_quad_data),
      H1(h1_fes),
      ir(ir),
      NE(H1.GetParMesh()->GetNE()),
      nqp(ir.GetNPoints()),
      shear_method(method)
   {
      /*************** CONFIGURATION OUTPUT ***************/
      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
      << "@         Elastic Class Configuration      @\n"
      << "@------------------------------------------@\n"
      << "@ ShearEnergyMethod : " << std::setw(20) << std::left;
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

   void setShearEnergyMethod(ShearEnergyMethod method) { shear_method = method; }

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
}; // End Elastic class
} // end ns hydroLO

} // end ns mfem

#endif // ELASTIC
