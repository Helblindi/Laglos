#include "laglos_assembly.hpp"

namespace mfem
{

namespace hydroLO
{

void DGNormalIntegrator::AssembleFaceMatrix(const FiniteElement &el1,
                                            const FiniteElement &el2,
                                            FaceElementTransformations &Trans,
                                            DenseMatrix &elmat)
{
   // cout << "DGNormalIntegrator::AssembleFaceMatrix\n";
   int nd1 = el1.GetDof(), nd2;
   int dim = el1.GetDim();
   Vector nor(dim);
   real_t w;

   if (Trans.Elem2No >= 0)
   {
      nd2 = el2.GetDof();
   }
   else
   {
      nd2 = 0;
   }

   shape1.SetSize(nd1);
   shape2.SetSize(nd2);
   elmat.SetSize(nd1 + nd2);
   elmat = 0.0;

   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      int order;
      // Assuming order(u)==order(mesh)
      if (Trans.Elem2No >= 0)
         order = (min(Trans.Elem1->OrderW(), Trans.Elem2->OrderW()) +
                  2*max(el1.GetOrder(), el2.GetOrder()));
      else
      {
         order = Trans.Elem1->OrderW() + 2*el1.GetOrder();
      }
      if (el1.Space() == FunctionSpace::Pk)
      {
         order++;
      }
      ir = &IntRules.Get(Trans.GetGeometryType(), order);
   }

   for (int p = 0; p < ir->GetNPoints(); p++)
   {
      const IntegrationPoint &ip = ir->IntPoint(p);

      // Set the integration point in the face and the neighboring elements
      Trans.SetAllIntPoints(&ip);

      // Access the neighboring elements' integration points
      // Note: eip2 will only contain valid data if Elem2 exists
      const IntegrationPoint &eip1 = Trans.GetElement1IntPoint();
      const IntegrationPoint &eip2 = Trans.GetElement2IntPoint();

      el1.CalcPhysShape(*Trans.Elem1, shape1);

      if (dim == 1)
      {
         nor(0) = 2*eip1.x - 1.0;
      }
      else
      {
         CalcOrtho(Trans.Jacobian(), nor);
      }

      if (xi > dim)
      {
         MFEM_ABORT("Component of normal should not exceed dimension of problem.\n");
      }
      w = 0.5 * coeff * ip.weight * nor[xi];

      if (w != 0.0)
      {
         for (int i = 0; i < nd1; i++)
         {
            for (int j = i; j < nd1; j++)
            {
               // symmetric, see (4.17) 
               // must accomodate '-'
               elmat(i, j) += w * shape1(i) * shape1(j);
               elmat(j, i) += w * shape1(i) * shape1(j);
            }
         }

         if (nd2)
         {
            el2.CalcPhysShape(*Trans.Elem2, shape2);

            for (int i = 0; i < nd2; i++)
               for (int j = 0; j < nd1; j++)
               {
                  elmat(nd1+i, j) += w * shape2(i) * shape1(j);
               }

            for (int i = 0; i < nd2; i++)
               for (int j = i; j < nd2; j++)
               {
                  elmat(nd1+i, nd1+j) += w * shape2(i) * shape2(j);
                  elmat(nd1+j, nd1+i) -= w * shape2(i) * shape2(j);
               }

            for (int i = 0; i < nd1; i++)
               for (int j = 0; j < nd2; j++)
               {
                  elmat(i, nd1+j) -= w * shape1(i) * shape2(j);
               }
         }
      }      

      // /* Set diagonal entries accordingly */
      // Vector diag;
      // elmat.GetDiag(diag);
      // if (w != 0.0 && diag.Norml2() < 1E-12)
      // {
      //    Vector row_sums;
      //    elmat.GetRowSums(row_sums);
      //    for (int i = 0; i < elmat.Height(); i++)
      //    {
      //       elmat(i,i) = -1. * row_sums[i];
      //    }
      // }

   }


}
} // end ns hydrodynamics

} // end ns mfem