// IDPLimiter is a class that carries out the limiting technique described in
//
// Jean-Luc Guermond, Zuodong Wang,
// Mass conservative limiting and applications to the approximation of the 
// steady-state radiation transport equations,
// Journal of Computational Physics,
// Volume 521, Part 1,
// 2025,
// 113531,
// ISSN 0021-9991,
// https://doi.org/10.1016/j.jcp.2024.113531.
// (https://www.sciencedirect.com/science/article/pii/S0021999124007794)


#ifndef LIMITER
#define LIMITER

#include "mfem.hpp"

using namespace mfem;
using namespace std;

namespace mfem
{

class IdentityInterpolatorMax : public IdentityInterpolator
{
   virtual void AssembleElementMatrix2(const FiniteElement &dom_fe,
      const FiniteElement &ran_fe,
      ElementTransformation &Trans,
      DenseMatrix &elmat) 
{ 
   ran_fe.Project(dom_fe, Trans, elmat); 
   Vector row;
   /* Threshold by row */
   for (int row_it = 0; row_it < elmat.NumRows(); row_it++)
   {
      elmat.GetRow(row_it, row);
      double _threshold = row.Max();
      if (_threshold < 0.)
      {
         MFEM_ABORT("Negative numbers not handled properly.\n");
      }

      for (int i = 0; i < row.Size(); i++)
      {
         if (std::abs(row.Elem(i)) < _threshold)
         {
            row.Elem(i) = 0.0;
         }
         else
         {
            row.Elem(i) = 1.;
         }
      }
      elmat.SetRow(row_it, row);
   }
}
};


class IDPLimiter
{
private:
   const ParFiniteElementSpace &pfes;
   const HypreParVector &mass_vec;
   const int NDofs;
   mutable ParGridFunction rho_min, rho_max;

   ParFiniteElementSpace H1_LO, H1_HO;
   ParGridFunction LO_proj_max, LO_proj_min;
   ParDiscreteLinearOperator vert_elem_oper; 
   HypreParMatrix *vert_elem;

public:
IDPLimiter(ParFiniteElementSpace & _pfes, ParFiniteElementSpace & H1FESpace_proj_LO,
           ParFiniteElementSpace & H1FESpace_proj_HO, const HypreParVector &_mass_vec) :
   pfes(_pfes),
   mass_vec(_mass_vec),
   NDofs(mass_vec.Size()),
   rho_min(&_pfes),
   rho_max(&_pfes),
   H1_LO(H1FESpace_proj_LO),
   H1_HO(H1FESpace_proj_HO),
   LO_proj_max(&H1_LO),
   LO_proj_min(&H1_LO),
   vert_elem_oper(&H1_HO, &_pfes)
{ 
   /* Verify that #DOFS are equal across LO and HO H1 spaces */
   if (H1_HO.GetNDofs() != H1_LO.GetNDofs())
   {
      MFEM_ABORT("Continuous approximation spaces are not compatible.\n");
   }

   // Try @ https://github.com/mfem/mfem/issues/436
   vert_elem_oper.AddDomainInterpolator(new IdentityInterpolatorMax);
   vert_elem_oper.Assemble();
   vert_elem_oper.Finalize();
   vert_elem = vert_elem_oper.ParallelAssemble();
}

~IDPLimiter() 
{
   // delete spm;
   delete vert_elem;
}

/**
 * @brief Given a low order and high order representation of a thermodynamical variables,
 * perform local mass conservative limiting on the high order solution using the low 
 * order IDP update.
 *
 * @param gf_lo the invariant-domain-preserving low order update
 * @param gf_ho the high order update, not yet IDP
 */
void LocalConservativeLimit(const ParGridFunction &gf_lo, ParGridFunction &gf_ho)
{
   cout << "IDPLimiter::LocalConservativeLimit\n";
   ComputeRhoMinMax(gf_lo);

   MFEM_WARNING("Possible iterations.\n");
   LimitMassMax(gf_ho);
   LimitMassMin(gf_ho);
}

/**
 * @brief Given the low order IDP solution, compute the local min and max values
 * for the high order solution. this computes u_i^min and u_i^max defined in 
 * equation 2.1.
 *
 * @param rho_gf_lo the invariant-domain-preserving low order update
 */
void ComputeRhoMinMax(const ParGridFunction &gf_lo)
{
   // Continuous projection
   GridFunctionCoefficient LO_rho_gf_coeff(&gf_lo);

   /* Project local adjacent max to continuous gf */
   LO_proj_max.ProjectDiscCoefficientMax(LO_rho_gf_coeff);
   LO_proj_min.ProjectDiscCoefficientMin(LO_rho_gf_coeff);
   
   /* 
   Convert continuous space to high order dg space 
   The use of this operator assumes that the dofs in 
   the H1_LO space and H1_HO space are the same
   */
   vert_elem->Mult(LO_proj_max, rho_max);
   vert_elem->Mult(LO_proj_min, rho_min);
}

/**
 * @brief Limit the high-order approximation using the local maximum
 * defined by the IDP update. This algorithm is described in equations
 * (2.2) - (2.3).
 *
 * @param rho_gf_ho the high order update, not yet IDP
 */
void LimitMassMax(ParGridFunction &gf_ho)
{
   cout << "IDPLimiter::LimitMassMax\n";
   if (gf_ho.Size() != NDofs)
   {
      MFEM_ABORT("Incompatible sizes of gridfunctions provided.\n");
   }

   MFEM_ABORT("How to collect eligible neighbors?");

   /* Iterate over DoFs */
   for (int dof_it = 0; dof_it < NDofs; dof_it++)
   {
      /* Compute ai+ */

      /* Compute bi+ */

      /* Compute li+ */

      /* Set y_i */
      /* Set y_j for all neighbors */
   }
}

/**
 * @brief Limit the high-order approximation using the local minimum
 * defined by the IDP update. This algorithm is described in equations
 * (2.4) - (2.5).
 *
 * @param rho_gf_ho the high order update, not yet IDP
 */
void LimitMassMin(ParGridFunction &gf_ho)
{
   
}
   
}; // IDPLimiter
} // ns mfem

#endif // LIMITER