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
   double glob_rho_max = -1., glob_rho_min = -1.;

   ParFiniteElementSpace H1_LO, H1_HO;
   ParGridFunction LO_proj_max, LO_proj_min;
   ParDiscreteLinearOperator vert_elem_oper; 
   HypreParMatrix *vert_elem, *vert_elem_trans;

   SparseMatrix _spm, _spm_trans;

   const int max_it = 5;

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
   // _spm & _spm_trans will be used to determine adjacency
   vert_elem_oper.AddDomainInterpolator(new IdentityInterpolatorMax);
   vert_elem_oper.Assemble();
   vert_elem_oper.Finalize();
   vert_elem = vert_elem_oper.ParallelAssemble();
   vert_elem_trans = vert_elem->Transpose();

   vert_elem->MergeDiagAndOffd(_spm);
   vert_elem_trans->MergeDiagAndOffd(_spm_trans);
}

~IDPLimiter() 
{
   // delete spm;
   delete vert_elem;
}


/**
 * @brief Applies a local conservative limiting process to ensure the high-order 
 *        solution remains within specified bounds.
 * 
 * This function iteratively adjusts the high-order solution (`gf_ho`) to ensure 
 * it satisfies local conservation properties while staying within the minimum 
 * and maximum bounds computed from the low-order solution (`gf_lo`). The process 
 * continues until no further adjustments are needed or a maximum number of 
 * iterations is reached.
 * 
 * @param gf_lo The low-order solution used to compute the minimum and maximum bounds.
 * @param gf_ho The high-order solution to be limited. This parameter is modified in-place.
 * 
 * The function performs the following steps:
 * 1. Computes the minimum and maximum bounds (`rho_min` and `rho_max`) from `gf_lo`.
 * 2. Iteratively adjusts `gf_ho` using `LimitMassMax` and `LimitMassMin` to ensure 
 *    it stays within the bounds.
 * 3. Stops iterating when no further adjustments are made or the maximum number 
 *    of iterations (`max_it`) is reached.
 * 
 * Debugging output is provided to indicate the percentage of the solution that 
 * was limited in each iteration. Uncommented sections of the code can be used 
 * for additional debugging to verify that all degrees of freedom (DoFs) in `gf_ho` 
 * remain within the computed bounds.
 */
void LocalConservativeLimit(const ParGridFunction &gf_lo, ParGridFunction &gf_ho)
{
   // cout << "IDPLimiter::LocalConservativeLimit\n";
   ComputeRhoMinMax(gf_lo);

   double pct_max_lim = 1., pct_min_lim = 1.;
   double pct_max_lim_old = 0., pct_min_lim_old = 0.;
   int num_it = 0;
   bool max_it_reached = false, done_iterating = false;
   while (!done_iterating && !max_it_reached)
   {
      pct_min_lim_old = pct_min_lim;
      pct_max_lim_old = pct_max_lim;

      LimitMassMax(gf_ho, pct_max_lim);
      LimitMassMin(gf_ho, pct_min_lim);

      /* Check if any adjustment has been made */
      if (pct_min_lim == pct_min_lim_old && pct_max_lim == pct_max_lim_old)
      {
         /* No improvement from last iteratior */
         done_iterating = true;
      }

      /* Check for max it reached */
      num_it++;
      if (num_it == max_it) { max_it_reached = true;}

      // if(pct_max_lim + pct_min_lim > 0.)
      // {
      //    cout << "pct max limited: " << pct_max_lim << endl;
      //    cout << "pct min limited: " << pct_min_lim << endl;
      // }
   }
   // for (int i = 0; i < NDofs; i++)
   // {
   //    if (gf_ho[i] < rho_min[i] || gf_ho[i] > rho_max[i])
   //    {
   //       cout << "dof: " << i << ", min: " << rho_min[i] << ", ho: " << gf_ho[i] << ", max: " << rho_max[i] << endl;
   //    }
   // }
   
}

void GlobalConservativeLimit(ParGridFunction &gf_ho)
{
   // cout << "===== Limiter::GlobalConservativeLimit =====\n";
   assert(gf_ho.Size() == NDofs);
   assert(glob_rho_min >= 0. && glob_rho_max >= glob_rho_min);

   /* Compute total mass */
   double M = 0.;
   for (int i = 0; i < NDofs; i++)
   {
      M += mass_vec[i] * gf_ho[i]; 
   }
   // cout << "M: " << M << endl;

   /* Compute y */
   Vector y(NDofs);
   for (int i = 0; i < NDofs; i++)
   {
      y[i] = std::min( std::max( gf_ho[i], glob_rho_min), glob_rho_max );
   }

   /* Compute alpha+ and alpha-*/
   double _sum_my = 0., _denom_max = 0., _denom_min = 0.;
   for (int i = 0; i < NDofs; i++)
   {
      double _m = mass_vec[i];
      _sum_my += _m * y[i];
      _denom_max += _m * (glob_rho_max - y[i]);
      _denom_min += _m * (glob_rho_min - y[i]);
   }
   const double _num = M - _sum_my;
   const double alpha_p = std::max(0., _num / _denom_max);
   const double alpha_n = std::max(0., _num / _denom_min);

   /* Finally, set z */
   Vector z(NDofs);
   for (int i = 0; i < NDofs; i++)
   {
      z[i] = y[i] + alpha_p * (glob_rho_max - y[i]) + alpha_n * (glob_rho_min - y[i]);
   }

   gf_ho = z;
}

void LimitGlobal(const ParGridFunction &gf_lo, ParGridFunction &gf_ho)
{
   // cout << "===== Limiter::LimitGlobal =====\n";
   LocalConservativeLimit(gf_lo, gf_ho);

   GlobalConservativeLimit(gf_ho);
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
   // cout << "IDPLimiter::ComputeRhoMinMax\n";
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

   /* Global bounds necessary for global conservative limiting, see Remark 2.4 */
   glob_rho_max = rho_max.Max();
   glob_rho_min = rho_min.Min();
}


/**
 * @brief Retrieves the adjacency information for a given degree of freedom (DOF).
 *
 * This function populates the provided array with the adjacent LO degrees of freedom
 * (DOFs) for a given input HO DOF. 
 *
 * @param dof The HO degree of freedom (DOF) for which adjacency information is required.
 * @param adj_dofs The set of HO degrees of freedom (DOF) that dof can exchange mass with.
 *
 *
 * Example:
 *
 *    *-----*-----*
 *    |6   7|10 11|
 *    |     |     |
 *    |4   5|8   9|
 *    *-----*-----*
 *    |2   3|14 15|
 *    |     |     |
 *    |0   1|12 13|
 *    *-----*-----*
 *
 *    ------------------------
 *    | dof | adj_dofs       |
 *    |  5  | 3,4,5,6,7,8,14 |
 *    |  13 | 12,13,14,15    |
 *    |  4  | 2,4,5,6,7      |
 *    |  10 | 7.8.9.10,11    |
 *    ------------------------ 
 */
void GetAdjacency(const int dof, Array<int> &adj_dofs)
{
   // cout << "----------IDPLimiter::GetAdjacency----------\n";

   /* first retrieve HO dofs that share the same max/min given by the high order space */
   Vector LO_vals, adj_vals; // vector isn't used
   Array<int> LO_verts;
   _spm.GetRow(dof, LO_verts, LO_vals); // val isn't used
   assert(LO_verts.Size() == 1); // Each HO dof should only correspond to one LO H1 dof
   _spm_trans.GetRow(LO_verts[0], adj_dofs, adj_vals);
   assert((adj_dofs.Size() <= 4) && "Cannot have more than 4 values sharing a vertex.");

   /* Next, append HO dofs that are in the same mesh element */
   int el_index = pfes.GetElementForDof(dof);
   Array<int> el_dofs;
   pfes.GetElementDofs(el_index, el_dofs);
   adj_dofs.Append(el_dofs);

   /* Remove dof from adj_dofs */
   adj_dofs.DeleteFirst(dof);
   adj_dofs.DeleteFirst(dof);
}

/**
 * @brief Limit the high-order approximation using the local maximum
 * defined by the IDP update. This algorithm is described in equations
 * (2.2) - (2.3).
 *
 * @param rho_gf_ho the high order update, not yet IDP
 */
void LimitMassMax(ParGridFunction &gf_ho, double &pct_limited)
{
   int num_limited = 0;
   // cout << "IDPLimiter::LimitMassMax\n";
   if (gf_ho.Size() != NDofs)
   {
      MFEM_ABORT("Incompatible sizes of gridfunctions provided.\n");
   }

   Array<int> adj_dofs;

   /* Iterate over DoFs */
   for (int dof_it = 0; dof_it < NDofs; dof_it++)
   {
      if (gf_ho[dof_it] <= rho_max[dof_it])
      {
         /* Maximum principle is satisfied */
         continue;
      }
      num_limited++;

      /* Get adjacency */
      GetAdjacency(dof_it, adj_dofs);

      /* Compute ai+ */
      double aip = 0.;
      for (int j = 0; j < adj_dofs.Size(); j++)
      {
         int dof_j = adj_dofs[j];
         aip += mass_vec[dof_j]*std::max(0., rho_max[dof_j] - gf_ho[dof_j]);
      }

      /* Compute bi+ */
      double bip = max(gf_ho[dof_it] - aip / mass_vec[dof_it], rho_max[dof_it]);

      /* Compute li+ */
      double lip = 0.;
      if (aip > 0.)
      {
         lip = mass_vec[dof_it] * (gf_ho[dof_it] - bip) / aip;
      }
      else if (aip < 0.)
      {
         MFEM_ABORT("aim <= 0.");
      }

      /* Set y_i */
      gf_ho[dof_it] = bip;

      /* Set y_j for all neighbors */
      for (int j = 0; j < adj_dofs.Size(); j++)
      {
         int dof_j = adj_dofs[j];
         gf_ho[dof_j] = gf_ho[dof_j] + lip * std::max(0., rho_max[dof_j] - gf_ho[dof_j]);
      }
   }
   pct_limited = double(num_limited) / NDofs;
}

/**
 * @brief Limit the high-order approximation using the local minimum
 * defined by the IDP update. This algorithm is described in equations
 * (2.4) - (2.5).
 *
 * @param rho_gf_ho the high order update, not yet IDP
 */
void LimitMassMin(ParGridFunction &gf_ho, double &pct_limited)
{
   int num_limited = 0;
   // cout << "IDPLimiter::LimitMassMax\n";
   if (gf_ho.Size() != NDofs)
   {
      MFEM_ABORT("Incompatible sizes of gridfunctions provided.\n");
   }

   Array<int> adj_dofs;

   /* Iterate over DoFs */
   for (int dof_it = 0; dof_it < NDofs; dof_it++)
   {
      if (gf_ho[dof_it] >= rho_min[dof_it])
      {
         /* Maximum principle is satisfied */
         continue;
      }
      num_limited++;

      /* Get adjacency */
      GetAdjacency(dof_it, adj_dofs);

      /* Compute ai- */
      double aim = 0.;
      for (int j = 0; j < adj_dofs.Size(); j++)
      {
         int dof_j = adj_dofs[j];
         aim += mass_vec[dof_j]*std::max(0., gf_ho[dof_j] - rho_min[dof_j]);
      }

      /* Compute bi+ */
      double bim = min(gf_ho[dof_it] + aim / mass_vec[dof_it], rho_min[dof_it]);

      /* Compute li+ */
      double lim = 0.;
      if (aim > 0.)
      {
         lim = mass_vec[dof_it] * (gf_ho[dof_it] - bim) / aim;
      }
      else if (aim < 0.)
      {
         MFEM_ABORT("aim <= 0.");
      }

      /* Set y_i */
      gf_ho[dof_it] = bim;

      /* Set y_j for all neighbors */
      for (int j = 0; j < adj_dofs.Size(); j++)
      {
         int dof_j = adj_dofs[j];
         gf_ho[dof_j] = gf_ho[dof_j] + lim * std::max(0., gf_ho[dof_j] - rho_max[dof_j]);
      }
   }
   pct_limited =  double(num_limited) / NDofs;
}
   
}; // IDPLimiter
} // ns mfem

#endif // LIMITER