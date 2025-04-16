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

enum IterationMethod {
   JACOBI,
   GAUSS_SEIDEL
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

   /* Misc options */
   const int max_it = 2;
   bool use_glob = false;
   bool suppress_warnings = false;
   IterationMethod iter_method = GAUSS_SEIDEL; // options: JACOBI, GAUSS_SEIDEL

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

void GetRhoMax(Vector &rho_max_out) const { rho_max_out = rho_max; }
void GetRhoMin(Vector &rho_min_out) const { rho_min_out = rho_min; }
void GetRhoMaxMin(int i, double &rho_max_i, double &rho_min_i) const
{
   if (use_glob)
   {
      rho_max_i = glob_rho_max;
      rho_min_i = glob_rho_min;
   }
   else
   {
      rho_max_i = rho_max[i];
      rho_min_i = rho_min[i];
   }
}

bool CheckLocalMassConservation(const Vector &x, const Vector &y) const
{
   // cout << "IDPLimiter::CheckLocalMassConservation\n";
   Array<int> adj_dofs;
   double sum = 0.;
   bool is_locally_conservative = true;
   for (int i = 0; i < NDofs; i++)
   {
      /* check that mi(yi-xi) + SIGMA_j mj(yj-xj) = 0 */
      // cout << "checking local mass at dof: " << i << endl;
      sum = mass_vec[i] * (y[i] - x[i]);
      GetAdjacency(i, adj_dofs);
      // cout << "adj dofs: ";
      // adj_dofs.Print(cout);
      for (int adj_it = 0; adj_it < adj_dofs.Size(); adj_it++)
      {
         int j = adj_dofs[adj_it];
         sum += mass_vec[j] * (y[j] - x[j]);
      }

      if (abs(sum) > 1.e-12)
      {
         cout << "local mass not conserved at dof: " << i << ", sum: " << sum << endl;
         // cout << "adj dofs: ";
         // adj_dofs.Print(cout);
         is_locally_conservative = false;
         // cout << "Not locally conservative at dof: " << i << ", val: " << sum << endl;
         // cout << "mass_vec[i]: " << mass_vec[i] << ", x[i]: " << x[i] << ", y[i]: " << y[i] << endl;
         // cout << "mass_vec[j]: " << mass_vec[adj_dofs[adj_it]] << ", x[j]: " << x[adj_dofs[adj_it]] << ", y[j]: " << y[adj_dofs[adj_it]] << endl;
         break; // Break out of loop as there is not point looking at other dofs
      }
   }
   return is_locally_conservative;
}

bool CheckGlobalMassConservation(const Vector &mass_old, const Vector &x, const Vector &y) const
{
   // cout << "IDPLimiter::CheckGlobalMassConservation\n";
   bool is_globally_conservative = true;
   double sum_old = 0., sum_new = 0.;
   for (int i = 0; i < NDofs; i++)
   {
      /* Check if mass has been preserved */
      sum_old += mass_old[i] * x[i];
      sum_new += mass_vec[i] * y[i];
   }
   if (abs(sum_old - sum_new) > 1.e-12)
   {
      cout << "global mass is not conserved, sum_old: " << sum_old << ", sum_new: " << sum_new << endl;
      is_globally_conservative = false;
   }
   return is_globally_conservative;
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
 * @param gf_ho The high-order solution to be limited. This parameter is modified in-place.
 * 
 * The function performs the following steps:
 * 2. Iteratively adjusts `gf_ho` to ensure it stays within the bounds.
 * 3. Stops iterating when no further adjustments are made or the maximum number 
 *    of iterations (`max_it`) is reached.
 *
 * @note Jacobi method
 */
void LocalConservativeLimit(ParGridFunction &gf_ho)
{
   // cout << "IDPLimiter::LocalConservativeLimit\n";
   int num_max_lim_old, num_min_lim_old;
   int num_max_lim = 0, num_min_lim = 0;
   int num_it = 0;
   bool max_it_reached = false, done_iterating = false;
   while (!done_iterating && !max_it_reached)
   {
      num_min_lim_old = num_min_lim;
      num_max_lim_old = num_max_lim;

      num_min_lim = 0, num_max_lim = 0;

      /* Actual limit procedure */
      Array<int> adj_dofs;
      Vector y(NDofs);
      y = 0.;
      for (int dof_it = 0; dof_it < NDofs; dof_it++)
      {
         double mi = mass_vec[dof_it];
         double rho_max_i = rho_max[dof_it], rho_min_i = rho_min[dof_it];
         double rho_i = gf_ho[dof_it];
         GetAdjacency(dof_it, adj_dofs);
         if (mi == 0.) { y[dof_it] = min( max( rho_i, rho_min_i), rho_max_i); }
         else if (rho_i > rho_max_i)
         {
            num_max_lim++;
            /* Compute ai+ */
            double aip = 0.;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               aip += mass_vec[dof_j]*max(0., rho_max[dof_j] - gf_ho[dof_j]);
            }

            /* Compute bi+ */
            double bip = max(rho_i - aip / mi, rho_max_i);

            /* Compute li+ */
            double lip = 0.;
            if (aip > 0.) { lip = mi * (rho_i - bip) / aip; }
            else if (aip < 0.) { MFEM_ABORT("aim <= 0."); }

            /* Set y */
            y[dof_it] = bip;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               y[dof_j] = gf_ho[dof_j] + lip * max(0., rho_max[dof_j] - gf_ho[dof_j]);
            }
         }
         else if (rho_i < rho_min_i)
         {
            num_min_lim++;
            /* Compute ai- */
            double aim = 0.;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               aim += mass_vec[dof_j]*std::max(0., gf_ho[dof_j] - rho_min[dof_j]);
            }

            /* Compute bi+ */
            double bim = min(rho_i + aim / mi, rho_min_i);

            /* Compute li+ */
            double lim = 0.;
            if (aim > 0.)
            {
               lim = mi * (rho_i - bim) / aim;
            }
            else if (aim < 0.)
            {
               MFEM_ABORT("aim <= 0.");
            }

            /* Set y */
            y[dof_it] = bim;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               y[dof_j] = gf_ho[dof_j] + lim * std::max(0., gf_ho[dof_j] - rho_min[dof_j]);
            }
         }
         else { y[dof_it] = rho_i; }
      }

      /* End actual limit procedure */

      /* Check if any adjustment has been made */
      if (num_min_lim == num_min_lim_old && num_max_lim == num_max_lim_old)
      {
         /* No improvement from last iteratior */
         done_iterating = true;
      }

      // /* Check for max it reached */
      num_it++;
      if (num_it == max_it) { max_it_reached = true;}

      // cout << "num min limited: " << num_min_lim << endl;
      // cout << "num max limited: " << num_max_lim << endl;

      bool is_locally_conservative = CheckLocalMassConservation(gf_ho, y);
      if (!is_locally_conservative && !suppress_warnings)
      {
         MFEM_WARNING("Jacobi limiter not locally conservative.\n");
      }

      gf_ho = y;
   }
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
 * @param gf_ho The high-order solution to be limited. This parameter is modified in-place.
 * 
 * The function performs the following steps:
 * 2. Iteratively adjusts `gf_ho` to ensure it stays within the bounds.
 * 3. Stops iterating when no further adjustments are made or the maximum number 
 *    of iterations (`max_it`) is reached.
 *
 * @note Gauss-Seidel method
 */
void LocalConservativeLimitGS(ParGridFunction &gf_ho)
{
   // cout << "IDPLimiter::LocalConservativeLimit\n";
   int num_max_lim_old, num_min_lim_old;
   int num_max_lim = 0, num_min_lim = 0;
   int num_it = 0;
   bool max_it_reached = false, done_iterating = false;

   /* Iterate max_it times or when there is no local budget for any of the limited quantities */
   while (!done_iterating && !max_it_reached)
   {
      num_min_lim_old = num_min_lim;
      num_max_lim_old = num_max_lim;

      num_min_lim = 0, num_max_lim = 0;
      int num_min = 0, num_max = 0;

      /* Actual limit procedure */
      Array<int> adj_dofs;
      Vector x_old(NDofs);
      x_old = gf_ho;
      for (int dof_it = 0; dof_it < NDofs; dof_it++)
      {
         double mi = mass_vec[dof_it];
         double rho_max_i = rho_max[dof_it], rho_min_i = rho_min[dof_it];
         double rho_i = gf_ho[dof_it];
         GetAdjacency(dof_it, adj_dofs);
         // cout << "i: " << dof_it << ", adj_dofs: ";
         // adj_dofs.Print(cout);
         if (mi == 0.) { gf_ho[dof_it] = min( max( rho_i, rho_min_i), rho_max_i); }
         else if (rho_i > rho_max_i)
         {
            num_max++;
            /* Compute ai+ */
            double aip = 0.;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               aip += mass_vec[dof_j]*max(0., rho_max[dof_j] - gf_ho[dof_j]);
            }

            if (aip > 0.)
            {
               // cout << "aip > 0\n";
               num_max_lim++;
               /* Compute bi+ and li+*/
               double bip = max(rho_i - aip / mi, rho_max_i);
               double lip = mi * (rho_i - bip) / aip;

               // cout << "aip: " << aip << ", bip: " << bip << ", lip: " << lip << endl;

               /* Check mass transfer is 0 */
               double mass_delta = mi * (bip - rho_i);

               /* Set gf_ho */
               gf_ho[dof_it] = bip;
               for (int j = 0; j < adj_dofs.Size(); j++)
               {
                  int dof_j = adj_dofs[j];
                  double _new_val = gf_ho[dof_j] + lip * max(0., rho_max[dof_j] - gf_ho[dof_j]);
                  mass_delta += mass_vec[dof_j] * (_new_val - gf_ho[dof_j]);
                  gf_ho[dof_j] = _new_val;
               }
               if (abs(mass_delta) > 1.e-12)
               {
                  MFEM_ABORT("Mass delta should be 0.");
               }
            }
            
         }
         else if (rho_i < rho_min_i)
         {
            num_min++;
            /* Compute ai- */
            double aim = 0.;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               aim += mass_vec[dof_j]*max(0., gf_ho[dof_j] - rho_min[dof_j]);
            }

            if (aim > 0.)
            {
               // cout << "aim > 0\n";
               num_min_lim++;

               /* Compute bi- and li- */
               double bim = min(rho_i + aim / mi, rho_min_i);
               double lim = mi * (rho_i - bim) / aim;

               // cout << "aim: " << aim << ", bim: " << bim << ", lim: " << lim << endl;

               /* Check mass transfer is 0 */
               double mass_delta = mi * (bim - rho_i);

               /* Set gf_ho */
               gf_ho[dof_it] = bim;
               for (int j = 0; j < adj_dofs.Size(); j++)
               {
                  int dof_j = adj_dofs[j];
                  double _new_val = gf_ho[dof_j] + lim * max(0., gf_ho[dof_j] - rho_min[dof_j]);
                  mass_delta += mass_vec[dof_j] * (_new_val - gf_ho[dof_j]);
                  gf_ho[dof_j] = _new_val;
               }

               if (abs(mass_delta) > 1.e-12)
               {
                  MFEM_ABORT("Mass delta should be 0.");
               }
            }
         }
      }

      /* End actual limit procedure */

      /* Check if any adjustment has been made */
      if (num_min_lim == num_min_lim_old && num_max_lim == num_max_lim_old)
      {
         /* No improvement from last iteratior */
         done_iterating = true;
      }

      // /* Check for max it reached */
      num_it++;
      if (num_it == max_it) { max_it_reached = true;}

      // cout << "ratio of min able to be limited: " << double(num_min_lim)/num_min << endl;
      // cout << "ratio of max able to be limited: " << double(num_max_lim)/num_max << endl;

      /* Check conservation */
      bool is_locally_conservative = CheckLocalMassConservation(x_old, gf_ho);
      if (!is_locally_conservative && !suppress_warnings)
      {
         MFEM_WARNING("Gauss-Seidel limiter not locally conservative.\n");
      }
      bool is_globally_conservative = CheckGlobalMassConservation(mass_vec, x_old, gf_ho);
      if (!is_globally_conservative && !suppress_warnings)
      {
         MFEM_WARNING("Gauss-Seidel limiter not globally conservative.\n");
      }
   }
}


/**
 * @brief Applies a global conservative limiter to the given ParGridFunction.
 *
 * This function ensures that the values of the input ParGridFunction `gf_ho` 
 * are limited to a specified range [glob_rho_min, glob_rho_max] while 
 * preserving the total mass. The limiting process involves computing 
 * intermediate values and applying correction factors to maintain 
 * conservation properties.
 *
 * @param gf_ho The high-order ParGridFunction to be limited. It is modified
 *              in-place to satisfy the limiting constraints.
 *
 * @pre `gf_ho.Size()` must be equal to `NDofs`.
 * @pre `glob_rho_min` must be non-negative, and `glob_rho_max` must be 
 *      greater than or equal to `glob_rho_min`.
 *
 * @post The values of `gf_ho` are within the range [glob_rho_min, glob_rho_max].
 * @post The total mass of `gf_ho` is conserved.
 *
 * @details
 * - The function first computes the total mass of the input function.
 * - It then calculates a limited version of the input values (`y`) by 
 *   clamping them to the range [rho_min_i, rho_max_i] for each degree of freedom.
 * - Correction factors (`alpha_p` and `alpha_n`) are computed to adjust the 
 *   limited values while preserving the total mass.
 * - The final limited values (`z`) are computed and assigned back to `gf_ho`.
 * - The function checks for local mass conservation and issues a warning 
 *   if it is not satisfied.
 *
 * @warning If the function detects that the global limiter is not locally 
 *          conservative, a warning is issued.
 *
 * @throws An assertion failure if `alpha_p` or `alpha_n` are not within 
 *         the range [0, 1].
 * @throws An assertion failure if the final values of `gf_ho` are not 
 *         within the range [glob_rho_min, glob_rho_max].
 */
void GlobalConservativeLimit(ParGridFunction &gf_ho)
{
   // cout << "===== Limiter::GlobalConservativeLimit =====\n";
   assert(gf_ho.Size() == NDofs);
   assert(glob_rho_min >= 0. && glob_rho_max >= glob_rho_min);

   /* Compute total mass */
   double M = 0.;
   double rho_max_i, rho_min_i;
   for (int i = 0; i < NDofs; i++)
   {
      M += mass_vec[i] * gf_ho[i]; 
   }
   // cout << "M: " << M << endl;

   /* Compute y */
   Vector y(NDofs);
   for (int i = 0; i < NDofs; i++)
   {
      GetRhoMaxMin(i, rho_max_i, rho_min_i);
      // cout << "i " << i << " rho max " << rho_max_i << " rho min " << rho_min_i << endl;
      y[i] = std::min( std::max( gf_ho[i], rho_min_i), rho_max_i);
   }

   /* Compute alpha+ and alpha-*/
   double _sum_my = 0., _denom_max = 0., _denom_min = 0.;
   for (int i = 0; i < NDofs; i++)
   {
      GetRhoMaxMin(i, rho_max_i, rho_min_i);
      double _m = mass_vec[i];
      _sum_my += _m * y[i];
      _denom_max += _m * (rho_max_i - y[i]);
      _denom_min += _m * (rho_min_i - y[i]);
   }
   const double _num = M - _sum_my;
   const double alpha_p = std::max(0., _num / _denom_max);
   const double alpha_n = std::max(0., _num / _denom_min);

   MFEM_ASSERT(alpha_p >= 0. && alpha_p <= 1. && 
               alpha_n >= 0. && alpha_n <= 1.,
               "alpha_p and alpha_n should be between 0 and 1.\n");

   /* Finally, set z */
   Vector z(NDofs);
   for (int i = 0; i < NDofs; i++)
   {
      GetRhoMaxMin(i, rho_max_i, rho_min_i);
      z[i] = y[i] + alpha_p * (rho_max_i - y[i]) + alpha_n * (rho_min_i - y[i]);
   }

   /* Checks */
   bool is_globally_conservative = CheckGlobalMassConservation(mass_vec, gf_ho, z);
   if (!is_globally_conservative && !suppress_warnings)
   {
      MFEM_WARNING("Global limiter not locally conservative.\n");
   }

   MFEM_ASSERT(gf_ho.Max() <= glob_rho_max && gf_ho.Min() >= glob_rho_min,
               "Global limiting failed.\n");
   /* end checks */

   gf_ho = z;
}

bool CheckEstimate(const ParGridFunction &gf_ho)
{
   // cout << "IDPLimiter::CheckEstimate\n";
   double M_max = 0., M_min = 0., M = 0., mi;
   for (int i = 0; i < NDofs; i++)
   {
      mi = mass_vec[i];
      M += mi * gf_ho[i];
      M_max += mi * rho_max[i];
      M_min += mi * rho_min[i];
   }
   // cout << "M_min: " << M_min << ", M: " << M << ", M_max: " << M_max << endl;
   if (M >= M_min && M <= M_max)
   {
      return true;
   }
   else
   {
      return false;
   }
}

void Limit(const ParGridFunction &gf_lo, ParGridFunction &gf_ho)
{
   // cout << "===== Limiter::LimitGlobal =====\n";
   ComputeRhoMinMax(gf_lo);

   bool satisfies_estimate = CheckEstimate(gf_ho);
   if (!satisfies_estimate)
   {
      MFEM_ABORT("Bounding estimate not satisfied");
   }

   switch(iter_method)
   {
      case JACOBI:
         LocalConservativeLimit(gf_ho);
         GlobalConservativeLimit(gf_ho);
         break;
      case GAUSS_SEIDEL:
         LocalConservativeLimitGS(gf_ho);
         GlobalConservativeLimit(gf_ho);
         break;
      default:
         MFEM_ABORT("Unknown iteration method.\n");
         break;
   }
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
void GetAdjacency(const int dof, Array<int> &adj_dofs) const
{
   // cout << "----------IDPLimiter::GetAdjacency----------\n";

   /* first retrieve HO dofs that share the same max/min given by the high order space */
   // Vector LO_vals, adj_vals; // vector isn't used
   // Array<int> LO_verts;
   // _spm.GetRow(dof, LO_verts, LO_vals); // val isn't used
   // assert(LO_verts.Size() == 1); // Each HO dof should only correspond to one LO H1 dof
   // _spm_trans.GetRow(LO_verts[0], adj_dofs, adj_vals);
   // assert((adj_dofs.Size() <= 4) && "Cannot have more than 4 values sharing a vertex.");

   /* Next, append HO dofs that are in the same mesh element */
   int el_index = pfes.GetElementForDof(dof);
   Array<int> el_dofs;
   pfes.GetElementDofs(el_index, el_dofs);
   // adj_dofs.Append(el_dofs);
   adj_dofs = el_dofs;

   /* Remove dof from adj_dofs */
   adj_dofs.DeleteFirst(dof);
   adj_dofs.DeleteFirst(dof);
}

   
}; // IDPLimiter
} // ns mfem

#endif // LIMITER