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

/**
 * @brief Interpolator used to find adjacency information between finite element spaces.
 *
 * The AdjacencyInterpolator class inherits from DiscreteInterpolator and overrides
 * the AssembleElementMatrix2 method to create a dense matrix where every entry is 1.0.
 * This can be used for adjacency-based interpolation or as a placeholder for more
 * complex interpolation schemes.
 *
 * Source: https://github.com/mfem/mfem/issues/436
 *
 * How to define operator:
 *    ParMesh *pmesh;
 *    // define pmesh ...
 *    L2_FECollection elem_fec(0, pmesh->Dimension());
 *    H1_FECollection vert_fec(1, pmesh->Dimension());
 *    ParFiniteElementSpace elem_fes(pmesh, &elem_fec);
 *    ParFiniteElementSpace vert_fes(pmesh, &vert_fec);
 *    ParDiscreteLinearOperator vert_elem_oper(&vert_fes, &elem_fes); // maps vert_fes to elem_fes
 *    vert_elem_oper.AddDomainInterpolator(new AdjacencyInterpolator);
 *    vert_elem_oper.Assemble();
 *    HypreParMatrix *vert_elem = vert_elem_oper.ParallelAssemble();
 *    // use vert_elem ...
 *    delete vert_elem;
 *
 * Example usage:
 *    Vector elem_marker(elem_fes.GetTrueVSize());
 *    elem_marker = 0.0;
 *    // set elem_marker(i) to 1.0 if element i is an element whose neighbors you want to find
 *    Vector vert_marker(vert_fes.GetTrueVSize());
 *    vert_elem->MultTranspose(elem_marker, vert_marker);
 *    vert_elem->Mult(vert_marker, elem_marker);
 *    // if elem_marker(j) is > 0.0 then j is a neighbor of an element you marked above, or
 *    // j itself was marked above
 *
 * Matrix usage:
 *   This operator can also be used to build an el adjacency matrix, as is done in the Limiter
 *   class, as follows:
 *
 *    vert_elem_oper.AddDomainInterpolator(new AdjacencyInterpolator);
 *    vert_elem_oper.Assemble();
 *    vert_elem_oper.Finalize();
 *    vert_elem = vert_elem_oper.ParallelAssemble();
 *    vert_elem_trans = vert_elem->Transpose();
 *    vert_elem->MergeDiagAndOffd(_spm);
 *    vert_elem_trans->MergeDiagAndOffd(_spm_trans);
 *    el_adj_mat = Mult(_spm, _spm_trans);
 *
 * @note
 * The returned elem_marker vector will contain:
 * - 1.0 for elements that are diagonal neighbors of the marked element
 * - 2.0 for elements that are face neighbords of the marked element
 * - 4.0 elements that were marked in the first place
 * - 0.0 for all other elements
 */
class AdjacencyInterpolator : public DiscreteInterpolator
{
public:
   virtual void AssembleElementMatrix2(const FiniteElement &dom_fe,
                                       const FiniteElement &ran_fe,
                                       ElementTransformation &Trans,
                                       DenseMatrix &elmat)
   { elmat.SetSize(ran_fe.GetDof(), dom_fe.GetDof()); elmat = 1.0; }
};

enum IterationMethod {
   JACOBI,
   GAUSS_SEIDEL
};


class IDPLimiter
{
private:
   ParFiniteElementSpace &L2_HO, &L2_LO;
   Vector *mass_vec;
   const int NDofs;
   ParGridFunction vol_vec;
   mutable ParGridFunction x_min, x_max;
   mutable ParGridFunction x_min_relaxed, x_max_relaxed;
   double glob_x_max = -1., glob_x_min = -1.;

   ParFiniteElementSpace H1_LO, H1_HO;
   ParGridFunction LO_proj_max, LO_proj_min;
   ParDiscreteLinearOperator vert_elem_oper; 
   HypreParMatrix *vert_elem, *vert_elem_trans;
   SparseMatrix _spm, _spm_trans, *el_adj_mat;
   const IntegrationRule &ir;

   /* stiffness */
   HypreParMatrix *stiffness_mat;
   SparseMatrix *stiffness_mat_sp;

   /* Misc options */
   const int dim;
   const int max_it = 2;
   bool use_global_conservative_limiting = false;
   bool use_global_bounds = true; // use global bounds for limiting
   const int relaxation_option = 2; // 0: no relaxation, 
                                    // 1: relaxation type 1, 
                                    // 2: Guermond-Wang A.3
   bool suppress_warnings = false;
   bool suppress_output = true;
   IterationMethod iter_method = GAUSS_SEIDEL; // options: JACOBI, GAUSS_SEIDEL

   mutable bool local_bounds_computed = false;
   mutable bool glob_bounds_computed = false;

public:
IDPLimiter(ParFiniteElementSpace & l2_ho, ParFiniteElementSpace & l2_lo, ParFiniteElementSpace & H1FESpace_proj_LO,
           ParFiniteElementSpace & H1FESpace_proj_HO, Vector &_mass_vec, const int &oq) :
   L2_HO(l2_ho),
   L2_LO(l2_lo),
   mass_vec(&_mass_vec),
   NDofs(mass_vec->Size()),
   vol_vec(&l2_ho),
   x_min(&l2_ho),
   x_max(&l2_ho),
   x_min_relaxed(&l2_ho),
   x_max_relaxed(&l2_ho),
   H1_LO(H1FESpace_proj_LO),
   H1_HO(H1FESpace_proj_HO),
   LO_proj_max(&H1_LO),
   LO_proj_min(&H1_LO),
   vert_elem_oper(&H1FESpace_proj_LO, &L2_LO),
   ir(IntRules.Get(L2_HO.GetParMesh()->GetElementBaseGeometry(0),
                   (oq > 0) ? oq : 3 * H1_HO.GetOrder(0) + L2_HO.GetOrder(0) - 1)),
   dim(L2_HO.GetParMesh()->Dimension())
{
   PrintOptions();
   vol_vec = 0.;
   /* Verify that #DOFS are equal across LO and HO H1 spaces */
   if (H1_HO.GetNDofs() != H1_LO.GetNDofs())
   {
      MFEM_ABORT("Continuous approximation spaces are not compatible.\n");
   }

   /* Construct adjacency operator */
   vert_elem_oper.AddDomainInterpolator(new AdjacencyInterpolator);
   vert_elem_oper.Assemble();
   vert_elem_oper.Finalize();
   vert_elem = vert_elem_oper.ParallelAssemble();
   vert_elem_trans = vert_elem->Transpose();
   vert_elem->MergeDiagAndOffd(_spm);
   vert_elem_trans->MergeDiagAndOffd(_spm_trans);
   el_adj_mat = Mult(_spm, _spm_trans);
   /* Can further use SparseMatrix::Threshold to select which neighbors would like to keep */

   /* Construct stiffness matrix for relaxation */
   if (relaxation_option == 2)
   {
      ParBilinearForm stiffness(&L2_HO);
      ConstantCoefficient one_coeff(1.0);
      stiffness.AddDomainIntegrator(new DiffusionIntegrator(one_coeff, &ir));
      stiffness.Assemble();
      stiffness.Finalize();
      stiffness_mat = stiffness.ParallelAssemble();
      MFEM_ASSERT(stiffness_mat != nullptr, "Stiffness matrix must be assembled before relaxing bounds.");
      stiffness_mat_sp = new SparseMatrix(NDofs, NDofs); 
      stiffness_mat->MergeDiagAndOffd(*stiffness_mat_sp);
   }
}

~IDPLimiter() 
{
   // delete spm;
   delete vert_elem;

   // stiffness
   if (stiffness_mat_sp)
   {
      delete stiffness_mat_sp;
      stiffness_mat_sp = nullptr;
   }
}

void GetVolVec(Vector &vol_vec_out) const { vol_vec_out = vol_vec;}
void GetRhoMax(Vector &x_max_out) const { x_max_out = x_max; }
void GetRhoMin(Vector &x_min_out) const { x_min_out = x_min; }
void GetRhoMaxMin(int i, double &x_max_i, double &x_min_i) const
{
   x_max_i = x_max[i];
   x_min_i = x_min[i];
}


void PrintOptions()
{
   cout << "-----IDPLimiter options-----\n";
   cout << "  max_it: " << max_it << endl;
   cout << "  use_global_conservative_limiting: " << use_global_conservative_limiting << endl;
   cout << "  use_global_bounds: " << use_global_bounds << endl;
   cout << "  relaxation_option: " << relaxation_option << endl;
   cout << "  suppress_warnings: " << suppress_warnings << endl;
   cout << "  suppress_output: " << suppress_output << endl;
   cout << "  iter_method: ";
   if (iter_method == JACOBI) { cout << "JACOBI" << endl; }
   else if (iter_method == GAUSS_SEIDEL) { cout << "GAUSS_SEIDEL" << endl; }
   else { cout << "Unknown iteration method!" << endl; }
   cout << "----------------------------\n";
}


bool CheckCellMassConservation(const Vector &_mi_vec, const Vector &x_old, const Vector &x_new)
{
   // cout << "IDPLimiter::CheckCellMassConservation\n";
   bool is_cell_mass_conserved = true;
   for (int i = 0; i < L2_HO.GetNE(); i++)
   {
      double old_mass = 0., new_mass = 0.;
      Array<int> dofs;
      L2_HO.GetElementDofs(i, dofs);
      for (int dof_it = 0; dof_it < dofs.Size(); dof_it++)
      {
         int dof = dofs[dof_it];
         old_mass += _mi_vec.Elem(dof) * x_old[dof];
         new_mass += _mi_vec.Elem(dof) * x_new[dof];
      }
      double t_val = abs(old_mass - new_mass);
      if (t_val > 1.e-12)
      {
         is_cell_mass_conserved = false;
         cout << "cell mass not conserved at element: " << i << ", diff: " << t_val << endl;
         cout << "old mass: " << old_mass << ", new mass: " << new_mass << endl;
      }
   }
   return is_cell_mass_conserved;
}


bool CheckGlobalMassConservation(const Vector &_mi_vec, const Vector &x, const Vector &y) const
{
   // cout << "IDPLimiter::CheckGlobalMassConservation\n";
   bool is_globally_conservative = true;
   double sum_old = 0., sum_new = 0.;
   for (int i = 0; i < NDofs; i++)
   {
      /* Check if mass has been preserved */
      sum_old += _mi_vec.Elem(i) * x[i];
      sum_new += _mi_vec.Elem(i) * y[i];
   }
   double t_val = abs(sum_old - sum_new) / sum_old;
   if (t_val > 1.e-12)
   {
      cout << "global mass is not conserved, diff: " << t_val << ", sum_old: " << sum_old << ", sum_new: " << sum_new << endl;
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
void LocalConservativeLimitJacobi(const Vector &_mi_vec, ParGridFunction &x)
{
   // cout << "IDPLimiter::LocalConservativeLimitJacobi\n";
   MFEM_ASSERT(local_bounds_computed, "Local bounds must be computed before applying local conservative limit.");
   MFEM_ABORT("Change adjacency to use cell dofs.\n");
   int num_max_lim_old, num_min_lim_old;
   int num_max_lim = 0, num_min_lim = 0;
   int num_it = 0;
   bool max_it_reached = false, done_iterating = false;
   while (!done_iterating && !max_it_reached)
   {
      num_min_lim_old = num_min_lim;
      num_max_lim_old = num_max_lim;

      num_min_lim = 0, num_max_lim = 0;
      int num_min = 0, num_max = 0;

      /* Actual limit procedure */
      Array<int> adj_dofs;
      Vector y(NDofs), dx(NDofs);
      y = x;
      dx = 0.;
      for (int dof_it = 0; dof_it < NDofs; dof_it++)
      {
         double _mi = _mi_vec.Elem(dof_it);
         double xi_max = x_max[dof_it], xi_min = x_min[dof_it];
         double xi = y[dof_it];
         GetAdjacency(dof_it, adj_dofs);

         if (_mi == 0.) { y[dof_it] = fmin( fmax( xi, xi_min), xi_max); }
         else if (xi > xi_max)
         {
            num_max++;
            /* Compute ai+ */
            double aip = 0.;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               aip += _mi_vec.Elem(dof_j)*fmax(0., x_max[dof_j] - y[dof_j]);
            }
            if (aip > 0.)
            {
               num_max_lim++;
               /* Compute bi+ and li+ */
               double bip = fmax(xi - aip / _mi, xi_max);
               double lip = _mi * (xi - bip) / aip;

               /* Set dx */
               dx[dof_it] += bip - xi;
               for (int j = 0; j < adj_dofs.Size(); j++)
               {
                  int dof_j = adj_dofs[j];
                  dx[dof_j] += lip * fmax(0., x_max[dof_j] - y[dof_j]);
               }
            }
         }
         else if (xi < xi_min)
         {
            num_min++;
            /* Compute ai- */
            double aim = 0.;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               aim += _mi_vec.Elem(dof_j)*fmax(0., y[dof_j] - x_min[dof_j]);
            }

            if (aim > 0.)
            {
               num_min_lim++;
               /* Compute bi- and li- */
               double bim = fmin(xi + aim / _mi, xi_min);
               double lim = _mi * (xi - bim) / aim;

               /* Set dx */
               dx[dof_it] += bim - xi;
               for (int j = 0; j < adj_dofs.Size(); j++)
               {
                  int dof_j = adj_dofs[j];
                  dx[dof_j] += lim * fmax(0., y[dof_j] - x_min[dof_j]);
               }
            }
         }
      } // End dof it

      /* End actual limit procedure */

      /* Check if any adjustment has been made */
      if (num_min_lim == num_min_lim_old && num_max_lim == num_max_lim_old)
      {
         /* No improvement from last iteratior */
         if (!suppress_output)
         {
            cout << "No improvement from last iteration.\n";
         }
         done_iterating = true;
      }

      // /* Check for max it reached */
      num_it++;
      if (num_it == max_it) { max_it_reached = true;}

      if (!suppress_output)
      {
         cout << "ratio of min able to be limited: " << double(num_min_lim)/num_min << ", num min: " << num_min << endl;
         cout << "ratio of max able to be limited: " << double(num_max_lim)/num_max << ", num max: " << num_max << endl;
      }

      /* Check conservation */
      y.Add(1., dx);
      bool is_cell_mass_conserved = CheckCellMassConservation(_mi_vec, x, y);
      bool is_globally_conservative = CheckGlobalMassConservation(_mi_vec, x, y);

      /* Output warning */
      if (!suppress_warnings)
      {
         if (!is_cell_mass_conserved)
         {
            MFEM_WARNING("Jacobi limiter not cellwise conservative.\n");
         }
         if (!is_globally_conservative)
         {
            MFEM_WARNING("Jacobi limiter not globally conservative.\n");
         }
      }
      x = y;
   } // IDP iteration
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
void LocalConservativeLimitGS(const Vector &_mi_vec, ParGridFunction &x)
{
   // cout << "IDPLimiter::LocalConservativeLimitGS\n";
   MFEM_ASSERT(local_bounds_computed, "Local bounds must be computed before applying local conservative limit.");

   int num_max_lim_old, num_min_lim_old;
   int num_max_lim = 0, num_min_lim = 0;
   int num_it = 0;
   bool max_it_reached = false, done_iterating = false;

   Vector x_old = x;

   /* Iterate max_it times or when there is no local budget for any of the limited quantities */
   while (!done_iterating && !max_it_reached)
   {
      num_min_lim_old = num_min_lim;
      num_max_lim_old = num_max_lim;

      num_min_lim = 0, num_max_lim = 0;
      int num_min = 0, num_max = 0;

      /* Actual limit procedure */
      Array<int> adj_dofs;

      for (int dof_it = 0; dof_it < NDofs; dof_it++)
      {
         real_t _mi = _mi_vec.Elem(dof_it);
         real_t x_max_i = x_max[dof_it], x_min_i = x_min[dof_it];
         real_t x_i = x[dof_it];
         // cout << "dof_it: " << dof_it << ", x_i: " << x_i << ", x_max_i: " << x_max_i << ", x_min_i: " << x_min_i << endl;
         // GetAdjacency(dof_it, adj_dofs);
         int el = L2_HO.GetElementForDof(dof_it);
         L2_HO.GetElementDofs(el, adj_dofs);
         
         if (_mi < 1.E-12) { 
            cout << "lumped mass is 0\n";
            x[dof_it] = fmin( fmax( x_i, x_min_i), x_max_i); 
         }
         else if (x_i > x_max_i)
         {
            num_max++;
            /* Compute ai+ */
            double aip = 0.;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               aip += _mi_vec.Elem(dof_j)*fmax(0., x_max[dof_j] - x[dof_j]);
            }

            if (aip > 0.)
            {
               num_max_lim++;
               /* Compute bi+ and li+*/
               double bip = fmax(x_i - aip / _mi, x_max_i);
               double lip = _mi * (x_i - bip) / aip;

               /* Check mass transfer is 0 */
               double mass_delta = _mi * (bip - x_i);

               /* Set gf_ho */
               x[dof_it] = bip;
               for (int j = 0; j < adj_dofs.Size(); j++)
               {
                  int dof_j = adj_dofs[j];
                  double _old_val = x[dof_j];
                  double _new_val = _old_val + lip * fmax(0., x_max[dof_j] - _old_val);
                  mass_delta += _mi_vec.Elem(dof_j) * (_new_val - _old_val);
                  x[dof_j] = _new_val;
               }
               if (abs(mass_delta) > 1.e-12)
               {
                  MFEM_ABORT("Mass delta should be 0.");
               }
            }
            
         }
         else if (x_i < x_min_i)
         {
            num_min++;

            /* Compute ai- */
            double aim = 0.;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               aim += _mi_vec.Elem(dof_j)*fmax(0., x[dof_j] - x_min[dof_j]);
            }

            if (aim > 0.)
            {
               num_min_lim++;

               /* Compute bi- and li- */
               double bim = fmin(x_i + aim / _mi, x_min_i);
               double lim = _mi * (x_i - bim) / aim;

               /* Check mass transfer is 0 */
               double mass_delta = _mi * (bim - x_i);

               /* Set gf_ho */
               x[dof_it] = bim;
               for (int j = 0; j < adj_dofs.Size(); j++)
               {
                  int dof_j = adj_dofs[j];
                  double _old_val = x[dof_j];
                  double _new_val = _old_val + lim * fmax(0., x[dof_j] - x_min[dof_j]);
                  mass_delta += _mi_vec.Elem(dof_j) * (_new_val - _old_val);
                  x[dof_j] = _new_val;
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
         if (!suppress_output)
         {
            cout << "No improvement from last iteration.\n";
         }
         done_iterating = true;
      }

      // /* Check for max it reached */
      num_it++;
      if (num_it == max_it) { max_it_reached = true;}

      if (!suppress_output)
      {
         cout << "ratio of min able to be limited: " << double(num_min_lim)/num_min << ", num min: " << num_min << endl;
         cout << "ratio of max able to be limited: " << double(num_max_lim)/num_max << ", num max: " << num_max << endl;
      }

      /* Check conservation */
      bool is_cell_mass_conserved = CheckCellMassConservation(_mi_vec, x_old, x);
      bool is_globally_conservative = CheckGlobalMassConservation(_mi_vec, x_old, x);

      /* Output warning */
      if (!suppress_warnings)
      {
         if (!is_cell_mass_conserved)
         {
            MFEM_WARNING("Gauss-Seidel limiter not cellwise conservative.\n");
         }
         if (!is_globally_conservative)
         {
            MFEM_WARNING("Gauss-Seidel limiter not globally conservative.\n");
         }
      }
   }
}


/**
 * @brief Applies a global conservative limiter to the given ParGridFunction.
 *
 * This function ensures that the values of the input ParGridFunction `gf_ho` 
 * are limited to a specified range [glob_x_min, glob_x_max] while 
 * preserving the total mass. The limiting process involves computing 
 * intermediate values and applying correction factors to maintain 
 * conservation properties.
 *
 * @param gf_ho The high-order ParGridFunction to be limited. It is modified
 *              in-place to satisfy the limiting constraints.
 *
 * @pre `gf_ho.Size()` must be equal to `NDofs`.
 * @pre `glob_x_min` must be non-negative, and `glob_x_max` must be 
 *      greater than or equal to `glob_x_min`.
 *
 * @post The values of `gf_ho` are within the range [glob_x_min, glob_x_max].
 * @post The total mass of `gf_ho` is conserved.
 *
 * @details
 * - The function first computes the total mass of the input function.
 * - It then calculates a limited version of the input values (`y`) by 
 *   clamping them to the range [x_min_i, x_max_i] for each degree of freedom.
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
 *         within the range [glob_x_min, glob_x_max].
 */
void GlobalConservativeLimit(const Vector &_mi_vec, ParGridFunction &gf_ho)
{
   // cout << "===== Limiter::GlobalConservativeLimit =====\n";
   double x_max_i, x_min_i;
   assert(gf_ho.Size() == NDofs);
   MFEM_ASSERT(local_bounds_computed, "Local bounds must be computed before applying global conservative limit.");
   if (use_global_bounds)
   {
      MFEM_ASSERT(glob_bounds_computed, "Global bounds must be computed before applying global conservative limit.");
   }
   
   if (!less_than_or_equal(0, glob_x_min) || !less_than_or_equal(glob_x_min, glob_x_max))
   {
      cout << "glob x min: " << glob_x_min << ", glob_x_max: " << glob_x_max << endl;
      MFEM_ABORT("x_min < 0 or x_max < x_min.");
   }

   /* Compute total mass */
   double M = 0.;
   for (int i = 0; i < NDofs; i++)
   {
      M += _mi_vec.Elem(i) * gf_ho[i]; 
   }
   if (!suppress_output)
   {
      cout << setprecision(15) << "M: " << M << endl;
   }

   /* Compute y */
   Vector y(NDofs); y = 0.;
   for (int i = 0; i < NDofs; i++)
   {
      if (use_global_bounds)
      {
         y[i] = fmin( fmax( gf_ho[i], glob_x_min), glob_x_max);
      }
      else
      {
         GetRhoMaxMin(i, x_max_i, x_min_i);
         y[i] = fmin( fmax( gf_ho[i], x_min_i), x_max_i);
      }
   }

   /* Compute alpha+ and alpha-*/
   double _sum_my = 0., _denom_max = 0., _denom_min = 0.;
   for (int i = 0; i < NDofs; i++)
   {
      double _m = _mi_vec.Elem(i);
      _sum_my += _m * y[i];

      if (use_global_bounds)
      {
         _denom_max += _m * (glob_x_max - y[i]);
         _denom_min += _m * (glob_x_min - y[i]);
      }
      else
      {
         GetRhoMaxMin(i, x_max_i, x_min_i);
         _denom_max += _m * (x_max_i - y[i]);
         _denom_min += _m * (x_min_i - y[i]);
      }
   }
   const double _num = M - _sum_my;
   const double alpha_p = fmax(0., _num / _denom_max);
   const double alpha_n = fmax(0., _num / _denom_min);

   MFEM_VERIFY(alpha_p >= 0. && alpha_p <= 1. && 
               alpha_n >= 0. && alpha_n <= 1.,
               "alpha_p and alpha_n should be between 0 and 1.\n");

   /* Finally, set z */
   Vector z(NDofs);
   for (int i = 0; i < NDofs; i++)
   {
      if (use_global_bounds)
      {
         z[i] = y[i] + alpha_p * (glob_x_max - y[i]) + alpha_n * (glob_x_min - y[i]);
      }
      else
      {
         GetRhoMaxMin(i, x_max_i, x_min_i);
         z[i] = y[i] + alpha_p * (x_max_i - y[i]) + alpha_n * (x_min_i - y[i]);
      }
   }

   /* Checks */
   bool is_globally_conservative = CheckGlobalMassConservation(_mi_vec, gf_ho, z);
   if (!is_globally_conservative && !suppress_warnings)
   {
      MFEM_WARNING("Global limiter not locally conservative.\n");
   }

   gf_ho = z;
   MFEM_VERIFY(less_than_or_equal(gf_ho.Max(), glob_x_max) && less_than_or_equal(glob_x_min, gf_ho.Min()), "Global limiting failed.\n");

   MFEM_WARNING("Global limiter is not locally mass conservative.");
   /* end checks */

}

bool less_than_or_equal(double a, double b) {
   // double eps = std::numeric_limits<double>::epsilon();
   double eps = 1.E-12;
   return a < b || std::abs(a - b) < eps * std::max({1.0, std::abs(a), std::abs(b)});
}

bool greater_than_or_equal(double a, double b) {
   double eps = std::numeric_limits<double>::epsilon();
   return a > b || std::abs(a - b) < eps * std::max({1.0, std::abs(a), std::abs(b)});
}

bool CheckEstimate(const Vector &_mi_vec, const ParGridFunction &gf_ho)
{
   // cout << "IDPLimiter::CheckEstimate\n";
   MFEM_ASSERT(local_bounds_computed, "Local bounds must be computed before checking estimate.");
   double M_max = 0., M_min = 0., M = 0.;
   for (int i = 0; i < NDofs; i++)
   {
      double _mi = _mi_vec[i];
      M += _mi * gf_ho[i];
      M_max += _mi * x_max[i];
      M_min += _mi * x_min[i];
   }
   if (less_than_or_equal(M_min, M) && less_than_or_equal(M, M_max))
   {
      return true;
   }
   else
   {
      if (less_than_or_equal(M, M_min)) { cout << "M < M_min\n"; }
      if (less_than_or_equal(M_max, M)) { cout << "M > M_max\n"; }
      cout << std::fixed << std::setprecision(20) << "M_min: " << M_min << ", M: " << M << ", M_max: " << M_max << endl;
      return false;
   }
}

void ComputeVolume(Vector &volumes)
{
   // cout << "IDPLimiter::ComputeVolume\n";
   volumes.SetSize(NDofs);
   ParLinearForm _lf(&L2_HO);
   ConstantCoefficient one(1.0);
   _lf.AddDomainIntegrator(new DomainLFIntegrator(one, &ir));
   _lf.Assemble();
   HypreParVector *_hpv = _lf.ParallelAssemble();
   volumes = *_hpv;
}


void Limit(const ParGridFunction &gf_lo, ParGridFunction &gf_ho)
{
   // cout << "===== IDPLimiter::Limit =====\n";
   ComputeLocalBounds(gf_lo);
   ComputeGlobalBounds();
   RelaxLocalBounds(gf_ho, x_min, x_max);
   Vector mi_vec;
   ComputeVolume(mi_vec);

   bool satisfies_estimate = CheckEstimate(mi_vec, gf_ho);
   if (!satisfies_estimate)
   {
      cout << "gf_lo: \n";
      gf_lo.Print(cout);
      MFEM_ABORT("Bounding estimate not satisfied");
   }

   switch(iter_method)
   {
      case JACOBI:
         LocalConservativeLimitJacobi(mi_vec, gf_ho);
         if (use_global_conservative_limiting)
         {
            GlobalConservativeLimit(mi_vec, gf_ho);
         }
         break;
      case GAUSS_SEIDEL:
         LocalConservativeLimitGS(mi_vec, gf_ho);
         if (use_global_conservative_limiting)
         {
            GlobalConservativeLimit(mi_vec, gf_ho);
         }
         break;
      default:
         MFEM_ABORT("Unknown iteration method.\n");
         break;
   }

   /* Reset flags */
   local_bounds_computed = false;
   glob_bounds_computed = false;
}

// Algorithm A.3 in Guermond, Wang 2025
void RelaxLocalBounds(const ParGridFunction &_gf_ho, ParGridFunction &_x_min, ParGridFunction &_x_max)
{
   switch (relaxation_option)
   {
      case 0:
         /* No relaxation */
         return;
      case 1:
      {
         /* Simple relaxation */
         RelaxBoundsMin(x_min, x_min_relaxed);
         RelaxBoundsMax(x_max, x_max_relaxed);

         x_min = x_min_relaxed;
         x_max = x_max_relaxed;

         return;
      }
      case 2:
         /* Stiffness-based relaxation */
         RelaxLocalBoundsStiffnessBased(_gf_ho, _x_min, _x_max);
         return;
      default:
         MFEM_ABORT("Unknown relaxation option.\n");
         break;
   }

}
void RelaxLocalBoundsStiffnessBased(const ParGridFunction &_gf_ho, ParGridFunction &_x_min, ParGridFunction &_x_max)
{
   MFEM_ASSERT(local_bounds_computed, "Local bounds must be computed before relaxing bounds.");
   MFEM_ASSERT(relaxation_option == 2, "Relaxation option must be 2 for stiffness based method.");
   MFEM_ASSERT(stiffness_mat_sp != nullptr, "Stiffness matrix must be assembled before relaxing bounds.");
   // If using global bounds, ensure that they have been computed
   if (use_global_bounds)
   {
      MFEM_ASSERT(glob_bounds_computed, "Global bounds must be computed before relaxing bounds with global bounds option.");
   }
   // cout << "IDPLimiter::RelaxLocalBoundsStiffnessBased\n";

   /* Compute alpha_i for all L2 dofs */
   Array<int> cols;
   Vector row_entries, alpha_vec(NDofs);
   for (int l2_dof_it = 0; l2_dof_it < NDofs; l2_dof_it++)
   {
      stiffness_mat_sp->GetRow(l2_dof_it, cols, row_entries);
      double alpha_i = 0.;
      /* Get ui */
      double ui = _gf_ho[l2_dof_it];
      // cout << "i: " << l2_dof_it << ", ui: " << ui << endl;
      for (int col_it = 0; col_it < cols.Size(); col_it++)
      {
         int j = cols[col_it];
         if (j == l2_dof_it) { continue; } // skip diagonal entry
         /* get uj */
         double uj = _gf_ho[j];
         double bij = row_entries[col_it];
         alpha_i += bij * (ui - uj);

         // cout << "j: " << cols[col_it] << ", uj: " << uj << ", beta_ij: " << bij << endl;
      }
      /* since bii = - \sum_{j\in I*(i)} bij*/
      // cout << "Denom: " << -row_entries[l2_dof_it] << endl;
      assert(cols[0] == l2_dof_it);
      alpha_i /= -row_entries[0];
      // cout << "alpha_i: " << alpha_i << endl;
      alpha_vec[l2_dof_it] = alpha_i;
   }

   // cout << "alpha_vec: \n";
   // alpha_vec.Print(cout);

   /* Set relaxed bounds */
   // Iterate over all dofs
   for (int l2_dof_it = 0; l2_dof_it < NDofs; l2_dof_it++)
   {
      double beta_i = alpha_vec[l2_dof_it];
      stiffness_mat_sp->GetRow(l2_dof_it, cols, row_entries);

      // cout << "i: " << l2_dof_it << ", alpha_i: " << alpha_vec[l2_dof_it] << endl;
      // cout << "neigbor alphas: ";

      // Iterate over neighbors to compute beta_i
      for (int col_it = 0; col_it < cols.Size(); col_it++)
      {
         if (cols[col_it] == l2_dof_it) { continue; } //
         double alpha_j = alpha_vec[cols[col_it]];

         // cout << alpha_j << ", ";

         // Opposite curvature, set relaxation to 0 to avoid oscillations
         if (beta_i * alpha_j <= 0.)
         {
            beta_i = 0.;
            break;
         }
         // Same curvature, take the smaller value
         else if (abs(beta_i) > abs(alpha_j))
         // if (abs(beta_i) > abs(alpha_j))
         {
            beta_i = alpha_j;
         }
      }
      // cout << endl;
      // cout << "dof: " << l2_dof_it << ", beta_i: " << beta_i << endl;
      if (fabs(beta_i) > 0.)
      {
         cout << "beta_i: " << beta_i << endl;
         _x_min[l2_dof_it] -= fabs(beta_i);
         _x_max[l2_dof_it] += fabs(beta_i);
         MFEM_ABORT("We have a nonzero beta_i");
      }

      // Ensure relaxed bounds do not exceed global bounds
      if (use_global_bounds)
      {
         _x_min[l2_dof_it] = fmax(_x_min[l2_dof_it], glob_x_min);
         _x_max[l2_dof_it] = fmin(_x_max[l2_dof_it], glob_x_max);
      }  
   }
}

void RelaxBoundsMin(const ParGridFunction &x_min_in, ParGridFunction &x_min_out)
{
   // cout << "IDPLimiter::RelaxBoundsMin\n";
   MFEM_ASSERT(local_bounds_computed, "Local bounds must be computed before relaxing bounds.");
   MFEM_ABORT("Check mass_vec");
   Array<int> adj_dofs;
   Vector x_min_2nd(NDofs);
   x_min_2nd = 0.;

   /* Compute 2nd order fromstencil*/
   for (int i = 0; i < NDofs; i++)
   {
      double x_min_i = x_min_in[i];
      GetAdjacency(i, adj_dofs);
      for (int j = 0; j < adj_dofs.Size(); j++)
      {
         int dof_j = adj_dofs[j];
         assert(dof_j != i);
         x_min_2nd[i] += x_min_i - x_min_in[dof_j];
      }
   }
   /* Compute avg 2nd order stencil and relaxed value */
   for (int i = 0; i < NDofs; i++)
   {
      const double x_min_i = x_min_in[i];
      const double x_min_2nd_i = x_min_2nd[i];
      GetAdjacency(i, adj_dofs);
      double _avg = 0.;

      for (int j = 0; j < adj_dofs.Size(); j++)
      {
         int dof_j = adj_dofs[j];
         _avg += (0.5*x_min_2nd_i + 0.5*x_min_2nd[dof_j]);
      }
      _avg /= 2. * adj_dofs.Size();
      
      double _rh = pow(mass_vec->Elem(i), 1.5 / dim);
      double _val1 = (1. - _rh) * x_min_i;
      double _val2 = x_min_i - abs(_avg);
      x_min_out[i] = fmin(_val1, _val2);
   }
}

void RelaxBoundsMax(const ParGridFunction &x_max_in, ParGridFunction &x_max_out)
{
   MFEM_ASSERT(local_bounds_computed, "Local bounds must be computed before relaxing bounds.");
   MFEM_ABORT("Check mass_vec");
   // cout << "IDPLimiter::RelaxBoundsMax\n";
   Array<int> adj_dofs;
   Vector x_max_2nd(NDofs);
   x_max_2nd = 0.;

   /* Compute 2nd order fromstencil*/
   for (int i = 0; i < NDofs; i++)
   {
      double x_max_i = x_max_in[i];
      GetAdjacency(i, adj_dofs);
      for (int j = 0; j < adj_dofs.Size(); j++)
      {
         int dof_j = adj_dofs[j];
         assert(dof_j != i);
         x_max_2nd[i] += x_max_i - x_max_in[dof_j];
      }
   }
   /* Compute avg 2nd order stencil and relaxed value */
   for (int i = 0; i < NDofs; i++)
   {
      const double x_max_i = x_max_in[i];
      const double x_max_2nd_i = x_max_2nd[i];
      GetAdjacency(i, adj_dofs);
      double _avg = 0.;

      for (int j = 0; j < adj_dofs.Size(); j++)
      {
         int dof_j = adj_dofs[j];
         _avg += (0.5*x_max_2nd_i + 0.5*x_max_2nd[dof_j]);
      }
      _avg /= 2. * adj_dofs.Size();


      double _rh = pow(mass_vec->Elem(i), 1.5 / dim);
      double _val1 = (1. + _rh) * x_max_i;
      double _val2 = x_max_i - abs(_avg);

      x_max_out[i] = fmax(_val1, _val2);
   }
}

/**
 * @brief Given the low order IDP solution, compute the local min and max values
 * for the high order solution. this computes u_i^min and u_i^max defined in 
 * equation 2.1.
 *
 * @param rho_gf_lo the invariant-domain-preserving low order update
 */
void ComputeLocalBounds(const ParGridFunction &gf_lo)
{
   // cout << "Limiter::ComputeLocalBounds\n";
   Array<int> adj_dofs;
   Vector sub_vec;
   for (int l2_dof_it = 0; l2_dof_it < NDofs; l2_dof_it++)
   {
      GetAdjacency(l2_dof_it, adj_dofs);
      gf_lo.GetSubVector(adj_dofs, sub_vec);
      double x_min_i = sub_vec.Min();
      double x_max_i = sub_vec.Max();
      x_min[l2_dof_it] = x_min_i;
      x_max[l2_dof_it] = x_max_i;
   }

   /* Set flag */
   local_bounds_computed = true;
}

void ComputeGlobalBounds()
{
   // cout << "Limiter::ComputeGlobalBounds\n";
   MFEM_ASSERT(local_bounds_computed, "Local bounds must be computed before computing global bounds.");

   glob_x_min = x_min.Min();
   glob_x_max = x_max.Max();
   glob_bounds_computed = true;
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
   Vector adj_dof_vals;
   el_adj_mat->GetRow(dof,adj_dofs,adj_dof_vals);
   adj_dofs.Sort();

   /* Next, append HO dofs that are in the same mesh element */
   // int el_index = L2_HO.GetElementForDof(dof);
   // Array<int> el_dofs;
   // L2_HO.GetElementDofs(el_index, el_dofs);
   // adj_dofs.Append(el_dofs);
   // adj_dofs = el_dofs;

   /* Remove dof from adj_dofs */
   // adj_dofs.DeleteFirst(dof);
   // adj_dofs.DeleteFirst(dof);
}

   
}; // IDPLimiter
} // ns mfem

#endif // LIMITER