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
   ParFiniteElementSpace &L2;
   Vector *mass_vec;
   const int NDofs;
   ParGridFunction mass_lumped;
   ParGridFunction vol_vec;
   mutable ParGridFunction rho_min, rho_max;
   mutable ParGridFunction rho_min_relaxed, rho_max_relaxed;
   double glob_rho_max = -1., glob_rho_min = -1.;

   ParFiniteElementSpace H1_LO, H1_HO;
   ParGridFunction LO_proj_max, LO_proj_min;
   ParDiscreteLinearOperator vert_elem_oper; 
   HypreParMatrix *vert_elem, *vert_elem_trans;
   SparseMatrix _spm, _spm_trans;

   /* Misc options */
   const int dim;
   const int max_it = 2;
   bool use_glob = false;
   const int relaxation_option = 1; // 0: no relaxation, 1: relaxation type 1, 2: relaxation type 2
   bool suppress_warnings = false;
   IterationMethod iter_method = GAUSS_SEIDEL; // options: JACOBI, GAUSS_SEIDEL

public:
IDPLimiter(ParFiniteElementSpace & l2, ParFiniteElementSpace & H1FESpace_proj_LO,
           ParFiniteElementSpace & H1FESpace_proj_HO, Vector &_mass_vec) :
   L2(l2),
   mass_vec(&_mass_vec),
   NDofs(mass_vec->Size()),
   mass_lumped(&l2),
   vol_vec(&l2),
   rho_min(&l2),
   rho_max(&l2),
   rho_min_relaxed(&l2),
   rho_max_relaxed(&l2),
   H1_LO(H1FESpace_proj_LO),
   H1_HO(H1FESpace_proj_HO),
   LO_proj_max(&H1_LO),
   LO_proj_min(&H1_LO),
   vert_elem_oper(&H1_HO, &l2),
   dim(L2.GetParMesh()->Dimension())
{
   mass_lumped = 0.;
   vol_vec = 0.;
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

   cout << "IDPLimiter configuration options:\n";
   cout << "  - max_it: " << max_it << endl;
   cout << "  - relaxation_option: " << relaxation_option << endl;
   cout << "  - suppress_warnings: " << suppress_warnings << endl;
}

~IDPLimiter() 
{
   // delete spm;
   delete vert_elem;
}

void GetVolVec(Vector &vol_vec_out) const { vol_vec_out = vol_vec;}
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
      sum = mass_lumped[i] * (y[i] - x[i]);
      GetAdjacency(i, adj_dofs);
      // cout << "adj dofs: ";
      // adj_dofs.Print(cout);
      for (int adj_it = 0; adj_it < adj_dofs.Size(); adj_it++)
      {
         int j = adj_dofs[adj_it];
         sum += mass_lumped[j] * (y[j] - x[j]);
      }

      if (abs(sum) > 1.e-12)
      {
         cout << "local mass not conserved at dof: " << i << ", sum: " << sum << endl;
         // cout << "adj dofs: ";
         // adj_dofs.Print(cout);
         is_locally_conservative = false;
         // cout << "Not locally conservative at dof: " << i << ", val: " << sum << endl;
         break; // Break out of loop as there is not point looking at other dofs
      }
   }
   return is_locally_conservative;
}

bool CheckLocalCellMassConservation(const Vector &x_old, const Vector &x_new)
{
   MFEM_ABORT("Function requires updated volume to work.");
   bool is_cell_mass_conserved = true;
   for (int i = 0; i < L2.GetNE(); i++)
   {
      double old_mass = 0., new_mass = 0.;
      Array<int> dofs;
      L2.GetElementDofs(i, dofs);
      for (int dof_it = 0; dof_it < dofs.Size(); dof_it++)
      {
         int dof = dofs[dof_it];
         old_mass += mass_lumped[dof] * x_old[dof];
         new_mass += mass_lumped[dof] * x_new[dof];
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

bool CheckGlobalMassConservation(const Vector &x, const Vector &y) const
{
   // cout << "IDPLimiter::CheckGlobalMassConservation\n";
   bool is_globally_conservative = true;
   double sum_old = 0., sum_new = 0.;
   for (int i = 0; i < NDofs; i++)
   {
      /* Check if mass has been preserved */
      sum_old += mass_lumped[i] * x[i];
      sum_new += mass_lumped[i] * y[i];
   }
   double t_val = abs(sum_old - sum_new);
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
void LocalConservativeLimitJacobi(ParGridFunction &x)
{
   cout << "IDPLimiter::LocalConservativeLimitJacobi\n";
   int num_it = 0;
   bool max_it_reached = false, done_iterating = false;
   while (!done_iterating && !max_it_reached)
   {
      /* Actual limit procedure */
      Array<int> adj_dofs;
      Vector y(NDofs), dx(NDofs);
      y = x;
      dx = 0.;
      for (int dof_it = 0; dof_it < NDofs; dof_it++)
      {
         double mi = mass_lumped[dof_it];
         double xi_max = rho_max[dof_it], xi_min = rho_min[dof_it];
         double xi = y[dof_it];
         GetAdjacency(dof_it, adj_dofs);

         if (mi == 0.) { y[dof_it] = min( max( xi, xi_min), xi_max); }
         else if (xi > xi_max)
         {
            /* Compute ai+ */
            double aip = 0.;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               aip += mass_lumped[dof_j]*max(0., rho_max[dof_j] - y[dof_j]);
            }

            /* Compute bi+ */
            double bip = max(xi - aip / mi, xi_max);

            /* Compute li+ */
            double lip = 0.;
            if (aip > 0.) { lip = mi * (xi - bip) / aip; }
            else if (aip < 0.) { MFEM_ABORT("aim <= 0."); }

            /* Set dx */
            dx[dof_it] += bip - xi;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               dx[dof_j] += lip * max(0., rho_max[dof_j] - y[dof_j]);
            }
         }
         else if (xi < xi_min)
         {
            /* Compute ai- */
            double aim = 0.;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               aim += mass_lumped[dof_j]*std::max(0., y[dof_j] - rho_min[dof_j]);
            }

            /* Compute bi+ */
            double bim = min(xi + aim / mi, xi_min);

            /* Compute li+ */
            double lim = 0.;
            if (aim > 0.)
            {
               lim = mi * (xi - bim) / aim;
            }
            else if (aim < 0.)
            {
               MFEM_ABORT("aim <= 0.");
            }

            /* Set y */
            dx[dof_it] += bim - xi;
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               dx[dof_j] += lim * std::max(0., y[dof_j] - rho_min[dof_j]);
            }
         }
      } // End dof it

      /* End actual limit procedure */

      // /* Check for max it reached */
      num_it++;
      if (num_it == max_it) { max_it_reached = true;}

      // cout << "num min limited: " << num_min_lim << endl;
      // cout << "num max limited: " << num_max_lim << endl;

      /* Check conservation */
      y.Add(1., dx);
      // cout << "checking local mass conservation\n";
      bool is_locally_conservative = CheckLocalMassConservation(x, y);
      if (!is_locally_conservative && !suppress_warnings)
      {
         MFEM_WARNING("Gauss-Seidel limiter not locally conservative.\n");
      }
      // cout << "checkout local cell mass conservation\n";
      // bool is_cell_mass_conserved = CheckLocalCellMassConservation(x, y);
      // if (!is_cell_mass_conserved && !suppress_warnings)
      // {
      //    MFEM_WARNING("Gauss-Seidel limiter not locally conservative on cell.\n");
      // }
      // cout << "checking global mass conservation\n";
      bool is_globally_conservative = CheckGlobalMassConservation(x, y);
      if (!is_globally_conservative && !suppress_warnings)
      {
         MFEM_WARNING("Gauss-Seidel limiter not globally conservative.\n");
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
void LocalConservativeLimitGS(const ParGridFunction &x, ParGridFunction &gf_ho)
{
   cout << "IDPLimiter::LocalConservativeLimit\n";

   cout << "gf_ho:\n";
   gf_ho.Print(cout);

   int num_max_lim_old, num_min_lim_old;
   int num_max_lim = 0, num_min_lim = 0;
   int num_it = 0;
   bool max_it_reached = false, done_iterating = false;

   // const IntegrationRule &ir = IntRules.Get(L2.GetFE(0)->GetGeomType(), 5);
   /* e vec to quadrature */
   const FiniteElement *fe = gf_ho.ParFESpace()->GetFE(0);
   auto *Tr = x.FESpace()->GetMesh()->GetElementTransformation(0);
   const IntegrationRule &ir = MassIntegrator::GetRule(*fe, *fe, *Tr);
   auto qi_u = L2.GetQuadratureInterpolator(ir);

   const int nqp = ir.GetNPoints();
   const int NE = L2.GetNE();

   /* Interpolate the gf to be limited */
   Vector u_qvals(nqp * NE);
   qi_u->Values(gf_ho, u_qvals);

   /* Use as our volume object */
   GeometricFactors geom(x, ir, GeometricFactors::DETERMINANTS);

   /* Interpolate the bounds*/
   Vector rho_max_qvals(nqp * NE), rho_min_qvals(nqp * NE);
   qi_u->Values(rho_max, rho_max_qvals);
   qi_u->Values(rho_min, rho_min_qvals);

   Vector x_old = u_qvals;
   Vector gf_ho_old = gf_ho;

   // cout << "\t---gf_ho first:---\n";
   // gf_ho.Print(cout);
   // cout << "\t end gh_ho\n";
   // cout << "\t===== uqvals first: =====\n";
   // u_qvals.Print(cout);
   // cout << "\t===== end uqvals =====\n";

   /* Iterate max_it times or when there is no local budget for any of the limited quantities */
   // while (!done_iterating && !max_it_reached)
   // {
   //    num_min_lim_old = num_min_lim;
   //    num_max_lim_old = num_max_lim;

   //    num_min_lim = 0, num_max_lim = 0;
   //    int num_min = 0, num_max = 0;

   //    /* Actual limit procedure */
   //    Array<int> adj_dofs;

   //    for (int dof_it = 0; dof_it < nqp * NE; dof_it++)
   //    {
   //       double mi = geom.detJ(dof_it);
   //       double rho_max_i = rho_max_qvals[dof_it], rho_min_i = rho_min_qvals[dof_it];
   //       double rho_i = u_qvals[dof_it];
   //       // cout << "dof_it: " << dof_it << ", rho_i: " << rho_i << ", rho_max_i: " << rho_max_i << ", rho_min_i: " << rho_min_i << endl;
   //       // GetAdjacency(dof_it, adj_dofs);
   //       adj_dofs.SetSize(nqp-1);
   //       int el = dof_it / nqp;
   //       // cout << "el: " << el << ", dof_it: " << dof_it << endl;
   //       int q = dof_it % nqp;
   //       for (int i = 0; i < nqp; i++)
   //       {
   //          if (i != q)
   //          {
   //             adj_dofs[i] = i + el * nqp;
   //          }
   //       }
   //       // cout << "qp: " << dof_it << ", adj dofs: ";
   //       // adj_dofs.Print(cout);
   //       // cout << std::setprecision (15) << "i: " << dof_it << ", adj_dofs: ";
   //       // adj_dofs.Print(cout);
   //       if (mi == 0.) { u_qvals[dof_it] = min( max( rho_i, rho_min_i), rho_max_i); }
   //       else if (rho_i > rho_max_i)
   //       {
   //          num_max++;
   //          /* Compute ai+ */
   //          double aip = 0.;
   //          for (int j = 0; j < adj_dofs.Size(); j++)
   //          {
   //             int dof_j = adj_dofs[j];
   //             aip += geom.detJ(dof_j)*max(0., rho_max_qvals[dof_j] - u_qvals[dof_j]);
   //          }

   //          if (aip > 0.)
   //          {
   //             // cout << "aip > 0\n";
   //             num_max_lim++;
   //             /* Compute bi+ and li+*/
   //             double bip = max(rho_i - aip / mi, rho_max_i);
   //             double lip = mi * (rho_i - bip) / aip;

   //             // cout << "aip: " << aip << ", bip: " << bip << ", lip: " << lip << endl;

   //             /* Check mass transfer is 0 */
   //             double mass_delta = mi * (bip - rho_i);

   //             /* Set gf_ho */
   //             u_qvals[dof_it] = bip;
   //             for (int j = 0; j < adj_dofs.Size(); j++)
   //             {
   //                int dof_j = adj_dofs[j];
   //                double _old_val = u_qvals[dof_j];
   //                double _new_val = _old_val + lip * max(0., rho_max_qvals[dof_j] - _old_val);
   //                mass_delta += geom.detJ(dof_j) * (_new_val - _old_val);
   //                u_qvals[dof_j] = _new_val;
   //             }
   //             if (abs(mass_delta) > 1.e-12)
   //             {
   //                MFEM_ABORT("Mass delta should be 0.");
   //             }
   //          }
            
   //       }
   //       else if (rho_i < rho_min_i)
   //       {
   //          num_min++;
   //          /* Compute ai- */
   //          double aim = 0.;
   //          for (int j = 0; j < adj_dofs.Size(); j++)
   //          {
   //             int dof_j = adj_dofs[j];
   //             aim += geom.detJ(dof_j)*max(0., u_qvals[dof_j] - rho_min_qvals[dof_j]);
   //          }

   //          if (aim > 0.)
   //          {
   //             // cout << "aim > 0\n";
   //             num_min_lim++;

   //             /* Compute bi- and li- */
   //             double bim = min(rho_i + aim / mi, rho_min_i);
   //             double lim = mi * (rho_i - bim) / aim;

   //             // cout << setprecision(12) << "\trho_i: " << rho_i << ", rho_min_i: " << rho_min_i << endl;
   //             // cout << "aim: " << aim << ", bim: " << bim << ", lim: " << lim << endl;
               

   //             /* Check mass transfer is 0 */
   //             double mass_delta = mi * (bim - rho_i);

   //             /* Set gf_ho */
   //             // cout << "i. old val: " << gf_ho[dof_it] << ", new val: " << bim << ", rho_min: " << rho_min[dof_it] << endl;
   //             u_qvals[dof_it] = bim;
   //             for (int j = 0; j < adj_dofs.Size(); j++)
   //             {
   //                int dof_j = adj_dofs[j];
   //                double _old_val = u_qvals[dof_j];
   //                double _new_val = _old_val + lim * max(0., u_qvals[dof_j] - rho_min_qvals[dof_j]);
   //                mass_delta += geom.detJ(dof_j) * (_new_val - _old_val);
   //                u_qvals[dof_j] = _new_val;
   //             }

   //             if (abs(mass_delta) > 1.e-12)
   //             {
   //                MFEM_ABORT("Mass delta should be 0.");
   //             }
   //          }
   //       }
   //    }

   //    /* End actual limit procedure */
   //    cout << "ratio of min able to be limited: " << double(num_min_lim)/num_min << ", num min: " << num_min << endl;
   //    cout << "ratio of max able to be limited: " << double(num_max_lim)/num_max << ", num max: " << num_max << endl;

   //    /* Check if any adjustment has been made */
   //    if (num_min_lim == num_min_lim_old && num_max_lim == num_max_lim_old)
   //    {
   //       /* No improvement from last iteratior */
   //       cout << "No improvement from last iteration.\n";
   //       done_iterating = true;
   //    }

   //    // /* Check for max it reached */
   //    num_it++;
   //    if (num_it == max_it) { max_it_reached = true;}

   //    /* Check conservation */
   //    // cout << "checking local mass conservation\n";
   //    // bool is_locally_conservative = CheckLocalMassConservation(x_old, u_qvals);
   //    // if (!is_locally_conservative && !suppress_warnings)
   //    // {
   //    //    MFEM_WARNING("Gauss-Seidel limiter not locally conservative.\n");
   //    // }
   //    // cout << "checkout local cell mass conservation\n";
   //    // bool is_cell_mass_conserved = CheckLocalCellMassConservation(x_old, gf_ho);
   //    // if (!is_cell_mass_conserved && !suppress_warnings)
   //    // {
   //    //    MFEM_WARNING("Gauss-Seidel limiter not locally conservative on cell.\n");
   //    // }
   //    // cout << "checking global mass conservation\n";
   //    // bool is_globally_conservative = CheckGlobalMassConservation(x_old, gf_ho);
   //    // if (!is_globally_conservative && !suppress_warnings)
   //    // {
   //    //    MFEM_WARNING("Gauss-Seidel limiter not globally conservative.\n");
   //    // }
   // }

   /* Check different from old to new vals */
   Vector _tvec(NE*nqp);
   _tvec = u_qvals;
   _tvec.Add(-1., x_old);
   cout << "change in u_qvals: " << _tvec.Sum() << endl;

   /* project back onto gf_ho */
   QuadratureSpace qspace(*L2.GetMesh(), ir);
   QuadratureFunction qf(qspace);
   qf = u_qvals;
   mfem::QuadratureFunctionCoefficient u_qvals_coeff(qf);
   gf_ho.ProjectCoefficient(u_qvals_coeff);

   _tvec.SetSize(gf_ho.Size());
   _tvec = gf_ho;
   _tvec.Add(-1., gf_ho_old);
   cout << "_tvec: \n";
   _tvec.Print(cout);
   cout << "change in gf_ho: " << _tvec.Sum() << endl;


   cout << "gf_ho:\n";
   gf_ho.Print(cout);

   // cout << "\t---gf_ho last:---\n";
   // gf_ho.Print(cout);
   // cout << "\t end gh_ho\n";
   // cout << "\t===== uqvals last: =====\n";
   // u_qvals.Print(cout);
   // cout << "\t===== end uqvals =====\n";
   cout << "IDPLimiter::LocalConservativeLimitGS done.\n";
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
      M += mass_lumped[i] * gf_ho[i]; 
   }
   cout << "M: " << M << endl;

   /* Compute y */
   Vector y(NDofs); y = 0.;
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
      double _m = mass_lumped[i];
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
   bool is_globally_conservative = CheckGlobalMassConservation(gf_ho, z);
   if (!is_globally_conservative && !suppress_warnings)
   {
      MFEM_WARNING("Global limiter not locally conservative.\n");
   }

   MFEM_ASSERT(gf_ho.Max() <= glob_rho_max && gf_ho.Min() >= glob_rho_min,
               "Global limiting failed.\n");
   /* end checks */

   gf_ho = z;
}

bool CheckEstimate(const Vector &_lumped_mass_vec, const ParGridFunction &gf_ho)
{
   // cout << "IDPLimiter::CheckEstimate\n";
   double M_max = 0., M_min = 0., M = 0., mi;
   for (int i = 0; i < NDofs; i++)
   {
      mi = _lumped_mass_vec[i];
      M += mi * gf_ho[i];
      M_max += mi * rho_max[i];
      M_min += mi * rho_min[i];
   }
   if (M >= M_min && M <= M_max)
   {
      return true;
   }
   else
   {
      cout << "M_min: " << M_min << ", M: " << M << ", M_max: " << M_max << endl;
      return false;
   }
}

void ComputeVolume(const ParGridFunction &gf_ho)
{
   cout << "IDPLimiter::ComputeVolume\n";
   // MFEM_ABORT("Need to compute volume based on quadrature");
   vol_vec = 0.;
   for (int i = 0; i < NDofs; i++)
   {
      vol_vec[i] = mass_vec->Elem(i) / gf_ho[i];
   }
}

void ComputeLumpedMasses(const ParGridFunction &gf_ho)
{
   LinearForm _lf(&L2);
   ConstantCoefficient one(1.0);
   _lf.AddDomainIntegrator(new DomainLFIntegrator(one));
   _lf.Assemble();

   mass_lumped = _lf;
}

void ComputeCellMassesAndVolumes(const ParGridFunction &gf_ho)
{
   // cout << "IDPLimiter::ComputeCellMasses\n";
   // MFEM_ABORT("Need to compute volume based on quadrature");
   const int NE = L2.GetNE();
   Vector cell_mass(NE);
   Vector cell_vol(NE);
   for (int i = 0; i < NE; i++)
   {
      cell_mass[i] = 0.;
      cell_vol[i] = 0.;
      Array<int> dofs;
      L2.GetElementDofs(i, dofs);
      for (int j = 0; j < dofs.Size(); j++)
      {
         int dof_j = dofs[j];
         cell_mass[i] += mass_lumped[dof_j] * gf_ho[dof_j];
         cell_vol[i] += mass_lumped[dof_j];
      }
   }
   cout << "cell masses: ";
   cell_mass.Print(cout);
   cout << "sum idp cell masses: " << cell_mass.Sum() << endl;
   cout << "cell volumes: ";
   cell_vol.Print(cout);
   cout << "sum idp cell volumes: " << cell_vol.Sum() << endl;
}

void ComputeCellMassesAndVolumesQUAD(const ParGridFunction &x, const ParGridFunction &gf_ho)
{
   cout << "IDPLimiter::ComputeCellMassesQUAD\n";
   // MFEM_ABORT("Need to compute volume based on quadrature");
   const IntegrationRule &ir = IntRules.Get(L2.GetFE(0)->GetGeomType(), 5);
   auto qi_u = L2.GetQuadratureInterpolator(ir);
   const int nqp = ir.GetNPoints();
   const int NE = L2.GetNE();
   cout << "nqp: " << nqp << ", NE: " << NE << endl;
   Vector u_qvals(nqp * NE);
   // mass_lumped_qvals(nqp * NE);
   qi_u->Values(gf_ho, u_qvals);
   // qi_u->Values(mass_lumped, mass_lumped_qvals);

   GeometricFactors geom(x, ir, GeometricFactors::DETERMINANTS);

   Vector el_mass(NE), el_vol(NE);

   for (int i = 0; i < NE; i++)
   {
      el_mass[i] = 0.;
      el_vol[i] = 0.;
      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         int dof_it = i * nqp + q;
         el_mass[i] +=  ip.weight * geom.detJ(dof_it) * u_qvals[dof_it];
         el_vol[i] += ip.weight * geom.detJ(dof_it);
      }

   }
   cout << "cell masses: ";
   el_mass.Print(cout);
   cout << "sum idp cell masses: " << el_mass.Sum() << endl;
   cout << "cell volumes: ";
   el_vol.Print(cout);
   cout << "sum idp cell volumes: " << el_vol.Sum() << endl;
}

void Limit(const ParGridFunction &x, const ParGridFunction &gf_lo, ParGridFunction &gf_ho)
{
   cout << "===== Limiter::Limit =====\n";
   ComputeRhoMinMax(gf_lo);
   // ComputeVolume(gf_ho);
   ComputeLumpedMasses(gf_ho);
   cout << "---Initial cell masses---\n";
   // ComputeCellMassesAndVolumes(gf_ho);
   ComputeCellMassesAndVolumesQUAD(x, gf_ho);
   // assert(false);

   bool satisfies_estimate = CheckEstimate(mass_lumped, gf_ho);
   if (!satisfies_estimate)
   {
      MFEM_ABORT("Bounding estimate not satisfied");
   }

   switch(iter_method)
   {
      case JACOBI:
         LocalConservativeLimitJacobi(gf_ho);
         // GlobalConservativeLimit(gf_ho);
         break;
      case GAUSS_SEIDEL:
         LocalConservativeLimitGS(x, gf_ho);
         // GlobalConservativeLimit(gf_ho);
         break;
      default:
         MFEM_ABORT("Unknown iteration method.\n");
         break;
   }

   cout << "---Final cell masses QUAD---\n";
   // ComputeCellMassesAndVolumes(gf_ho);
   ComputeCellMassesAndVolumesQUAD(x, gf_ho);

   cout << "===== END Limiter::Limit =====\n";
   assert(false);
}

void RelaxBoundsMin(const ParGridFunction &rho_min_in, ParGridFunction &rho_min_out)
{
   // cout << "IDPLimiter::RelaxBoundsMin\n";
   Array<int> adj_dofs;
   Vector rho_min_2nd(NDofs);
   rho_min_2nd = 0.;

   /* Compute 2nd order fromstencil*/
   for (int i = 0; i < NDofs; i++)
   {
      double rho_min_i = rho_min_in[i];
      GetAdjacency(i, adj_dofs);
      for (int j = 0; j < adj_dofs.Size(); j++)
      {
         int dof_j = adj_dofs[j];
         assert(dof_j != i);
         rho_min_2nd[i] += rho_min_i - rho_min_in[dof_j];
      }
   }
   /* Compute avg 2nd order stencil and relaxed value */
   for (int i = 0; i < NDofs; i++)
   {
      const double rho_min_i = rho_min_in[i];
      const double rho_min_2nd_i = rho_min_2nd[i];
      GetAdjacency(i, adj_dofs);
      double _avg = 0.;
      switch(relaxation_option)
      {
         case 1:
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               _avg += (0.5*rho_min_2nd_i + 0.5*rho_min_2nd[dof_j]);
            }
            _avg /= 2. * adj_dofs.Size();
            break;
         case 2:
            MFEM_ABORT("Not implemented.");
            break;
         default:
            MFEM_ABORT("Unknown relaxation option.\n");
            break;
      }
      
      double _rh = pow(mass_lumped[i], 1.5 / dim);
      double _val1 = (1. - _rh) * rho_min_i;
      double _val2 = rho_min_i - abs(_avg);
      // cout << "_val1: " << _val1 << ", _val2: " << _val2 << endl;
      rho_min_out[i] = min(_val1, _val2);
   }
}

void RelaxBoundsMax(const ParGridFunction &rho_max_in, ParGridFunction &rho_max_out)
{
   // cout << "IDPLimiter::RelaxBoundsMax\n";
   Array<int> adj_dofs;
   Vector rho_max_2nd(NDofs);
   rho_max_2nd = 0.;

   /* Compute 2nd order fromstencil*/
   for (int i = 0; i < NDofs; i++)
   {
      double rho_max_i = rho_max_in[i];
      GetAdjacency(i, adj_dofs);
      for (int j = 0; j < adj_dofs.Size(); j++)
      {
         int dof_j = adj_dofs[j];
         assert(dof_j != i);
         rho_max_2nd[i] += rho_max_i - rho_max_in[dof_j];
      }
   }
   /* Compute avg 2nd order stencil and relaxed value */
   for (int i = 0; i < NDofs; i++)
   {
      const double rho_max_i = rho_max_in[i];
      const double rho_max_2nd_i = rho_max_2nd[i];
      GetAdjacency(i, adj_dofs);
      double _avg = 0.;
      switch(relaxation_option)
      {
         case 1:
            for (int j = 0; j < adj_dofs.Size(); j++)
            {
               int dof_j = adj_dofs[j];
               _avg += (0.5*rho_max_2nd_i + 0.5*rho_max_2nd[dof_j]);
            }
            _avg /= 2. * adj_dofs.Size();
            break;
         case 2:
            MFEM_ABORT("Not implemented.");
            break;
         default:
            MFEM_ABORT("Unknown relaxation option.\n");
            break;
      }
      
      double _rh = pow(mass_lumped[i], 1.5 / dim);
      double _val1 = (1. + _rh) * rho_max_i;
      double _val2 = rho_max_i - abs(_avg);
      // cout << "_val1: " << _val1 << ", _val2: " << _val2 << endl;
      rho_max_out[i] = max(_val1, _val2);
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
   
   /* Optionally, relax the bounds */
   if (relaxation_option > 0)
   {
      RelaxBoundsMin(rho_min, rho_min_relaxed);
      RelaxBoundsMax(rho_max, rho_max_relaxed);

      /* compare bounds to relaxed bounds */
      // for (int i = 0; i < NDofs; i++)
      // {
      //    cout << "i: " << i << ", rho_min: " << rho_min[i] << ", rho_min_relaxed: " << rho_min_relaxed[i] << endl;
      //    cout << "rho_max: " << rho_max[i] << ", rho_max_relaxed: " << rho_max_relaxed[i] << endl;
      // }

      rho_min = rho_min_relaxed;
      rho_max = rho_max_relaxed;
   }

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
   int el_index = L2.GetElementForDof(dof);
   Array<int> el_dofs;
   L2.GetElementDofs(el_index, el_dofs);
   // adj_dofs.Append(el_dofs);
   adj_dofs = el_dofs;

   /* Remove dof from adj_dofs */
   adj_dofs.DeleteFirst(dof);
   // adj_dofs.DeleteFirst(dof);
}

   
}; // IDPLimiter
} // ns mfem

#endif // LIMITER