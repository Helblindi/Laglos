#ifndef LAGRANGE_MULTIPLIER
#define LAGRANGE_MULTIPLIER

#include "mfem.hpp"
#include "geometry.hpp"

using namespace std;

namespace mfem
{

namespace hydrodynamics
{

/**
 * Class to represent the local mass conservation
 * constraints on our Lagrange Multiplier problem
 * C(V)
 */
template <int dim>
class LocalMassConservationOperator : public Operator
{
private:
   const ParGridFunction &X;
   const ParMesh *pmesh;
   const int num_cells;
   const Geometric<dim> &geom;
   double dt; // TODO: Add in dt here to compute updated nodal locations
   mutable DenseMatrix grad;


public:
   LocalMassConservationOperator(const Geometric<dim> &_geom, const ParGridFunction &_X, 
                                 const ParMesh &_pmesh, const int _num_cells, const double &_dt)
      : Operator(_num_cells, _X.Size()),
        X(_X), pmesh(&_pmesh), 
        num_cells(_num_cells),
        geom(_geom),
        dt(_dt),
        grad(_num_cells, _X.Size())
   {
      cout << "LMCOp constructor\n";
      cout << "op width: " << _X.Size() << endl;
      cout << "LMC num rows: " << this->NumRows() << endl;
      cout << "dt: " << dt << endl;

      /* Initialize gradient */
      grad = 0.;

      /* 
      * Set sparsite pattern for gradient matrix 
      * Hiop only allows for DenseMatrix grad object???
      */
      // grad = new SparseMatrix(num_cells, _X.Size(), 8); // Hardcoded dim * num_corners = 8
   }

   virtual void Mult(const Vector &x, Vector &y) const
   {
      cout << "LMC::Mult()\n";
      Vector x1(dim), x2(dim), x3(dim), x4(dim);
      Vector v1n(dim), v2n(dim), v3n(dim), v4n(dim);
      Vector diag1(dim), diag2(dim), temp_vec(dim);

      Array<int> verts;

      for (int cell_it = 0; cell_it < num_cells; cell_it++)
      {
         /* Get adjacent vertices */
         pmesh->GetElementVertices(cell_it, verts);

         /* Get current node positions */
         geom.GetNodePosition(X, verts[0], x1);
         geom.GetNodePosition(X, verts[1], x2);
         geom.GetNodePosition(X, verts[2], x3);
         geom.GetNodePosition(X, verts[3], x4);

         /* Get node velocities */
         geom.GetNodeVelocity(x, verts[0], v1n);
         geom.GetNodeVelocity(x, verts[1], v2n);
         geom.GetNodeVelocity(x, verts[2], v3n);
         geom.GetNodeVelocity(x, verts[3], v4n);

         /* Compute updated node locations */
         add(x1, dt, v1n, x1);
         add(x2, dt, v2n, x2);
         add(x3, dt, v3n, x3);
         add(x4, dt, v4n, x4);

         /* Compute volume at time tn+1 */
         subtract(x3, x1, diag1);
         subtract(x2, x4, diag2);
         geom.Orthogonal(diag1);
         double val = diag1 * diag2;

         y[cell_it] = val;

         /* Rotate diagonals for gradient */
         // geom.Orthogonal(diag1);
         geom.Orthogonal(diag2);

         /* Fill in gradient */
         double coeff = dt / 2;
         grad.Elem(cell_it, verts[0]) = -coeff * diag1[1]; 
         grad.Elem(cell_it, verts[0] + X.Size() / 2) = coeff * diag1[0];
         grad.Elem(cell_it, verts[1]) = coeff * diag2[1];
         grad.Elem(cell_it, verts[1] + X.Size() / 2) = - coeff * diag2[0];
         grad.Elem(cell_it, verts[2]) = coeff * diag1[1];
         grad.Elem(cell_it, verts[2] + X.Size() / 2) = - coeff * diag1[0];
         grad.Elem(cell_it, verts[3]) = -coeff * diag2[1];
         grad.Elem(cell_it, verts[3] + X.Size() / 2) = coeff * diag2[0];
      }
      y *= 0.5;

   }

   /* Evaluate the gradient operator at the point x. */
   virtual Operator &GetGradient(const Vector &x) const
   {
      // TODO: Fix here
      cout << "LMVGetGradient: \n";
      // grad.Print(cout);
      return grad;
   }
};

template<int dim>
class zeroDenseMatrix : public Operator
{
private:
   DenseMatrix dm;
   mutable DenseMatrix grad;

public:
   zeroDenseMatrix(const int &height, const int &width) : 
      Operator(height, width),
      dm(height, width),
      grad(dim, width)
   {
      grad = 0.;
      dm = 0.;
   }
   zeroDenseMatrix(const int &n) :
      Operator(n,n),
      dm(n),
      grad(n,n)
   {
      grad = 0.;
      dm = 0.;
   }

   virtual void Mult(const Vector &x, Vector &y) const
   {
      dm.Mult(x, y);
   }

   virtual Operator &GetGradient(const Vector &x) const 
   {
      if (dm.Height() == 1)
      {
         return grad;
      }
      else 
      {
         MFEM_ABORT("Should not be using this operator class with hieght > 1.\n");
      }
   } 
};

/**
 * Conservative a posteriori correction to mesh velocity to guarantee local mass conservation:
 * Find V that minimizes || Vi - \bar{V}i ||^2, subject to
 *   Kc(n+1)   Kcn
 *   ------- - --- = 0  for all cells c in the mesh,
 *   Tc(n+1)   Tcn
 * 
 * OR (as implemented here)
 *                     Kcn
 *   Kc(n+1) = Tc(n+1) ---
 *                     Tcn 
 *
 * Here, \bar{V}_i is a weighted average of the adjacent cell velocities.
 * 
 * Note that both velocity vectors are of size dim * NVDofsH1
 */
template <int dim>
class OptimizedMeshVelocityProblem : public OptimizationProblem
{
private:
   const Vector &V_bar;
   Vector massvec, d_lo, d_hi;
   const LocalMassConservationOperator<dim> LMCoper;
   zeroDenseMatrix<dim> zDMoper;

public:
   OptimizedMeshVelocityProblem(const Geometric<dim> &_geom, const Vector &_Vbar, const Vector &_massvec, 
                                const ParGridFunction &_X, const ParMesh &_pmesh, const int num_cells,
                                const double &dt, const Vector &xmin, const Vector &xmax) 
      : OptimizationProblem(_Vbar.Size(), NULL, NULL),
        LMCoper(_geom, _X, _pmesh, num_cells, dt),
        zDMoper(1, input_size),
        V_bar(_Vbar), massvec(_massvec), d_lo(1), d_hi(1)
   {
      cout << "input size: " << input_size << endl;
      cout << "OMVProblem constructor\n";
      cout << "num_cells: " << num_cells << endl;
      C = &LMCoper;
      SetEqualityConstraint(massvec);
      // Hiop must have C and D DenseMatrices

      D = &zDMoper;
      d_lo(0) = -1.E4;
      d_hi(0) = 1.E4;
      SetInequalityConstraint(d_lo, d_hi);

      cout << "xmin: ";
      xmin.Print(cout);
      cout << "xmax: ";
      xmax.Print(cout);
      // MFEM_ABORT(false);

      SetSolutionBounds(xmin, xmax);
   }

   virtual double CalcObjective(const Vector &V) const
   {
      cout << "OMV calc objective\n";
      cout << "inpute size: " << input_size << endl;
      double res = 0.0;
      for (int i = 0; i < input_size; i++)
      {
         const double d = V(i) - V_bar(i);
         res += d * d;
      }
      return res;
   }
   int get_input_size() 
   {
      return input_size;
   }

   virtual void CalcObjectiveGrad(const Vector &V, Vector &grad) const 
   {
      cout << "OMV CalcObjectiveGrad\n";
      if (V.Size() != input_size)
      {
         MFEM_ABORT("Vectors must be of same size\n");
      }

      for (int i = 0; i < input_size; i++)
      {
         grad(i) = (V(i) - V_bar(i));
      }
      grad *= 2.;
      cout << "End OMV CalcObjectiveGrad\n";
   }
};

} // End ns hydrodynamics

} // End ns mfem

#endif // LAGRANGE_MULTIPLIER