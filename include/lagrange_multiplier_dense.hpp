#ifndef LAGRANGE_MULTIPLIER_DENSE
#define LAGRANGE_MULTIPLIER_DENSE

#include "mfem.hpp"
#include "geometry.hpp"

using namespace std;

namespace mfem
{

namespace hydroLO
{
/**
 * Class to represent the local mass conservation
 * constraints on our Lagrange Multiplier problem
 */
class LocalMassConservationOperatorDense : public Operator
{
private:
   const int dim;
   const ParGridFunction &X;
   const int num_cells, input_size;
   const Geometric &geom;
   double dt; 
   mutable DenseMatrix grad;

public:
   LocalMassConservationOperatorDense(const int &_dim, const Geometric &_geom, const ParGridFunction &_X, 
                                 const int &_num_cells, const int &_input_size, const double &_dt)
      : Operator(_num_cells, _input_size),
        dim(_dim),
        num_cells(_num_cells),
        input_size(_input_size),
        grad(_num_cells, _input_size),
        X(_X),  
        geom(_geom),
        dt(_dt)
   {
      // cout << "LocalMassConservationOperatorDense::Non-default constructor\n";
      // cout << "input_size: " << input_size << endl;
      /* Initialize gradient */
      grad = 0.;
   }

   virtual void Mult(const Vector &v, Vector &y) const
   {
      // cout << "LMC::Mult()\n";
      Vector x1(dim), x2(dim), x3(dim), x4(dim);
      Vector v1n(dim), v2n(dim), v3n(dim), v4n(dim);
      Vector diag1(dim), diag2(dim), temp_vec(dim);

      Array<int> verts;

      for (int cell_it = 0; cell_it < num_cells; cell_it++)
      {
         /* Get adjacent vertices */
         geom.GetParMesh()->GetElementVertices(cell_it, verts);

         /* Get current node positions */
         geom.GetNodePosition(X, verts[0], x1);
         geom.GetNodePosition(X, verts[1], x2);
         geom.GetNodePosition(X, verts[2], x3);
         geom.GetNodePosition(X, verts[3], x4);

         /* Get node velocities */
         geom.GetNodeVelocityVecL(v, verts[0], v1n);
         geom.GetNodeVelocityVecL(v, verts[1], v2n);
         geom.GetNodeVelocityVecL(v, verts[2], v3n);
         geom.GetNodeVelocityVecL(v, verts[3], v4n);

         /* Compute updated node locations */
         add(x1, dt, v1n, x1);
         add(x2, dt, v2n, x2);
         add(x3, dt, v3n, x3);
         add(x4, dt, v4n, x4);

         /* Compute volume at time tn+1 */
         subtract(x3, x1, diag1);
         subtract(x2, x4, diag2);
         geom.Perpendicular(diag1);
         double val = diag1 * diag2;

         y[cell_it] = val;
         // if (val <= 0.)
         // {
         //    MFEM_ABORT("Cell volume should be positive.\n");
         // }
      }
      y *= 0.5;
   }

   void ComputeGradient(const Vector &v) const
   {
      // cout << "LMC::ComputeGradient\n";
      Vector x1(dim), x2(dim), x3(dim), x4(dim);
      Vector v1n(dim), v2n(dim), v3n(dim), v4n(dim);
      Vector diag1(dim), diag2(dim), temp_vec(dim);

      Array<int> verts;

      for (int cell_it = 0; cell_it < num_cells; cell_it++)
      {
         /* Get adjacent vertices */
         geom.GetParMesh()->GetElementVertices(cell_it, verts);

         /* Get current node positions */
         geom.GetNodePosition(X, verts[0], x1);
         geom.GetNodePosition(X, verts[1], x2);
         geom.GetNodePosition(X, verts[2], x3);
         geom.GetNodePosition(X, verts[3], x4);

         /* Get node velocities */
         geom.GetNodeVelocityVecL(v, verts[0], v1n);
         geom.GetNodeVelocityVecL(v, verts[1], v2n);
         geom.GetNodeVelocityVecL(v, verts[2], v3n);
         geom.GetNodeVelocityVecL(v, verts[3], v4n);

         /* Compute updated node locations */
         add(x1, dt, v1n, x1);
         add(x2, dt, v2n, x2);
         add(x3, dt, v3n, x3);
         add(x4, dt, v4n, x4);

         /* Compute volume at time tn+1 */
         subtract(x3, x1, diag1);
         subtract(x2, x4, diag2);

         /* Rotate diagonals for gradient */
         geom.Perpendicular(diag1);
         geom.Perpendicular(diag2);

         /* Fill in gradient */
         double coeff = dt / 2;
         // V1
         grad.Elem(cell_it, verts[0]) = coeff * diag2[0]; 
         grad.Elem(cell_it, verts[0] + input_size / 2) = coeff * diag2[1];
         // V2
         grad.Elem(cell_it, verts[1]) = coeff * diag1[0];
         grad.Elem(cell_it, verts[1] + input_size / 2) = coeff * diag1[1];
         // V3
         grad.Elem(cell_it, verts[2]) = - coeff * diag2[0];
         grad.Elem(cell_it, verts[2] + input_size / 2) = - coeff * diag2[1];
         // V4
         grad.Elem(cell_it, verts[3]) = -coeff * diag1[0];
         grad.Elem(cell_it, verts[3] + input_size / 2) = -coeff * diag1[1];
      }
   }

   /* Evaluate the gradient operator at the point v. */
   virtual Operator &GetGradient(const Vector &v) const
   {
      // cout << "LMC::GetGradient\n";
      ComputeGradient(v);

      return grad;
   }
};

class zeroDenseMatrixDense : public Operator
{
private:
   const int dim;
   DenseMatrix dm;
   mutable DenseMatrix grad;

public:
   zeroDenseMatrixDense(const int &_dim, const int &height, const int &width) : 
      Operator(height, width),
      dim(_dim),
      dm(height, width),
      grad(dim, width)
   {
      grad = 0.;
      dm = 0.;
   }
   zeroDenseMatrixDense(const int &_dim, const int &n) :
      Operator(n,n),
      dim(_dim),
      dm(n),
      grad(n,n)
   {
      grad = 0.;
      dm = 0.;
   }

   virtual void Mult(const Vector &v, Vector &y) const
   {
      dm.Mult(v, y);
   }

   virtual Operator &GetGradient(const Vector &v) const 
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
class OptimizedMeshVelocityProblemDense : public OptimizationProblem
{
private:
   const int dim;
   const Vector &V_target;
   Vector massvec, d_lo, d_hi;
   const LocalMassConservationOperatorDense LMCoper;
   zeroDenseMatrixDense zDMoper;

public:
   OptimizedMeshVelocityProblemDense(const int &_dim, const Geometric &_geom, const Vector &_V_target, const Vector &_massvec, 
                                const ParGridFunction &_X, const int num_cells,
                                const double &dt, const Vector &xmin, const Vector &xmax) 
      : OptimizationProblem(_V_target.Size(), NULL, NULL),
        dim(_dim),
        LMCoper(dim, _geom, _X, num_cells, input_size, dt),
        zDMoper(dim, 1, input_size),
      //   OptimizationProblem(_V_target.Size(), &LMCoper, &zDMoper),
        V_target(_V_target), massvec(_massvec), d_lo(1), d_hi(1)
   {
      cout << "OMVProblem constructor\n";

      C = &LMCoper;
      SetEqualityConstraint(massvec);

      // Hiop must have C and D DenseMatrices
      D = &zDMoper;
      d_lo(0) = -1.E4;
      d_hi(0) = 1.E4;
      SetInequalityConstraint(d_lo, d_hi);

      SetSolutionBounds(xmin, xmax);
   }

   virtual double CalcObjective(const Vector &V) const
   {
      // cout << "OMV calc objective\n";
      double res = 0.0;
      for (int i = 0; i < input_size; i++)
      {
         const double d = V(i) - V_target(i);
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
      // cout << "OMV CalcObjectiveGrad\n";
      if (V.Size() != input_size)
      {
         MFEM_ABORT("Vectors must be of same size\n");
      }

      for (int i = 0; i < input_size; i++)
      {
         grad(i) = (V(i) - V_target(i));
      }
      grad *= 2.;
   }
};

} // End ns hydroLO

} // End ns mfem

#endif // LAGRANGE_MULTIPLIER_DENSE