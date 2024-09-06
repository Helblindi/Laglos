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
 */
template <int dim>
class LocalMassConservationOperator : public Operator
{
private:
   const ParGridFunction &X;
   const int num_cells, input_size;
   const int nnzSparse = 8;
   const Geometric<dim> &geom;
   double dt; 
   Array<int> GradCIArr, GradCJArr;
   Array<double> GradDataArr;
   SparseMatrix * grad;

public:
   LocalMassConservationOperator(const Geometric<dim> &_geom, const ParGridFunction &_X, 
                                 const int &_num_cells, const int &_input_size, const double &_dt,
                                 const Array<int> &GradCI, const Array<int> &GradCJ)
      : Operator(_num_cells, _input_size),
        num_cells(_num_cells),
        input_size(_input_size),
        X(_X),  
        geom(_geom),
        dt(_dt),
        GradCIArr(GradCI),
        GradCJArr(GradCJ)
   {
      // cout << "LocalMassConservationOperator::Non-default constructor\n";

      /* Initialize gradient */
      GradDataArr.SetSize(GradCJArr.Size());
      GradDataArr = 0.;
      grad = new SparseMatrix(GradCIArr.GetData(), GradCJArr.GetData(), GradDataArr.GetData(), num_cells, input_size);  
      // cout << "height: " << num_cells << ", width: " << input_size << endl;
   }
   ~LocalMassConservationOperator()
   {
      // delete grad; 
      // grad = nullptr;
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

      // cout << "nnz: " << nnzSparse << endl;

      for (int cell_it = 0; cell_it < num_cells; cell_it++)
      {
         /* Get adjacent vertices */
         geom.GetParMesh()->GetElementVertices(cell_it, verts);

         // cout << "cell: " << cell_it << ", verts: ";
         // verts.Print(cout);
         // cout << "----------\n";

         /* Get current node positions */
         geom.GetNodePosition(X, verts[0], x1);
         geom.GetNodePosition(X, verts[1], x2);
         geom.GetNodePosition(X, verts[2], x3);
         geom.GetNodePosition(X, verts[3], x4);

         // cout << "\tNode Positions:\n";
         // cout << "x1n: ";
         // x1.Print(cout);
         // cout << "x2n: ";
         // x2.Print(cout);
         // cout << "x3n: ";
         // x3.Print(cout);
         // cout << "x4n: ";
         // x4.Print(cout);
         // cout << "----------\n";

         /* Get node velocities */
         geom.GetNodeVelocityVecL(v, verts[0], v1n);
         geom.GetNodeVelocityVecL(v, verts[1], v2n);
         geom.GetNodeVelocityVecL(v, verts[2], v3n);
         geom.GetNodeVelocityVecL(v, verts[3], v4n);

         // cout << "\tNode Velocities:\n";
         // cout << "v1np1: ";
         // v1n.Print(cout);
         // cout << "v2np1: ";
         // v2n.Print(cout);
         // cout << "v3np1: ";
         // v3n.Print(cout);
         // cout << "v4np1: ";
         // v4n.Print(cout);
         // cout << "----------\n";

         /* Compute updated node locations */
         add(x1, dt, v1n, x1);
         add(x2, dt, v2n, x2);
         add(x3, dt, v3n, x3);
         add(x4, dt, v4n, x4);

         // cout << "\tMoved Node Positions:\n";
         // cout << "x1np1: ";
         // x1.Print(cout);
         // cout << "x2np1: ";
         // x2.Print(cout);
         // cout << "x3np1: ";
         // x3.Print(cout);
         // cout << "x4np1: ";
         // x4.Print(cout);
         // cout << "----------\n";

         /* Compute volume at time tn+1 */
         subtract(x3, x1, diag1);
         subtract(x2, x4, diag2);

         /* Rotate diagonals for gradient */
         geom.Perpendicular(diag1);
         geom.Perpendicular(diag2);

         /* Fill in gradient */
         int nnzSparseHalf = nnzSparse/2;
         Array<int> cols_arr(nnzSparse);
         for (int i = 0; i < nnzSparse/2; i++)
         {
            cols_arr[i] = verts[i];
            cols_arr[i+nnzSparse/2] = verts[i] + input_size / 2;
         }
         Array<int> rows_arr(1); rows_arr[0] = cell_it;

         /* Fill DM */
         DenseMatrix dm(1, nnzSparse);
         dm(0,0) = diag2[0]; 
         dm(0,nnzSparseHalf) = diag2[1];
         // V2
         dm(0,1) = diag1[0];
         dm(0,1+nnzSparseHalf) = diag1[1];
         // V3
         dm(0,2) = -diag2[0];
         dm(0,2+nnzSparseHalf) = -diag2[1];
         // V4
         dm(0,3) = -diag1[0];
         dm(0, 3+nnzSparseHalf) = -diag1[1];
         double coeff = dt / 2;
         dm *= coeff;

         /* Set row */
         // cout << "rows_arr: ";
         // rows_arr.Print(cout);
         // cout << "cols_arr: ";
         // cols_arr.Print(cout);

         // cout << "rows_arr: ";
         // rows_arr.Print(cout);
         // cout << "cols_arr: ";
         // cols_arr.Print(cout);
         // cout << "dt: " << dt << endl;
         // cout << "dm: ";
         // dm.Print(cout);
         grad->SetSubMatrix(rows_arr, cols_arr, dm);
         // cout << endl;
         
         /* Dense Matrix implementation */
         // V1
         // grad->Elem(cell_it, verts[0]) = coeff * diag2[0]; 
         // grad->Elem(cell_it, verts[0] + input_size / 2) = coeff * diag2[1];
         // // V2
         // grad->Elem(cell_it, verts[1]) = coeff * diag1[0];
         // grad->Elem(cell_it, verts[1] + input_size / 2) = coeff * diag1[1];
         // // V3
         // grad->Elem(cell_it, verts[2]) = - coeff * diag2[0];
         // grad->Elem(cell_it, verts[2] + input_size / 2) = - coeff * diag2[1];
         // // V4
         // grad->Elem(cell_it, verts[3]) = -coeff * diag1[0];
         // grad->Elem(cell_it, verts[3] + input_size / 2) = -coeff * diag1[1];
      }
      grad->Finalize();

      // MFEM_ABORT("End compute gradient.\n");
   }

   /* Evaluate the gradient operator at the point v. */
   virtual Operator &GetGradient(const Vector &v) const
   {
      // cout << "LMC::GetGradient\n";
      ComputeGradient(v);

      const int *In = grad->GetI(), *Jn = grad->GetJ();
      
      for (int i = 0, k=0; i < num_cells; i++)
      {
         // std::cout << "i: " << i << ", In[i]: " << In[i] << std::endl;
         // grad_C->GetRow(i, cols, row);
         // std::cout << "i: " << i << ", row: ";
         // row.Print(std::cout);
         // std::cout << "cols: ";
         // cols.Print(std::cout);
         // constr_grads->SetRow(i, cols, row);
         for (int end = In[i+1]; k < end; k++)
         {
            int j = Jn[k];
            // std::cout << "sparse entry [i: " << i << ", j: " << j << "]: " << grad->Elem(i,j) << std::endl;
            // constr_grads->Elem(i, j) = grad->Elem(i, j);
         }
      }

      // cout << "LMC::EndGradient\n";
      return *grad;
   }
};

template<int dim>
class zeroDenseMatrix : public Operator
{
private:
   DenseMatrix dm;
   SparseMatrix * grad;

public:
   zeroDenseMatrix(const int &height, const int &width) : 
      Operator(height, width),
      dm(height, width)
   {
      // cout << "zeroDenseMatrix constructor\n";
      // cout << "height: " << height << ", width: " << width << endl;
      grad = new SparseMatrix(height, width);

      Array<int> cols_arr(8);
      Array<int> rows_arr(1);
      for (int i = 0; i < 8; i++)
      {
         cols_arr[i] = i;
      }
      rows_arr[0] = 0;
      // Vector entries(8);
      DenseMatrix entries(1,8);
      entries = 1.;

      // grad->SetRow(0, cols_arr, entries);
      grad->SetSubMatrix(rows_arr, cols_arr, entries);

      grad->Finalize();
      dm = 0.;
   }
   ~zeroDenseMatrix()
   {
      delete grad;
      grad = nullptr;
   }

   virtual void Mult(const Vector &v, Vector &y) const
   {
      dm.Mult(v, y);
   }

   virtual Operator &GetGradient(const Vector &v) const 
   {
      if (dm.Height() == 1)
      {
         return *grad;
      }
      else 
      {
         MFEM_ABORT("Should not be using this operator class with hieght > 1.\n");
      }
      MFEM_ABORT("Should not have reached here.\n");
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
   const Geometric<dim> &geom;
   const int num_cells;
   const Vector &V_bar;
   Vector massvec, d_lo, d_hi;
   const LocalMassConservationOperator<dim> LMCoper;
   // const zeroDenseMatrix<dim> zDMoper;
   Array<int> HessIArr;
   Array<int> HessJArr;
   Array<double> HessData;
   SparseMatrix *hess;
   DenseMatrix block;

public:
   OptimizedMeshVelocityProblem(const Geometric<dim> &_geom, const Vector &_Vbar, const Vector &_massvec, 
                                const ParGridFunction &_X, const int _num_cells,
                                const double &dt, const Vector &xmin, const Vector &xmax,
                                const Array<int> &I, const Array<int> &J,
                                const Array<int> &GradCI, const Array<int> &GradCJ) 
      : geom(_geom),
        num_cells(_num_cells),
        OptimizationProblem(_Vbar.Size(), NULL, NULL),
        V_bar(_Vbar), massvec(_massvec), d_lo(1), d_hi(1),
        LMCoper(_geom, _X, _num_cells, input_size, dt, GradCI, GradCJ),
      //   zDMoper(1, input_size),
        HessIArr(I), HessJArr(J),
        block(4)
   {
      // cout << "OMVProblem constructor\n";
      // cout << "inputsize: " << input_size << endl;
      C = &LMCoper;
      SetEqualityConstraint(massvec);

      // Hiop must have C and D DenseMatrices
      // D = &zDMoper;
      // d_lo(0) = -1.E4;
      // d_hi(0) = 1.E4;
      // SetInequalityConstraint(d_lo, d_hi);

      SetSolutionBounds(xmin, xmax);

      HessData.SetSize(J.Size());
      HessData = 0.;
      hess = new SparseMatrix(HessIArr.GetData(), HessJArr.GetData(), HessData.GetData(), input_size, input_size); // input_size x input_size sparse matrix.

      block(0,1) = -0.5, block(0,3) = 0.5;
      block(1,0) = 0.5, block(1,2) = -0.5;
      block(2,1) = 0.5, block(2,3) = -0.5;
      block(3,0) = -0.5, block(3,2) = 0.5;

      block *= pow(dt, 2);
   }

   ~OptimizedMeshVelocityProblem()
   {
      // delete hess;
      // hess = nullptr;
   }

   virtual double CalcObjective(const Vector &V) const
   {
      // cout << "OMV calc objective\n";
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
      // cout << "OMV CalcObjectiveGrad\n";
      if (V.Size() != input_size)
      {
         MFEM_ABORT("Vectors must be of same size\n");
      }

      for (int i = 0; i < input_size; i++)
      {
         grad(i) = (V(i) - V_bar(i));
      }
      grad *= 2.;
   }

   void ComputeObjectiveHessData(const Vector &x) const
   {
      // std::cout << "OptimizedMeshVelocityProblem::ComputeObjectiveHessData\n";
      Array<int> verts, Vy_arr;
      DenseMatrix dm;

      for (int cell_it = 0; cell_it < num_cells; cell_it++)
      {
         /* Only add the contribution to the hessian if the LM is nonzero */
         double coeff = 0.5 * this->lambda[cell_it];
         if (coeff != 0.)
         {
            // std::cout << "cell_it: " << cell_it << std::endl;
            /* Get adjacent vertices */
            geom.GetParMesh()->GetElementVertices(cell_it, verts);

            Vy_arr.SetSize(verts.Size());
            for (int i = 0; i < verts.Size(); i++)
            {
               Vy_arr[i] = verts[i] + 0.5 * input_size;
            }

            // std::cout << "verts: ";
            // verts.Print(std::cout);
            // std::cout << "Vy_arr: ";
            // Vy_arr.Print(std::cout);

            // double val = hess->SearchRow(7,45);
            // val = hess->SearchRow(24,45);

            /* top right block */
            // hess->Elem(verts[0], Vy_arr[1]) = - coeff;
            // hess->Elem(verts[0], Vy_arr[3]) = coeff;
            // hess->Elem(verts[1], Vy_arr[0]) = coeff;
            // hess->Elem(verts[1], Vy_arr[2]) = - coeff;
            // hess->Elem(verts[2], Vy_arr[1]) = coeff;
            // hess->Elem(verts[2], Vy_arr[3]) = - coeff;
            // hess->Elem(verts[3], Vy_arr[0]) = - coeff;
            // hess->Elem(verts[3], Vy_arr[2]) = coeff;
            dm = block;
            dm *= this->lambda[cell_it];
            // std::cout << "lambda: " << this->lambda[cell_it] << std::endl;
            // std::cout << "dm to be added to sparse matrix: ";
            // dm.Print(std::cout);
            DenseMatrix dm2;
            // hess->GetSubMatrix(verts, Vy_arr, dm2);
            // std::cout << "dm pre add: ";
            hess->AddSubMatrix(verts, Vy_arr, dm, 2/*skip_zeros*/);

            // /* bottom left block */
            // cout << "dm^T to be added: ";
            dm.Transpose();
            // dm.Print(cout);
            hess->AddSubMatrix(Vy_arr, verts, dm, 2/*skip_zeros*/);

            // /* Check if mat was added correctly */
            // hess->GetSubMatrix(verts, Vy_arr, dm);
            // std::cout << "getsubmatrix: ";
            // dm.Print(std::cout);
         }
         
      }

      /* Don't forget the diagonal entries */
      for (int i = 0; i < input_size; i++)
      {
         hess->Elem(i,i) = 2.;
      }

      hess->Finalize();

      // std::cout << "OptimizedMeshVelocityProblem::ComputeObjectiveHessData - DONE\n";
   }

   virtual Operator & CalcObjectiveHess(const Vector &x) const 
   {
      // std::cout << "OptimizedMeshVelocityProblem::CalcObjectiveHessian\n";
      ComputeObjectiveHessData(x);

      // std::cout << "OptimizedMeshVelocityProblem::CalcObjectiveHessian - DONE\n";
      return *hess;
   }
};

} // End ns hydrodynamics

} // End ns mfem

#endif // LAGRANGE_MULTIPLIER