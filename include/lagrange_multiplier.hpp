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
 * Helper functions to set the sparsity pattern that pertains to this implementation
 */
/****************************************************************************************************
* Function: SetHiopConstraintGradSparsityPattern
*
* Purpose:
*  Set the sparsity pattern for the Jacobian matrix of the constraints vector to be used in the 
*  Hiop mesh velocity solve.  Recall the the Gradient matrix is a SparseMatrix
*  of size |eta_cel| x 2*eta_geo
****************************************************************************************************/
inline void SetHiopConstraintGradSparsityPattern(const ParMesh *pmesh, const int &num_cells, const int &num_vertices, Array<int> &I, Array<int> &J)
{
   // TODO: Remove this function in LagrangianLOOperator
   Array<int> verts;

   I.SetSize(num_cells + 1);
   J.SetSize(8*num_cells); /// TODO: Remove hardcoded parameter representing the number of nonzeros per row
   I[0] = 0;

   for (int cell_it = 0; cell_it < num_cells; cell_it++)
   {
      /* Set the (i+1)th entry of I */
      I[cell_it+1] = 8*(cell_it+1); /// TODO: Remove hardcoded parameter representing the number of nonzeros per row

      /* Fill J */
      pmesh->GetElementVertices(cell_it, verts);
      verts.Sort();
      int verts_length = verts.Size();
      assert(verts_length == 4);
      for (int j = 0; j < verts_length; j++)
      {
         // cout << "J at " << 8*cell_it + j << ": " << verts[j] << endl;
         J[8*cell_it+j] = verts[j];
         J[8*cell_it +j+4] = verts[j] + num_vertices;
         // cout << "J at " << 8*cell_it + j +4<< ": " << verts[j] + num_vertices << endl;
      }
   }
}


/****************************************************************************************************
* Function: SetHiopHessianSparsityPattern
*
* Purpose:
*  Set the sparsity pattern for the Hessian matrix to be used in the 
*  Hiop mesh velocity solve.  Recall the the Hessian matrix is a SparseMatrix
*  of size 2*eta_geo x 2*eta_geo
*
*  Note: According to hiop/src/Linalg/README.md, this matrix is upper triangular.
*
*        1x  2x  3x ... 1y 2y ...
*    1x  *   0
*    2x  0   *            *
*    ..           *  ...
*    1y       0          * 
*    2y                     *
*    ..
*
*  Upper triangular, CSR.
****************************************************************************************************/

inline void SetHiopHessianSparsityPattern(const ParMesh *pmesh, const ParFiniteElementSpace &H1, const int &num_vertices, Array<int> &I, Array<int> &J)
{
   // std::cout << "SetHiopHessianSparsityPattern\n";
   Array<int> faces_row, face_dofs_row;
   I.SetSize(num_vertices+1);
   J.SetSize(0);
   int faces_length;
   mfem::Mesh::FaceInformation FI;

   Table *edge_vertex(pmesh->GetEdgeVertexTable());
   Table vertex_edge;
   Transpose(*edge_vertex, vertex_edge);
   
   for (int node = 0; node < num_vertices; node++)
   {
      vertex_edge.GetRow(node, faces_row);
      faces_length = faces_row.Size();
      Array<int> cols_for_node(faces_length+1);

      /* Fill I */
      if (node <= num_vertices) {
         I[node+1] = I[node] + faces_length + 1; 
      }

      for (int face_it = 0; face_it < faces_length; face_it++)
      {
         /* Get adj node index*/
         int face = faces_row[face_it];
         FI = pmesh->GetFaceInformation(face);
         H1.GetFaceDofs(face, face_dofs_row);
         int adj_node = (node == face_dofs_row[0]) ? face_dofs_row[1] : face_dofs_row[0];

         /* Fill adjacency array */
         cols_for_node[face_it] = adj_node + num_vertices;
      }
      /* Append diagonal element and sort */
      cols_for_node[faces_length] = node;
      cols_for_node.Sort();

      /* Append onto J */
      J.Append(cols_for_node);
   }

   /* 
   * Fill rest of I and J for y components. Upper Triangular 
   * Notice that the nonzero entries for the y components of the
   * velocity will be lower triangular and thus are omitted.
   * */
   const int base = I[num_vertices];
   for (int i = 0; i < num_vertices; i++)
   {
      I.Append(base + i + 1);
      J.Append(num_vertices + i);
   }

   /* Validate propert of csr matrices */
   assert(I[2*num_vertices] == J.Size());
}


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
      // cout << "grad height: " << grad->Height() << ", grad width: " << grad->Width() << endl;
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

         grad->SetSubMatrix(rows_arr, cols_arr, dm);
      }
      grad->Finalize();
   }

   /* Evaluate the gradient operator at the point v. */
   virtual Operator &GetGradient(const Vector &v) const
   {
      // cout << "LMC::GetGradient\n";
      ComputeGradient(v);

      const int *In = grad->GetI(), *Jn = grad->GetJ();
      
      for (int i = 0, k=0; i < num_cells; i++)
      {
         for (int end = In[i+1]; k < end; k++)
         {
            int j = Jn[k];
         }
      }

      return *grad;
   }
};

template<int dim>
class zeroSparseMatrix : public Operator
{
private:
   Array<int> arrI, arrJ;
   Array<double> arrData;
   SparseMatrix * grad;

public:
   zeroSparseMatrix(const int &height, const int &width) : 
      Operator(height, width),
      arrI(2), arrJ(8), arrData(8)
   {
      // cout << "zeroSparseMatrix constructor\n";

      arrI[0] = 0, arrI[1] = 8;

      for (int i = 0; i < 8; i++)
      {
         arrJ[i] = i;
         arrData[i] = 0.;
      }

      grad = new SparseMatrix(arrI.GetData(), arrJ.GetData(), arrData.GetData(), height, width);  

      grad->Finalize();
   }
   // ~zeroSparseMatrix()
   // {
   //    delete grad;
   //    grad = nullptr;
   // }

   virtual void Mult(const Vector &v, Vector &y) const
   {
      y = 0.;
   }

   virtual Operator &GetGradient(const Vector &v) const 
   {
      return *grad;
   } 
};

/**
 * Conservative a posteriori correction to mesh velocity to guarantee local mass conservation:
 * Find V that minimizes || Vi - {V_target}_i ||^2, subject to
 *   Kc(n+1)   Kcn
 *   ------- - --- = 0  for all cells c in the mesh,
 *   Tc(n+1)   Tcn
 * 
 * OR (as implemented here)
 *                     Kcn
 *   Kc(n+1) = Tc(n+1) ---
 *                     Tcn 
 *
 * Here, {V_target}_i is some target velocity that approximates the mesh velocity well.
 * That is to say that moving the mesh with V_target results in minimal necessary correction 
 * at the faces to guarantee local mass conservation.
 * 
 * Note that both velocity vectors are of size dim * NVDofsH1
 */
template <int dim>
class TargetOptimizedMeshVelocityProblem : public OptimizationProblem
{
private:
   const Geometric<dim> &geom;
   const int num_cells, nnz_sparse_jaceq, nnz_sparse_Hess_Lagr;
   const Vector &V_target;
   Vector massvec, d_lo, d_hi;
   const LocalMassConservationOperator<dim> LMCoper;
   // zeroSparseMatrix<dim> zSMoper;
   Array<int> HessIArr;
   Array<int> HessJArr;
   Array<double> HessData;
   SparseMatrix *hess;
   DenseMatrix block;

public:
   TargetOptimizedMeshVelocityProblem(
      const Geometric<dim> &_geom, const Vector &_V_target, const Vector &_massvec, 
      const ParGridFunction &_X, const int _num_cells,
      const double &dt, const Vector &xmin, const Vector &xmax,
      const Array<int> &I, const Array<int> &J,
      const Array<int> &GradCI, const Array<int> &GradCJ) 
      : geom(_geom),
        num_cells(_num_cells),
        OptimizationProblem(dim*_geom.GetNVDofs_H1(), NULL, NULL),
        V_target(_V_target), massvec(_massvec), d_lo(1), d_hi(1),
        LMCoper(_geom, _X, _num_cells, input_size, dt, GradCI, GradCJ),
        nnz_sparse_jaceq(GradCJ.Size()),
        nnz_sparse_Hess_Lagr(J.Size()),
      //   zSMoper(1, input_size),
        HessIArr(I), HessJArr(J),
        block(4)
   {
      // cout << "TOMVProblem constructor\n";

      // Problem assumes that vectors are of the correct size
      assert(_massvec.Size() == _num_cells);
      assert(_V_target.Size() == input_size);

      C = &LMCoper;
      SetEqualityConstraint(massvec);

      // D = &zSMoper;
      // d_lo(0) = -1.E4;
      // d_hi(0) = 1.E4;
      // SetInequalityConstraint(d_lo, d_hi);

      SetSolutionBounds(xmin, xmax);

      HessData.SetSize(J.Size());
      HessData = 0.;
      hess = new SparseMatrix(HessIArr.GetData(), HessJArr.GetData(), HessData.GetData(), input_size, input_size); // input_size x input_size sparse matrix.
      hess->Finalize();

      block(0,1) = -0.5, block(0,3) = 0.5;
      block(1,0) = 0.5, block(1,2) = -0.5;
      block(2,1) = 0.5, block(2,3) = -0.5;
      block(3,0) = -0.5, block(3,2) = 0.5;

      block *= pow(dt, 2);
   }

   ~TargetOptimizedMeshVelocityProblem()
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
      assert(grad.Size() == input_size);
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

   void ComputeObjectiveHessData(const Vector &x) const
   {
      // std::cout << "OptimizedMeshVelocityProblem::ComputeObjectiveHessData\n";
      Array<int> verts, Vy_arr;
      DenseMatrix dm;

      /* Reset data in sparse Hessian since we build it through addition */
      double * _data = hess->GetData();
      for (int d = 0; d < HessJArr.Size(); d++)
      {
         _data[d] = 0.;
      }

      /* Verify our Hessian is starting out zeroed out */
      if (hess->MaxNorm() > 0.)
      {
         std::cout << "max norm: " << hess->MaxNorm() << std::endl;
         MFEM_ABORT("Hessian has not been reset properly.\n")
      }

      /* Fill data array cell by cell */
      for (int cell_it = 0; cell_it < num_cells; cell_it++)
      {
         /* Only add the contribution to the hessian if the LM is nonzero */
         double coeff = 0.5 * this->lambda[cell_it];
         if (coeff != 0.)
         {
            /* Get adjacent vertices */
            geom.GetParMesh()->GetElementVertices(cell_it, verts);

            Vy_arr.SetSize(verts.Size());
            for (int i = 0; i < verts.Size(); i++)
            {
               Vy_arr[i] = verts[i] + 0.5 * input_size;
            }

            /* top right block */
            dm = block;
            dm *= this->lambda[cell_it];

            hess->AddSubMatrix(verts, Vy_arr, dm, 2/*skip_zeros*/);
         }
      } // End cell it

      /* Set the diagonal entries */
      for (int i = 0; i < input_size; i++)
      {
         hess->Elem(i,i) = 2.;
      }
   }

   virtual Operator & CalcObjectiveHess(const Vector &x) const 
   {
      // std::cout << "OptimizedMeshVelocityProblem::CalcObjectiveHessian\n";
      ComputeObjectiveHessData(x);

      // std::cout << "OptimizedMeshVelocityProblem::CalcObjectiveHessian - DONE\n";
      return *hess;
   }

   virtual int get_nnz_sparse_Jaceq() const { return nnz_sparse_jaceq; }
   virtual int get_nnz_sparse_Jacineq() const { return 0; }
   virtual int get_nnz_sparse_Hess_Lagr() const { return nnz_sparse_Hess_Lagr; }
};


/**
 * Conservative a posteriori correction to mesh velocity to guarantee local mass conservation:
 * Find V that minimizes SUM_i SUM_{j\in nbr{(i)} || Vi - Vj ||^2, subject to
 *   Kc(n+1)   Kcn
 *   ------- - --- = 0  for all cells c in the mesh,
 *   Tc(n+1)   Tcn
 *
 * OR (as implemented here)
 *                     Kcn
 *   Kc(n+1) = Tc(n+1) ---
 *                     Tcn
 * 
 * Note that both velocity vectors are of size dim * NVDofsH1
 */
template <int dim>
class ViscousOptimizedMeshVelocityProblem : public OptimizationProblem
{
private:
   const Geometric<dim> &geom;
   const int num_cells;
   const int nnz_sparse_jaceq, nnz_sparse_Hess_Lagr;
   Vector massvec, d_lo, d_hi;
   const LocalMassConservationOperator<dim> LMCoper;
   // const zeroSparseMatrix<dim> zSMoper;
   Array<int> HessIArr;
   Array<int> HessJArr;
   Array<double> HessData;
   SparseMatrix *hess;
   DenseMatrix block;

public:
   ViscousOptimizedMeshVelocityProblem(const Geometric<dim> &_geom, const Vector &_massvec,
                                       const ParGridFunction &_X, const int _num_cells,
                                       const double &dt, const Vector &xmin, const Vector &xmax,
                                       const Array<int> &I, const Array<int> &J,
                                       const Array<int> &GradCI, const Array<int> &GradCJ)
      : geom(_geom),
        num_cells(_num_cells),
        OptimizationProblem(_geom.GetNVDofs_H1(), NULL, NULL),
        massvec(_massvec), d_lo(1), d_hi(1),
        LMCoper(_geom, _X, _num_cells, input_size, dt, GradCI, GradCJ),
        nnz_sparse_jaceq(GradCJ.Size()),
        nnz_sparse_Hess_Lagr(J.Size()),
      //   zSMoper(1, input_size),
        HessIArr(I), HessJArr(J),
        block(4)
   {
      MFEM_ABORT("OptimizationProblem is not yet implemented.\n");
   }

   ~ViscousOptimizedMeshVelocityProblem()
   {
      // delete hess;
      // hess = nullptr;
   }

   virtual double CalcObjective(const Vector &V) const
   {
      MFEM_ABORT("CalcObjective must be implemented.\n");
   }
   int get_input_size()
   {
      return input_size;
   }

   virtual void CalcObjectiveGrad(const Vector &V, Vector &grad) const
   {
      if (V.Size() != input_size)
      {
         MFEM_ABORT("Vectors must be of same size\n");
      }
      MFEM_ABORT("CalcObjectiveGrad must be implemented.\n");
   }

   void ComputeObjectiveHessData(const Vector &x) const
   {
   }

   virtual Operator & CalcObjectiveHess(const Vector &x) const
   {
      ComputeObjectiveHessData(x);

      return *hess;
   }
   virtual int get_nnz_sparse_Jaceq() const { return nnz_sparse_jaceq; }
   virtual int get_nnz_sparse_Jacineq() const { return 0; }
   virtual int get_nnz_sparse_Hess_Lagr() const { return nnz_sparse_Hess_Lagr; }
};

} // End ns hydrodynamics

} // End ns mfem

#endif // LAGRANGE_MULTIPLIER