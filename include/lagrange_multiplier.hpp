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
* Function: SetHiopConstraintGradSparsityPattern
*
* Purpose:
*  Set the sparsity pattern for the Jacobian matrix of the boundary constraints vector to be 
*  used in the Hiop mesh velocity solve.  Recall the the Gradient matrix is a SparseMatrix
*  of size |n_bdr_dofs| x 2*eta_geo
****************************************************************************************************/
inline void SetHiopBoundaryConstraintGradSparsityPattern(const Array<int> ess_tdofs, Array<int> &I, Array<int> &J, Array<double> &Data)
{
   // std::cout << "SetHiopBoundaryConstraintGradSparsityPattern\n";
   const int n_bdr_dofs = ess_tdofs.Size();

   I.SetSize(n_bdr_dofs + 1);
   J.SetSize(n_bdr_dofs);
   Data.SetSize(n_bdr_dofs);
   Data = 1;
   I[0] = 0;

   for (int bdr_it = 0; bdr_it < n_bdr_dofs; bdr_it++)
   {
      /* Set the (i+1)th entry of I and ith entry of J */
      I[bdr_it+1] = bdr_it+1;
      J[bdr_it] = ess_tdofs[bdr_it];
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
   I[0] = 0;
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

   /* Validate property of csr matrices */
   assert(I[2*num_vertices] == J.Size());
}


/****************************************************************************************************
* Function: SetHiopHessianSparsityPatternViscous
*
* Purpose:
*  Set the sparsity pattern for the Hessian matrix to be used in the 
*  Hiop mesh velocity solve with a Viscous function to minimize.  
*  Recall the the Hessian matrix is a SparseMatrix of size 
*  2*eta_geo x 2*eta_geo
*
*  For the viscous case, the sparsity pattern includes all neighbors
*  and neighbors of neighbors.
*
*  Note: According to hiop/src/Linalg/README.md, this matrix is upper triangular.
*
*        1x  2x  3x ... 1y 2y ...
*    1x  *
*    2x      *            *
*    ..           *  ...
*    1y       0          * 
*    2y                     *
*    ..
*
*  Upper triangular, CSR.
****************************************************************************************************/
template<int dim>
inline void SetHiopHessianSparsityPatternViscous(const ParMesh *pmesh, const Geometric<dim> &geom, const ParFiniteElementSpace &H1, const int &num_vertices, Array<int> &I, Array<int> &J)
{
   // std::cout << "SetHiopHessianSparsityPatternViscous\n";
   Array<int> adj_verts, adj_adj_verts, cols_for_node, cols_upper_diag;
   int size_adj_verts, size_adj_adj_verts;
   I.SetSize(num_vertices+1);
   I[0] = 0;
   J.SetSize(0);
   
   /* Iterate over nodes and build x portion of sparse matrix objects */
   for (int ix = 0; ix < num_vertices; ix++)
   {
      cols_for_node.SetSize(0);
      geom.VertexGetAdjacentVertices(ix, adj_verts);
      cols_for_node.Append(adj_verts);
      size_adj_verts = adj_verts.Size();

      /* Iterate over adjacent faces and append adjacent vertices with higher x index, always append y index as it is bigger by design */
      for (int j_it = 0; j_it < size_adj_verts; j_it++)
      {
         /* Get adj node index*/
         int jx = adj_verts[j_it];
         geom.VertexGetAdjacentVertices(jx, adj_adj_verts);
         cols_for_node.Append(adj_adj_verts);

         /* Fill adjacency array append y of adj_verts from constraints */
         cols_for_node.Append(jx + num_vertices); // always append y
      }

      /* Append diagonal element and sort */
      cols_for_node.Append(ix);
      cols_for_node.Sort();
      cols_for_node.Unique();
      int col_node_index = cols_for_node.Find(ix);
      cols_for_node.GetSubArray(col_node_index, cols_for_node.Size() - col_node_index, cols_upper_diag);

      // cout << "node " << ix;// << " has nonzero col: ";
      // cols_for_node.Print(cout);
      // cout << " has upper diag: ";
      // cols_upper_diag.Print(cout);

      /* Fill I */
      I[ix+1] = I[ix] + cols_upper_diag.Size(); 

      /* Append onto J */
      J.Append(cols_upper_diag);
   }

   // MFEM_ABORT("NOT FULLY IMPLEMENTED");

   /* 
   * Fill rest of I and J for y components. Upper Triangular 
   * Notice that the nonzero entries for the y components of the
   * velocity will be lower triangular and thus are omitted.
   * */
   const int base = I[num_vertices];

   for (int ix = 0; ix < num_vertices; ix++)
   {      
      /// TODO: Remove duplicated code.  Is there a way to put this code into the stuff above?
      geom.VertexGetAdjacentVertices(ix, adj_verts);
      size_adj_verts = adj_verts.Size();

      cols_for_node.SetSize(0);
      cols_for_node.Append(num_vertices + ix);

      /* iterate over faces and append onto cols_for_node those indices who are greated than i */
      for (int j_it = 0; j_it < size_adj_verts; j_it++)
      {
         int jx = adj_verts[j_it];
         geom.VertexGetAdjacentVertices(jx, adj_adj_verts);
         size_adj_adj_verts = adj_adj_verts.Size();
         for (int h_it = 0; h_it < size_adj_adj_verts; h_it++)
         {
            int hx = adj_adj_verts[h_it];
            if (hx > ix) { cols_for_node.Append(hx + num_vertices); }
         }

         /* Append if adj index is greater than one we are iterating on, upper triangular */
         if (jx > ix) { cols_for_node.Append(jx + num_vertices); }
      }

      /* Sort cols and append to J */
      cols_for_node.Sort();
      cols_for_node.Unique();
      assert(cols_for_node.Size() <= 13);
      J.Append(cols_for_node);
      // cout << "node " << ix << " has nonzero col: ";
      // cols_for_node.Print(cout);

      /* Fill i */
      I.Append(I[num_vertices + ix] + cols_for_node.Size());
   }

   /* Validate property of csr matrices */
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
   mutable SparseMatrix grad;

public:
   LocalMassConservationOperator(const Geometric<dim> &_geom, const ParGridFunction &_X, 
                                 const int &_num_cells, const int &_input_size, const double &_dt,
                                 const Array<int> GradCI, const Array<int> GradCJ)
      : Operator(_num_cells, _input_size),
        num_cells(_num_cells),
        input_size(_input_size),
        X(_X),  
        geom(_geom),
        dt(_dt),
        GradCIArr(GradCI),
        GradCJArr(GradCJ),
        GradDataArr(GradCJArr.Size()),
        grad(GradCIArr.GetData(), GradCJArr.GetData(), GradDataArr.GetData(), num_cells, input_size,
             false/*ownij*/, false/*owna*/, true/*is_sorted*/)
   {
      // cout << "LocalMassConservationOperator::Non-default constructor\n";
   }

   ~LocalMassConservationOperator()
   {
      // std::cout << "LocalMassConservationOperator - destructor\n";
   }

   void Mult(const Vector &v, Vector &y) const
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

         grad.SetSubMatrix(rows_arr, cols_arr, dm);
      }
      grad.Finalize();
   }

   /* Evaluate the gradient operator at the point v. */
   Operator &GetGradient(const Vector &v) const
   {
      // cout << "LMC::GetGradient\n";
      ComputeGradient(v);

      const int *In = grad.GetI(), *Jn = grad.GetJ();
      
      for (int i = 0, k=0; i < num_cells; i++)
      {
         for (int end = In[i+1]; k < end; k++)
         {
            int j = Jn[k];
         }
      }

      return grad;
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
         arrData[i] = 1.;
      }

      grad = new SparseMatrix(arrI.GetData(), arrJ.GetData(), arrData.GetData(), height, width, false, false, true);  

      // grad->Finalize();
   }
   ~zeroSparseMatrix()
   {
      delete grad;
      grad = nullptr;
   }

   void Mult(const Vector &v, Vector &y) const
   {
      y = 0.;
   }

   Operator &GetGradient(const Vector &v) const 
   {
      return *grad;
   } 
};


template<int dim>
class BoundaryConditionsOperator : public Operator
{
private:
   const int n_bdr_dofs, input_size;
   const int nnzSparse = 1;
   Array<int> ess_tdofs;
   Array<int> GradDIArr, GradDJArr;
   Array<double> GradDataArr;
   mutable SparseMatrix grad;

public:
   BoundaryConditionsOperator(const int &_n_bdr_dofs, const int &_input_size, const Array<int> _ess_tdofs,
                              const Array<int> GradDI, const Array<int> GradDJ, const Array<double> GradDData) :
      n_bdr_dofs(_n_bdr_dofs),
      input_size(_input_size),
      Operator(_n_bdr_dofs, _input_size),
      ess_tdofs(_ess_tdofs),
      GradDIArr(GradDI),
      GradDJArr(GradDJ),
      GradDataArr(GradDData),
      grad(GradDIArr.GetData(), GradDJArr.GetData(), GradDataArr.GetData(), n_bdr_dofs, input_size, false, false, true)
   {
      // std::cout << "BCoper constructor\n";
      assert(ess_tdofs.Size() == n_bdr_dofs);
   }

   void Mult(const Vector &v, Vector &y) const
   {
      y.SetSize(n_bdr_dofs);
      for (int i = 0; i < n_bdr_dofs; i++)
      {
         y[i] = v[ess_tdofs[i]];
      }
   }

   Operator &GetGradient(const Vector &v) const
   {
      return grad;
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
   const int n_bdr_dofs, num_cells, num_vertices;
   const int nnz_sparse_jaceq, nnz_sparse_jacineq, nnz_sparse_Hess_Lagr;
   const Vector &V_target;
   const Vector xmin, xmax;
   Vector massvec, bdr_vals, d_lo, d_hi;
   const LocalMassConservationOperator<dim> LMCoper;
   // const BoundaryConditionsOperator<dim> BCoper;
   Array<int> HessIArr, HessJArr;
   Array<double> HessData;
   Array<int> ess_tdofs;
   Array<int> BdrVertexIndexingArray;
   mutable SparseMatrix hess, hessf;
   DenseMatrix block;
   const ParGridFunction &X;
   const double dt;
   // const double interior_multiplier = 1.E-2;
   const double interior_multiplier = 1.;
   const bool add_obj_visc;

public:
   TargetOptimizedMeshVelocityProblem(
      const Geometric<dim> &_geom, const Vector &_V_target, const Vector &_massvec, 
      const ParGridFunction &_X, const int _num_cells,
      const double &dt, const Vector &_xmin, const Vector &_xmax,
      Array<int> _HessI, Array<int> _HessJ,
      Array<int> _GradCI, Array<int> _GradCJ,
      Array<int> _GradDI, Array<int> _GradDJ,
      Array<double> _GradDData, Array<double> &_bdr_vals, 
      Array<int> _ess_tdofs, Array<int> _BdrVertexIndexingArray,
      const bool &_add_obj_visc) 
      : geom(_geom),
        X(_X),
        dt(dt),
        num_cells(_num_cells),
        num_vertices(_geom.GetNVDofs_H1()),
        xmin(_xmin), xmax(_xmax),
        OptimizationProblem(dim*_geom.GetNVDofs_H1(), NULL, NULL),
        V_target(_V_target), 
        massvec(_massvec), 
        bdr_vals(_bdr_vals, _bdr_vals.Size()),
      //   d_lo(_massvec), d_hi(_massvec),
        n_bdr_dofs(_ess_tdofs.Size()),
      //   d_lo(bdr_vals), d_hi(bdr_vals),
      //   BCoper(n_bdr_dofs, input_size, _ess_tdofs, _GradDI, _GradDJ, _GradDData),
        LMCoper(_geom, _X, _num_cells, input_size, dt, _GradCI, _GradCJ),
        HessIArr(_HessI), HessJArr(_HessJ), HessData(_HessJ.Size()),
        ess_tdofs(_ess_tdofs),
        BdrVertexIndexingArray(_BdrVertexIndexingArray),
        nnz_sparse_jaceq(_GradCJ.Size()),
      //   nnz_sparse_jacineq(n_bdr_dofs),
        nnz_sparse_jacineq(0),
        nnz_sparse_Hess_Lagr(_HessJ.Size()),
        hess(HessIArr.GetData(), HessJArr.GetData(), HessData.GetData(), input_size, 
             input_size, false/*ownij*/, false/*owna*/, true/*is_sorted*/),
        hessf(HessIArr.GetData(), HessJArr.GetData(), HessData.GetData(), input_size, 
             input_size, false/*ownij*/, false/*owna*/, true/*is_sorted*/),
        block(4),
        add_obj_visc(_add_obj_visc)
   {
      // cout << "TOMVProblem constructor\n";

      // Problem assumes that vectors are of the correct size
      assert(_massvec.Size() == _num_cells);
      assert(_V_target.Size() == input_size);
      assert(_bdr_vals.Size() == n_bdr_dofs);

      C = &LMCoper;
      SetEqualityConstraint(massvec);

      // D = &LMCoper;
      // double tol = 1.;
      // d_lo -= tol;
      // d_hi += tol;
      // SetInequalityConstraint(d_lo, d_hi);
      // D = &BCoper;
      // d_lo -= 1.E-6, d_hi += 1.E-6;
      // SetInequalityConstraint(d_lo, d_hi);

      SetSolutionBounds(xmin, xmax);

      /* Construct block object used in building the Hessian matrix */
      block(0,1) = -0.5, block(0,3) = 0.5;
      block(1,0) = 0.5, block(1,2) = -0.5;
      block(2,1) = 0.5, block(2,3) = -0.5;
      block(3,0) = -0.5, block(3,2) = 0.5;

      block *= pow(dt, 2);

      /* Compute hess of objective function since this does not change in time */
      ComputeObjectiveHessFData();
   }

   ~TargetOptimizedMeshVelocityProblem()
   {
      // std::cout << "TargetOptimizedMeshVelocityProblem - destructor\n";
   }

   double CalcObjective(const Vector &V) const
   {
      // cout << "TargetOptimizedMeshVelocityProblem::CalcObjective\n";
      double val = 0., visc = 0.;
      Vector Vi(dim), Vi_hat(dim), Vi_target(dim), Vj(dim), temp_vec(dim);
      Array<int> adj_verts;

      /* Iterate over geometric vertices */
      for (int ix = 0; ix < num_vertices; ix++)
      {
         /* Grab Vi */
         geom.GetNodeVelocityVecL(V, ix, Vi);
         geom.GetNodeVelocityVecL(V_target, ix, Vi_target);
         subtract(Vi, Vi_target, temp_vec);
   
         /* contribution from x coord */
         double ix_val = temp_vec[ix] * temp_vec[ix];
         if (ess_tdofs.Find(ix) != -1) {
            val += ix_val; 
         } else { // interior node
            val += interior_multiplier * ix_val;
         }
         /* contribution from y coord */
         int iy = ix + num_vertices;
         double iy_val = temp_vec[iy] * temp_vec[iy];
         if (ess_tdofs.Find(iy) != -1) {
            val += iy_val; 
         } else { // interior node
            val += interior_multiplier * iy_val;
         }

         /* Optionally add in viscosity at interior geometric nodes */
         if (add_obj_visc && BdrVertexIndexingArray[ix] == -1)
         {
            /* Grab Vi */
            geom.GetNodeVelocityVecL(V, ix, Vi);

            /* Compute Vi_hat */
            Vi_hat = 0.;
            geom.VertexGetAdjacentVertices(ix, adj_verts);
            int adj_verts_size = adj_verts.Size();
            assert(adj_verts_size == 4);

            for (int j_it = 0; j_it < adj_verts_size; j_it++)
            {
               int jx = adj_verts[j_it];
               geom.GetNodeVelocityVecL(V, jx, Vj);
               Vi_hat += Vj;
            }

            /* Compute norm and add to val */
            Vi_hat *= .25;
            subtract(Vi, Vi_hat, temp_vec);
            visc += std::pow(temp_vec.Norml2(),2);
         } // End interior node viscosity
      } // End geometric vertices loops

      // std::cout //<< "TargetOptimizedMeshVelocityProblem::CalcObjective\n"
      //           << "Target val: " << std::setw(8) << val << std::endl
      //           << "Interior visc: " << std::setw(8) << visc << std::endl;

      val += visc;
      return val;
   }

   int get_input_size() 
   {
      return input_size;
   }

   void CalcObjectiveGrad(const Vector &V, Vector &grad) const 
   {
      // cout << "TargetOptimizedMeshVelocityProblem::CalcObjectiveGrad\n";
      assert(grad.Size() == input_size);
      assert(V.Size() == input_size);

      Vector Vi(dim), Vi_hat(dim), Vj(dim), temp_vec(dim);
      Array<int> adj_verts;
      int jx, iy, jy;

      /* Reset gradient */
      grad = 0.;

      /* Iterate over geometric vertices */
      for (int ix = 0; ix < num_vertices; ix++)
      {
         /* contribution from x coord */
         double ix_val = 2.*(V(ix) - V_target(ix));
         if (ess_tdofs.Find(ix) != -1) {
            grad(ix) += ix_val;
         } else {
            grad(ix) += interior_multiplier * ix_val;
         }
         /* contribution from x coord */
         iy = ix + num_vertices;
         double iy_val = 2.*(V(iy) - V_target(iy));
         if (ess_tdofs.Find(iy) != -1) {
            grad(iy) += iy_val;
         } else {
            grad(iy) += interior_multiplier * iy_val;
         }

         /* For viscosity, we only care about interior vertices */
         if (add_obj_visc && BdrVertexIndexingArray[ix] == -1)
         {
            geom.GetNodeVelocityVecL(V, ix, Vi);
            /* Compute Vi_hat */
            Vi_hat = 0.;
            geom.VertexGetAdjacentVertices(ix, adj_verts);
            int adj_verts_size = adj_verts.Size();
            assert(adj_verts_size == 4);

            for (int j_it = 0; j_it < adj_verts_size; j_it++)
            {
               jx = adj_verts[j_it];
               geom.GetNodeVelocityVecL(V, jx, Vj);
               Vi_hat += Vj;
            }

            /* Add interior portion to d/dvix, d/dviy */
            add(2, Vi, -.5, Vi_hat, temp_vec); // 32Vi - 8Vi_hat
            grad[ix] += temp_vec[0];
            grad[iy] += temp_vec[1];

            /* iterate over adjacent vertices */
            for (int j_it = 0; j_it < adj_verts_size; j_it++)
            {  
               /* add adjacency portion to d/dvjx, d/dvjy */
               jx = adj_verts[j_it];
               jy = jx + num_vertices;
               add(-.5, Vi, .125, Vi_hat, temp_vec);
               grad[jx] += temp_vec[0];
               grad[jy] += temp_vec[1];
            } // End interior vertices if statement
         } // End add_obj_visc
      } // End geometric vertices loops
   }

   void ComputeObjectiveHessFData() const
   {
      Array<int> adj_verts;
      int size_adj_verts;

      /* Reset data in sparse Hessian since we build it through addition */
      double * _dataf = hessf.GetData();
      for (int d = 0; d < HessJArr.Size(); d++)
      {
         _dataf[d] = 0.;
      }

      /* This target function imposes nonzero quantities dependent on the adjacency */
      for (int ix = 0; ix < num_vertices; ix++)
      {
         /* Target velocity portion */
         /* contribution from x coord */
         if (ess_tdofs.Find(ix) != -1) {
            hessf.Elem(ix,ix) += 2.; 
         } else { // interior node
            hessf.Elem(ix,ix) += 2. * interior_multiplier;
         }
         /* contribution from y coord */
         int iy = ix + num_vertices;
         if (ess_tdofs.Find(iy) != -1) {
            hessf.Elem(iy,iy) += 2.;
         } else { // interior node
            hessf.Elem(iy,iy) += 2. * interior_multiplier;
         }

         /* For the viscosity component, we only care about the interior vertices */
         if (add_obj_visc && BdrVertexIndexingArray[ix] == -1)
         {
            /* Add to diagonal elements for both x and y coords */
            hessf.Elem(ix,ix) += 2.;
            hessf.Elem(iy,iy) += 2.;

            /* Get adjacent vertices */
            geom.VertexGetAdjacentVertices(ix, adj_verts);
            int adj_verts_size = adj_verts.Size();
            assert(adj_verts_size == 4);

            /* Iterate over adjacent vertices */
            for (int j_it = 0; j_it < adj_verts_size; j_it++)
            {
               int jx = adj_verts[j_it];
               assert(jx < num_vertices);
               int jy = jx + num_vertices;
               assert(jy < input_size);

               /* Contribution to off diagonal that depends on i */
               if (jx > ix)
               {
                  hessf.Elem(ix,jx) -= .5;
                  hessf.Elem(iy,jy) -= .5;
               } 
               else
               {
                  hessf.Elem(jx,ix) -= .5;
                  hessf.Elem(jy,iy) -= .5;
               }

               /* Contribution to entries not dependent on i (from adjacency) */
               for (int h_it = 0; h_it < adj_verts_size; h_it++)
               {
                  int hx = adj_verts[h_it];
                  assert(hx < num_vertices);
                  int hy = hx + num_vertices;
                  assert(hy < input_size);

                  /* Only add in this contribution once per pair */
                  if (hx >= jx)
                  {
                     hessf.Elem(jx,hx) += .125;
                     hessf.Elem(jy,hy) += .125;
                  }
               } // end nested adj verts loop
            } // End interior vertices/adj_verts loop if statement
         } // End add_obj_visc
      } // End geometric vertices loops
   }

   void ComputeObjectiveHessData(const Vector &x) const
   {
      // std::cout << "OptimizedMeshVelocityProblem::ComputeObjectiveHessData\n";
      Array<int> verts, Vy_arr;
      DenseMatrix dm;

      /* Reset data in sparse Hessian since we build it through addition */
      double * _data = hess.GetData();
      double * _dataf = hessf.GetData();
      for (int d = 0; d < HessJArr.Size(); d++)
      {
         _data[d] = _dataf[d];
      }

      /* Verify our Hessian is starting out zeroed out */
      if (hess.MaxNorm() > hessf.MaxNorm())
      {
         std::cout << "max norm: " << hess.MaxNorm() << std::endl;
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

            hess.AddSubMatrix(verts, Vy_arr, dm, 2/*skip_zeros*/);
         }
      } // End cell it
   }

   virtual Operator & CalcObjectiveHess(const Vector &x) const 
   {
      // std::cout << "OptimizedMeshVelocityProblem::CalcObjectiveHessian\n";
      ComputeObjectiveHessData(x);

      // std::cout << "OptimizedMeshVelocityProblem::CalcObjectiveHessian - DONE\n";
      return hess;
   }

   virtual int get_nnz_sparse_Jaceq() const { return nnz_sparse_jaceq; }
   virtual int get_nnz_sparse_Jacineq() const { return nnz_sparse_jacineq; }
   // virtual int get_nnz_sparse_Jacineq() const { return 0; }
   virtual int get_nnz_sparse_Hess_Lagr() const { return nnz_sparse_Hess_Lagr; }
};


/**
 * Conservative a posteriori correction to mesh velocity to guarantee local mass conservation:
 * Find V that minimizes SUM_F || VF1 - VF2 ||^2, subject to
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
   const int num_cells, num_faces, num_vertices;
   const int nnz_sparse_jaceq, nnz_sparse_Hess_Lagr;
   Vector massvec, d_lo, d_hi;
   const LocalMassConservationOperator<dim> LMCoper;
   // const zeroSparseMatrix<dim> zSMoper;
   Array<int> HessIArr, HessJArr;
   Array<double> HessData;
   Array<int> ess_tdofs;
   Array<int> BdrVertexIndexingArray;
   mutable SparseMatrix hessf, hess;
   DenseMatrix block;

public:
   ViscousOptimizedMeshVelocityProblem(const Geometric<dim> &_geom, const Vector &_massvec,
                                       const ParGridFunction &_X, const int _num_cells,
                                       const double &dt, const Vector &xmin, const Vector &xmax,
                                       const Array<int> &_HessI, const Array<int> &_HessJ,
                                       const Array<int> &GradCI, const Array<int> &GradCJ,
                                       const Array<int> _ess_tdofs, const Array<int> _BdrVertexIndexingArray)
      : geom(_geom),
        num_cells(_num_cells),
        num_faces(geom.GetNumFaces()),
        OptimizationProblem(dim*_geom.GetNVDofs_H1(), NULL, NULL),
        massvec(_massvec), d_lo(1), d_hi(1),
        LMCoper(_geom, _X, _num_cells, input_size, dt, GradCI, GradCJ),
        nnz_sparse_jaceq(GradCJ.Size()),
        nnz_sparse_Hess_Lagr(_HessJ.Size()),
        num_vertices(_geom.GetNVDofs_H1()),
        HessIArr(_HessI), HessJArr(_HessJ), HessData(_HessJ.Size()),
        ess_tdofs(_ess_tdofs),
        BdrVertexIndexingArray(_BdrVertexIndexingArray),
        hessf(HessIArr.GetData(), HessJArr.GetData(), HessData.GetData(), input_size, 
             input_size, false/*ownij*/, false/*owna*/, true/*is_sorted*/),
        hess(HessIArr.GetData(), HessJArr.GetData(), HessData.GetData(), input_size, 
             input_size, false/*ownij*/, false/*owna*/, true/*is_sorted*/),
        block(4)
   {
      assert(_massvec.Size() == _num_cells);
      // std::cout << "mass vec size: " << _massvec.Size() << std::endl;
      // std::cout << "xmin size: " << xmin.Size() << ", xmax size: " << xmax.Size() << std::endl;
      // std::cout << "input size: " << input_size << endl;

      C = &LMCoper;
      SetEqualityConstraint(massvec);

      SetSolutionBounds(xmin, xmax);

      /* Construct block object used in building the Hessian matrix */
      block(0,1) = -0.5, block(0,3) = 0.5;
      block(1,0) = 0.5, block(1,2) = -0.5;
      block(2,1) = 0.5, block(2,3) = -0.5;
      block(3,0) = -0.5, block(3,2) = 0.5;

      block *= pow(dt, 2);

      /* Compute hessf since this does not change in time */
      ComputeObjectiveHessFData();
   }

   ~ViscousOptimizedMeshVelocityProblem()
   {
   }

   double CalcObjective(const Vector &V) const
   {
      double val = 0.;
      Vector Vi(dim), Vi_hat(dim), Vj(dim), temp_vec(dim);
      Array<int> adj_verts;

      /* Iterate over geometric vertices */
      for (int ix = 0; ix < num_vertices; ix++)
      {
         /* We only care about interior vertices */
         if (BdrVertexIndexingArray[ix] == -1)
         {
            /* Grab Vi */
            geom.GetNodeVelocityVecL(V, ix, Vi);

            /* Compute Vi_hat */
            Vi_hat = 0.;
            geom.VertexGetAdjacentVertices(ix, adj_verts);
            int adj_verts_size = adj_verts.Size();
            assert(adj_verts_size == 4);

            for (int j_it = 0; j_it < adj_verts_size; j_it++)
            {
               int jx = adj_verts[j_it];
               geom.GetNodeVelocityVecL(V, jx, Vj);
               Vi_hat += Vj;
            }

            /* Compute norm and add to val */
            Vi_hat *= .25;
            subtract(Vi, Vi_hat, temp_vec);
            val += std::pow(temp_vec.Norml2(),2);
         } // End interior vertices if statement
      } // End geometric vertices loops

      return val;
   }

   void CalcObjectiveGrad(const Vector &V, Vector &grad) const
   {
      assert(grad.Size() == input_size);
      assert(V.Size() == input_size);

      Vector Vi(dim), Vi_hat(dim), Vj(dim), temp_vec(dim);
      Array<int> adj_verts;
      int jx, iy, jy;

      /* Reset gradient */
      grad = 0.;

      /* Iterate over geometric vertices */
      for (int ix = 0; ix < num_vertices; ix++)
      {
         /* We only care about interior vertices */
         if (BdrVertexIndexingArray[ix] == -1)
         {
            geom.GetNodeVelocityVecL(V, ix, Vi);
            /* Compute Vi_hat */
            Vi_hat = 0.;
            geom.VertexGetAdjacentVertices(ix, adj_verts);
            int adj_verts_size = adj_verts.Size();
            assert(adj_verts_size == 4);

            for (int j_it = 0; j_it < adj_verts_size; j_it++)
            {
               jx = adj_verts[j_it];
               geom.GetNodeVelocityVecL(V, jx, Vj);
               Vi_hat += Vj;
            }

            /* Add interior portion to d/dvix, d/dviy */
            iy = ix + num_vertices;
            add(2, Vi, -.5, Vi_hat, temp_vec); // 32Vi - 8Vi_hat
            grad[ix] += temp_vec[0];
            grad[iy] += temp_vec[1];

            /* iterate over adjacent vertices */
            for (int j_it = 0; j_it < adj_verts_size; j_it++)
            {  
               /* add adjacency portion to d/dvjx, d/dvjy */
               jx = adj_verts[j_it];
               jy = jx + num_vertices;
               add(-.5, Vi, .125, Vi_hat, temp_vec);
               grad[jx] += temp_vec[0];
               grad[jy] += temp_vec[1];
            }
         } // End interior vertices if statement
      } // End geometric vertices loops

      // cout << "grad: ";
      // grad.Print(cout);
      // MFEM_ABORT("NOT FULLY IMPLEMENTED\n");
   }

   void ComputeObjectiveHessFData() const
   {
      Array<int> adj_verts;
      int size_adj_verts;

      /* Reset data in sparse Hessian since we build it through addition */
      double * _dataf = hessf.GetData();
      for (int d = 0; d < HessJArr.Size(); d++)
      {
         _dataf[d] = 0.;
      }

      /* This target function imposes nonzero quantities dependent on the adjacency */
      for (int ix = 0; ix < num_vertices; ix++)
      {
         /* We only care about interior vertices */
         if (BdrVertexIndexingArray[ix] == -1)
         {
            /* Add to diagonal elements for both x and y coords */
            int iy = ix + num_vertices;
            hessf.Elem(ix,ix) += 2.;
            hessf.Elem(iy,iy) += 2.;

            /* Get adjacent vertices */
            geom.VertexGetAdjacentVertices(ix, adj_verts);
            int adj_verts_size = adj_verts.Size();
            assert(adj_verts_size == 4);

            /* Iterate over adjacent vertices */
            for (int j_it = 0; j_it < adj_verts_size; j_it++)
            {
               int jx = adj_verts[j_it];
               assert(jx < num_vertices);
               int jy = jx + num_vertices;
               assert(jy < input_size);

               /* Contribution to off diagonal that depends on i */
               if (jx > ix)
               {
                  hessf.Elem(ix,jx) -= .5;
                  hessf.Elem(iy,jy) -= .5;
               } 
               else
               {
                  hessf.Elem(jx,ix) -= .5;
                  hessf.Elem(jy,iy) -= .5;
               }

               /* Contribution to entries not dependent on i (from adjacency) */
               for (int h_it = 0; h_it < adj_verts_size; h_it++)
               {
                  int hx = adj_verts[h_it];
                  assert(hx < num_vertices);
                  int hy = hx + num_vertices;
                  assert(hy < input_size);

                  /* Only add in this contribution once per pair */
                  if (hx >= jx)
                  {
                     hessf.Elem(jx,hx) += .125;
                     hessf.Elem(jy,hy) += .125;
                  }
               } // end nested adj verts loop
            } // end adj verts loop
         } // End interior vertices if statement
      } // End geometric vertices loops
   }

   void ComputeObjectiveHessData(const Vector &x) const
   {
      // std::cout << "ViscousOptimizedMeshVelocityProblem::ComputeObjectiveHessData\n";
      Array<int> adj_verts, adj_verts_y;
      DenseMatrix dm;
      int size_adj_verts;

      /* Reset data in sparse Hessian since we build it through addition */
      double * _data = hess.GetData();
      double * _dataf = hessf.GetData();
      for (int d = 0; d < HessJArr.Size(); d++)
      {
         _data[d] = _dataf[d];
      }

      /* Verify our Hessian is starting out zeroed out */
      if (hess.MaxNorm() > hessf.MaxNorm())
      {
         std::cout << "max norm: " << hess.MaxNorm() << std::endl;
         MFEM_ABORT("Hessian has not been reset properly.\n")
      }

      /* Fill data array cell by cell from local mass constraints */
      for (int cell_it = 0; cell_it < num_cells; cell_it++)
      {
         /* Only add the contribution to the hessian if the LM is nonzero */
         double coeff = 0.5 * this->lambda[cell_it];
         if (coeff != 0.)
         {
            /* Get adjacent vertices */
            geom.GetParMesh()->GetElementVertices(cell_it, adj_verts);
            size_adj_verts = adj_verts.Size();

            adj_verts_y.SetSize(size_adj_verts);
            for (int i = 0; i < size_adj_verts; i++)
            {
               adj_verts_y[i] = adj_verts[i] + num_vertices;
            }

            /* top right block */
            dm = block;
            dm *= this->lambda[cell_it];

            hess.AddSubMatrix(adj_verts, adj_verts_y, dm, 2/*skip_zeros*/);
         }
      } // End cell it
   }

   Operator & CalcObjectiveHess(const Vector &x) const
   {
      ComputeObjectiveHessData(x);

      return hess;
   }
   int get_input_size() const { return input_size; }
   int get_nnz_sparse_Jaceq() const { return nnz_sparse_jaceq; }
   int get_nnz_sparse_Jacineq() const { return 0; }
   int get_nnz_sparse_Hess_Lagr() const { return nnz_sparse_Hess_Lagr; }
};

} // End ns hydrodynamics

} // End ns mfem

#endif // LAGRANGE_MULTIPLIER