#include "geometry.hpp"


namespace mfem 
{

namespace hydroLO
{

/**
* Constructor 
*/
template<int dim>
Geometric<dim>::Geometric(Array<int> offsets, ParFiniteElementSpace &h1, ParFiniteElementSpace &l2) :
   H1(h1),
   L2(l2),
   pmesh(h1.GetParMesh()),
   block_offsets(offsets),
   num_faces(L2.GetNF()),
   edge_vertex(pmesh->GetEdgeVertexTable())
{
   cout << "Geometric::Non-default constructor\n";

   NDofs_H1 = h1.GetNDofs();
   NVDofs_H1 = h1.GetNVDofs();
   Transpose(*edge_vertex, vertex_edge);
   cout << "NDofs_H1: " << NDofs_H1 << endl;
}

/****************************************************************************************************
* Function: UpdateNodeVelocity
* Parameters:
*   dSdt - BlockVector representing time derivative of FiniteElement information
*   node - Global index of node in question.
*   vel  - Velocity computed at node to be used to move the mesh
* Purpose:
*   This function is used to update the mv_gf ParGridFunction which is used to move the mesh.
****************************************************************************************************/
template<int dim>
void Geometric<dim>::UpdateNodeVelocity(Vector &dSdt, const int & node, const Vector & vel) const
{
   ParGridFunction dxdt;
   MFEM_ABORT("Need to remove mv_gf implementation\n");
   dxdt.MakeRef(&H1, dSdt, block_offsets[0]);
   assert(dxdt.Size() == dim*NDofs_H1);

   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      dxdt[index] = vel[i];
   }
}


/****************************************************************************************************
* Function: UpdateNodeVelocity
* Parameters:
*   v_gf - ParGridFunction for the velocities
*   node - Global index of node in question.
*   vel  - Velocity computed at node to be used to move the mesh
* Purpose:
*   This function is used to update the mv_gf ParGridFunction which is used to move the mesh.
****************************************************************************************************/
template<int dim>
void Geometric<dim>::UpdateNodeVelocity(ParGridFunction &dxdt, const int &node, const Vector &vel) const
{
   assert(dxdt.Size() == dim*NDofs_H1);
   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      dxdt[index] = vel[i];
   }
}


/****************************************************************************************************
* Function: UpdateNodeVelocityVecL
* Parameters:
*   dxdt_l  - Vector on the first order FESpace representing velocity on geometric
*             corner nodes of the mesh
*   node    - Global index of node in question.
*   vel     - Velocity at given node
* Purpose:
*   This function is used to update the dxdt ParGridFunction which is used to move the mesh.
****************************************************************************************************/
template<int dim>
void Geometric<dim>::UpdateNodeVelocityVecL(Vector &dxdt_l, const int & node, const Vector &vel) const
{
   assert(dxdt_l.Size() == dim*NVDofs_H1);
   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NVDofs_H1;
      dxdt_l[index] = vel[i];
   }

}


/****************************************************************************************************
* Function: GetNodeVelocity
* Parameters:
*   dSdt - BlockVector representing time derivative of FiniteElement information
*   node - Global index of node in question.
*   vel  - Velocity at given node
*
* Purpose:
*  This function returns the velocity at the given global node.
****************************************************************************************************/
template<int dim>
void Geometric<dim>::GetNodeVelocity(const Vector &dSdt, const int & node, Vector & vel) const
{
   Vector* sptr = const_cast<Vector*>(&dSdt);
   ParGridFunction dxdt;
   dxdt.MakeRef(&H1, *sptr, block_offsets[0]);
   assert(dxdt.Size() == dim*NDofs_H1);

   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      vel[i] = dxdt[index];
   }
}  


/****************************************************************************************************
* Function: GetNodeVelocityGF
* Parameters:
*   dxdt  - Vector of mesh nodal velocities 
*   node  - Global index of node in question.
*   vel   - Velocity at given node
*
* Purpose:
*  This function returns the velocity at the given global node.
****************************************************************************************************/
template<int dim>
void Geometric<dim>::GetNodeVelocity(const ParGridFunction &dxdt, const int & node, Vector & vel) const
{
   assert(dxdt.Size() == dim*NDofs_H1);
   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      vel[i] = dxdt[index];
   }
}

/****************************************************************************************************
* Function: GetNodeVelocityVecL
* Parameters:
*   dxdt_l  - Vector on the first order FESpace representing velocity on geometric
*             corner nodes of the mesh
*   node    - Global index of node in question.
*   vel     - Velocity at given node
*
* Purpose:
*  This function returns the velocity at the given global node.
* Note: L indicates grid function is the low order form.
****************************************************************************************************/
template<int dim>
void Geometric<dim>::GetNodeVelocityVecL(const Vector &dxdt_l, const int & node, Vector & vel) const
{
   assert(dxdt_l.Size() == dim*NVDofs_H1);
   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NVDofs_H1;
      vel[i] = dxdt_l[index];
   }
}


/****************************************************************************************************
* Function: UpdateNodePosition
* Parameters:
*   S    - BlockVector representing FiniteElement information
*   node - Global index of node in question.
*   x    - Coordinate position at given node.
*
* Purpose:
*  This function updates the cartesian location corresponding to a global node.
****************************************************************************************************/
template<int dim>
void Geometric<dim>::UpdateNodePosition(Vector &S, const int & node, const Vector &x) const
{
   ParGridFunction x_gf;
   x_gf.MakeRef(&H1, S, block_offsets[0]);
   assert(x_gf.Size() == dim*NDofs_H1);

   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      x_gf[index] = x[i];
   }
}


/****************************************************************************************************
* Function: GetNodePosition
* Parameters:
*   S    - BlockVector representing FiniteElement information
*   node - Global index of node in question.
*   x    - Coordinate position at given node.
*
* Purpose:
*  This function returns the cartesian location corresponding to a global node.
****************************************************************************************************/
template<int dim>
void Geometric<dim>::GetNodePositionFromBV(const Vector &S, const int & node, Vector & x) const
{
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction x_gf;
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);
   assert(x_gf.Size() == dim*NDofs_H1);

   GetNodePosition(x_gf, node, x);
}


/****************************************************************************************************
* Function: GetNodePosition
* Parameters:
*   x_gf - ParGridFunction representing nodal locations
*   node - Global index of node in question.
*   x    - Coordinate position at given node.
*
* Purpose:
*  This function returns the cartesian location corresponding to a global node.
****************************************************************************************************/
template<int dim>
void Geometric<dim>::GetNodePosition(const ParGridFunction &x_gf, const int & node, Vector &x) const
{
   assert(x_gf.Size() == dim*NDofs_H1);
   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      x[i] = x_gf[index];
   }
}


/****************************************************************************************************
* Function: GetFaceDofs
* Parameters:
*   face      - integer representing the face that dofs are requested from
*   face_dofs - Integer array of face dofs 
*
* Purpose:
*    This function returns the H1 DOFS on a given face.
****************************************************************************************************/
template<int dim>
void Geometric<dim>::GetFaceDofs(const int &face, Array<int> &face_dofs) const
{
   H1.GetFaceDofs(face, face_dofs);
}


/****************************************************************************************************
* Function: VertexGetAdjacentFaces
* Parameters:
*   vertex       - integer representing the vertex
*   adj_vertices - array of adjacent face indices
*
* Purpose:
*    This function returns the H1 adjacent faces
****************************************************************************************************/
template<int dim>
void Geometric<dim>::VertexGetAdjacentFaces(const int &vertex, Array<int> &adj_faces) const
{
   vertex_edge.GetRow(vertex, adj_faces);
}


/****************************************************************************************************
* Function: GetNumAdjFaces
* Parameters:
*   vertex       - integer representing the vertex
*
* Purpose:
*    This function returns the number of adjacent faces
****************************************************************************************************/
template<int dim>
int Geometric<dim>::GetNumAdjFaces(const int &vertex) const
{
   Array<int> adj_faces;
   vertex_edge.GetRow(vertex, adj_faces);
   return adj_faces.Size();
}


/****************************************************************************************************
* Function: VertexGetAdjacentVertices
* Parameters:
*   vertex       - integer representing the vertex
*   adj_vertices - array of adjacent vertex indices
*
* Purpose:
*    This function returns the H1L adjacent vertices
****************************************************************************************************/
template<int dim>
void Geometric<dim>::VertexGetAdjacentVertices(const int &vertex, Array<int> &adj_vertices) const
{
   Array<int> faces_row, face_dofs;
   vertex_edge.GetRow(vertex, faces_row);
   int num_faces = faces_row.Size();
   adj_vertices.SetSize(num_faces);

   for (int face_it = 0; face_it < num_faces; face_it++)
   {
      int face = faces_row[face_it];
      H1.GetFaceDofs(face, face_dofs);
      int face_dof1 = face_dofs[0], face_dof2 = face_dofs[1];
      if (vertex == face_dof1) {
         adj_vertices[face_it] = face_dof2;
      } else {
         adj_vertices[face_it] = face_dof1;
      }
   }
}


/****************************************************************************************************
* Function: Orthogonal
* Parameters:
*     v - Vector to be rotated
*
* Purpose:
*     Rotate a vector counter-clockwise of angle pi/2.
*
* Example:
*  (-1,1) ---> (-1,-1)
*  (x, y) ---> (-y, x)
****************************************************************************************************/
template<int dim>
void Geometric<dim>::Orthogonal(Vector &v) const
{
   assert(v.Size() == dim);
   if (dim == 2)
   {
      double x = v(0), y = v(1);
      v(0) = -1. * y;
      v(1) = x;
   }
   else
   {
      MFEM_ABORT("Can only rotate 2D vectors\n");
   }
}


/****************************************************************************************************
* Function: Perpendicular
* Parameters:
*     v - Vector to be rotated
*
* Purpose:
*     Rotate a vector clockwise of angle pi/2.
*
* Example:
*  (1,0) ---> (0,-1)
****************************************************************************************************/
template<int dim>
void Geometric<dim>::Perpendicular(Vector &v) const
{
   assert(v.Size() == dim);
   if (dim == 2)
   {
      double x = v(0), y = v(1);
      v(0) = y;
      v(1) = -1. * x;
   }
   else
   {
      MFEM_ABORT("Can only rotate 2D vectors\n");
   }
}

/* Explicit instantiation */
template class Geometric<1>;
template class Geometric<2>;
template class Geometric<3>;

} // ns hydroLO
} // ns mfem