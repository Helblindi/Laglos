#include "geometry.hpp"


namespace mfem 
{

namespace hydrodynamics
{

/**
* Constructor 
*/
template<int dim>
Geometric<dim>::Geometric(Array<int> offsets, ParFiniteElementSpace &h1) :
   H1(h1),
   pmesh(H1.GetParMesh()),
   block_offsets(offsets),
   NDofs_H1(H1.GetNDofs()) 
   {}

/****************************************************************************************************
* Function: UpdateNodeVelocity
* Parameters:
*   S    - BlockVector representing FiniteElement information
*   node - Global index of node in question.
*   vel  - Velocity computed at node to be used to move the mesh
* Purpose:
*   This function is used to update the mv_gf ParGridFunction which is used to move the mesh.
****************************************************************************************************/
template<int dim>
void Geometric<dim>::UpdateNodeVelocity(Vector &S, const int & node, const Vector & vel) const
{
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction mv_gf;
   mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);

   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      mv_gf[index] = vel[i];
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
void Geometric<dim>::UpdateNodeVelocity(ParGridFunction &v_gf, const int &node, const Vector &vel) const
{
   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      v_gf[index] = vel[i];
   }
}


/****************************************************************************************************
* Function: GetNodeVelocity
* Parameters:
*   S    - BlockVector representing FiniteElement information
*   node - Global index of node in question.
*   vel  - Velocity at given node
*
* Purpose:
*  This function returns the velocity at the given global node.
****************************************************************************************************/
template<int dim>
void Geometric<dim>::GetNodeVelocity(const Vector &S, const int & node, Vector & vel) const
{
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction mv_gf;
   mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);

   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      vel[i] = mv_gf[index];
   }
}  


/****************************************************************************************************
* Function: GetNodeVelocityGF
* Parameters:
*   S    - BlockVector representing FiniteElement information
*   node - Global index of node in question.
*   vel  - Velocity at given node
*
* Purpose:
*  This function returns the velocity at the given global node.
****************************************************************************************************/
template<int dim>
void Geometric<dim>::GetNodeVelocity(const ParGridFunction &mv_gf, const int & node, Vector & vel) const
{
   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      vel[i] = mv_gf[index];
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
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction x_gf;
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);

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
   // cout << "GetNodePosition\n";
   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NDofs_H1;
      x[i] = x_gf[index];

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

} // ns hydrodynamics
} // ns mfem