#ifndef GEOMETRY
#define GEOMETRY

#include "mfem.hpp"

using namespace std;

namespace mfem
{

namespace hydrodynamics
{

template <int dim>
class Geometric
{

private:
   Array<int> block_offsets;
   ParFiniteElementSpace &H1, &L2;
   ParMesh *pmesh;
   int NDofs_H1;
   int NVDofs_H1;
   const int num_faces;
   Table vertex_edge, *edge_vertex;

public:
   Geometric(Array<int> offsets, ParFiniteElementSpace &h1, ParFiniteElementSpace &l2);
   ~Geometric() {}

   inline Array<int> GetBlockOffsets() const { return block_offsets; }
   inline ParFiniteElementSpace& GetPFES() const { return H1; }
   inline ParMesh *GetParMesh() const { return pmesh; }
   inline int GetNDofs_H1() const { return NDofs_H1; }
   inline int GetNVDofs_H1() const { return NVDofs_H1; }
   inline int GetNumFaces() const { return num_faces; }

   void UpdateNodeVelocity(Vector &S, const int & node, const Vector & vel) const;
   void UpdateNodeVelocity(ParGridFunction &mv_gf, const int &node, const Vector &vel) const;
   void UpdateNodeVelocityVecL(Vector &mv_gf_l, const int & node, const Vector &x) const;
   void GetNodeVelocity(const Vector &S, const int & node, Vector & vel) const;
   void GetNodeVelocity(const ParGridFunction &mv_gf, const int & node, Vector & vel) const;
   void GetNodeVelocityVecL(const Vector &mv_gf_l, const int & node, Vector & vel) const;
   void UpdateNodePosition(Vector &S, const int & node, const Vector &x) const;
   void GetNodePositionFromBV(const Vector &S, const int & node, Vector & x) const;
   void GetNodePosition(const ParGridFunction &x_gf, const int & node, Vector & x) const;

   /* Dofs */
   void GetFaceDofs(const int &face, Array<int> &face_dofs) const;
   void VertexGetAdjacentVertices(const int &vertex, Array<int> &adj_vertices) const;
   void VertexGetAdjacentFaces(const int &vertex, Array<int> &adj_faces) const;
   int GetNumAdjFaces(const int &vertex) const;

   /* Rotate a vector */
   void Orthogonal(Vector &v) const;
   void Perpendicular(Vector &v) const;

}; // End class Geometric

} // ns hydrodynamic

} // ns mfem

#endif // GEOMETRIC