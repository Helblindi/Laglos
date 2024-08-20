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
   ParFiniteElementSpace &H1;
   ParMesh *pmesh;
   int NDofs_H1;
   int NVDofs_H1;

public:
   Geometric(Array<int> offsets, ParFiniteElementSpace &h1);
   ~Geometric() {}

   inline Array<int> GetBlockOffsets() const { return block_offsets; }
   inline ParFiniteElementSpace GetPFES() const { return H1; }
   inline ParMesh *GetParMesh() const { return pmesh; }
   inline int GetNDofs_H1() const { return NDofs_H1; }

   void UpdateNodeVelocity(Vector &S, const int & node, const Vector & vel) const;
   void UpdateNodeVelocity(ParGridFunction &mv_gf, const int &node, const Vector &vel) const;
   void GetNodeVelocity(const Vector &S, const int & node, Vector & vel) const;
   void GetNodeVelocity(const ParGridFunction &mv_gf, const int & node, Vector & vel) const;
   void GetNodeVelocityVecL(const Vector &mv_gf, const int & node, Vector & vel) const;
   void UpdateNodePosition(Vector &S, const int & node, const Vector &x) const;
   void GetNodePositionFromBV(const Vector &S, const int & node, Vector & x) const;
   void GetNodePosition(const ParGridFunction &x_gf, const int & node, Vector & x) const;

   /* Rotate a vector */
   void Orthogonal(Vector &v) const;
   void Perpendicular(Vector &v) const;

}; // End class Geometric

} // ns hydrodynamic

} // ns mfem

#endif // GEOMETRIC