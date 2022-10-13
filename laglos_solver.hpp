#ifndef LAGLOS_SOLVER
#define LAGLOS_SOLVER

#include "mfem.hpp"

namespace mfem
{

//
class LagrangianLOOperator : public TimeDependentOperator
{
protected:
   ParFiniteElementSpace &H1, &L2, &L2V;
   ParMesh *pmesh;
   // FE spaces local and global sizes
   const int Vsize_H1;
   const int TVSize_H1;
   const HYPRE_Int GTVSize_H1;
   const int NDofs_H1;

   const int Vsize_L2;
   const int TVSize_L2;
   const HYPRE_Int GTVSize_L2;
   const int NDofs_L2;

   const int Vsize_L2V;
   const int TVSize_L2V;
   const HYPRE_Int GTVSize_L2V;
   const int NDofs_L2V;

   Array<int> block_offsets;
   // Reference to the current mesh configuration.
   mutable ParGridFunction x_gf;
   // const Array<int> &ess_tdofs;
   int dim, NE, l2dofs_cnt, h1dofs_cnt, l2vdofs_cnt;

public:
   LagrangianLOOperator(const int size,
                        ParFiniteElementSpace &h1,
                        ParFiniteElementSpace &l2,
                        ParFiniteElementSpace &l2v);
   ~LagrangianLOOperator();

   virtual void Mult(const Vector &S, Vector &S_new) const;

   void GetCellStateVector(const Vector &S, const int cell, Vector &U);

   void SetCellStateVector(Vector &S_new, const int cell, const Vector &U);

   double ComputeViscosity();

};

}

#endif // LAGLOS_SOLVER