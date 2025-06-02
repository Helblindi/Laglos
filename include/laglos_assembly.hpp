#ifndef LAGLOS_ASSEMBLY
#define LAGLOS_ASSEMBLY

#include "mfem.hpp"

using namespace std;

namespace mfem
{

namespace hydroLO
{

class DGNormalIntegrator : public BilinearFormIntegrator
{
private:
   real_t coeff;
   int xi; // Represents coordinate of normal
   Vector shape1, shape2;

public:
DGNormalIntegrator(const real_t &_coeff, const int &i) : coeff(_coeff), xi(i) {}
   virtual void AssembleFaceMatrix(const FiniteElement &el1,
                                   const FiniteElement &el2,
                                   FaceElementTransformations &trans,
                                   DenseMatrix &elmat);
};

}// end ns hydroLO

}// end ns mfem
#endif // LAGLOS_ASSEMBLY