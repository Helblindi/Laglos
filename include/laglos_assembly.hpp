#ifndef LAGLOS_ASSEMBLY
#define LAGLOS_ASSEMBLY

#include "mfem.hpp"

using namespace std;

namespace mfem
{

namespace hydroLO
{

// Container for all data needed at quadrature points.
struct QuadratureData
{
   // TODO: use QuadratureFunctions?

   // Reference to physical Jacobian for the initial mesh. These are computed
   // only at time zero and stored here.
   DenseTensor Jac0inv;

   // Quadrature data used for full/partial assembly of the force operator. At
   // each quadrature point, it combines the stress, inverse Jacobian,
   // determinant of the Jacobian and the integration weight. It must be
   // recomputed in every time step.
   DenseTensor stressJinvT;

   // Quadrature data used for full/partial assembly of the mass matrices. At
   // time zero, we compute and store (rho0 * det(J0) * qp_weight) at each
   // quadrature point. Note the at any other time, we can compute
   // rho = rho0 * det(J0) / det(J), representing the notion of pointwise mass
   // conservation.
   Vector rho0DetJ0w;

   // Initial length scale. This represents a notion of local mesh size. We
   // assume that all initial zones have similar size.
   double h0;

   // Estimate of the minimum time step over all quadrature points. This is
   // recomputed at every time step to achieve adaptive time stepping.
   double dt_est;

   QuadratureData(int dim, int nzones, int quads_per_zone)
      : Jac0inv(dim, dim, nzones * quads_per_zone),
        stressJinvT(nzones * quads_per_zone, dim, dim),
        rho0DetJ0w(nzones * quads_per_zone) { }

   void Resize(int dim, int nzones, int quads_per_zone)
   {
      Jac0inv.SetSize(dim, dim, nzones * quads_per_zone);
      stressJinvT.SetSize(nzones * quads_per_zone, dim, dim);
      rho0DetJ0w.SetSize(nzones * quads_per_zone);
   }
};


// This class is used only for visualization. It assembles (rho, phi) in each
// zone, which is used by LagrangianHydroOperator::ComputeDensity to do an L2
// projection of the density.
class DensityIntegrator : public LinearFormIntegrator
{
private:
   const QuadratureData &quad_data;

public:
   DensityIntegrator(QuadratureData &quad_data_) : quad_data(quad_data_) { }

   virtual void AssembleRHSElementVect(const FiniteElement &fe,
                                       ElementTransformation &Tr,
                                       Vector &elvect);
};


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