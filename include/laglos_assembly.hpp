#ifndef MFEM_LAGLOS_ASSEMBLY
#define MFEM_LAGLOS_ASSEMBLY

#include "mfem.hpp"

namespace mfem
{

// Container for all data needed at quadrature points.
struct QuadratureData
{
   // Reference to physical Jacobian for the initial mesh.
   // These are computed only at time zero and stored here.
   DenseTensor Jac0inv;

   // Quadrature data used for full/partial assembly of the force operator.
   // At each quadrature point, it combines the stress, inverse Jacobian,
   // determinant of the Jacobian and the integration weight.
   // It must be recomputed in every time step.
   DenseTensor stressJinvT;

   // Quadrature data used for full/partial assembly of the mass matrices.
   // At time zero, we compute and store (rho0 * det(J0) * qp_weight) at each
   // quadrature point. Note the at any other time, we can compute
   // rho = rho0 * det(J0) / det(J), representing the notion of pointwise mass
   // conservation.
   Vector rho0DetJ0w;

   // Initial length scale. This represents a notion of local mesh size.
   // We assume that all initial zones have similar size.
   double h0;

   // Estimate of the minimum time step over all quadrature points. This is
   // recomputed at every time step to achieve adaptive time stepping.
   double dt_est;

   QuadratureData(int dim, int NE, int quads_per_el)
      : Jac0inv(dim, dim, NE * quads_per_el),
        stressJinvT(NE * quads_per_el, dim, dim),
        rho0DetJ0w(NE * quads_per_el) { }
};

// This class is used only for visualization. It assembles (rho, phi) in each
// zone, which is used by LagrangianHydroOperator::ComputeDensity to do an L2
// projection of the density.
class DensityIntegrator : public LinearFormIntegrator
{
   using LinearFormIntegrator::AssembleRHSElementVect;
private:
   const QuadratureData &qdata;

public:
   DensityIntegrator(QuadratureData &qdata) : qdata(qdata) { }
   virtual void AssembleRHSElementVect(const FiniteElement &fe,
                                       ElementTransformation &Tr,
                                       Vector &elvect);
};

}
#endif // MFEM_LAGLOS_ASSEMBLY