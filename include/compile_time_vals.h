#ifndef COMPILE_TIME_VALS
#define COMPILE_TIME_VALS

namespace mfem
{
namespace hydrodynamics
{
class CompileTimeVals
{
public:
/*
*     Case 0: Constant density and specific internal energy. 0 velocity.
*             In this case, we shouldn't see any values changing from
*             the prescribed initial conditions.
*     Case 1: The above but with normal vector with equal components for
*             the cell velocity.
*     Case 2: Isentropic vortex as described in (6.1) with a stationary 
*             center at the origin.
*     Case 3: Isentropic vortex (see above) with moving center.
*     Case 4: Noh problem as described in (6.3)
*     Case 5: 1D Horizontal Movement on 2D mesh
*     Case 6: Shocktubes 1 - Sod, 2 - Lax, 3 - Leblanc
*     Case 7: Saltzman problem. See https://people.tamu.edu/~guermond/PUBLICATIONS/guermond_popov_Saavedra_JCP_2020.pdf.
*             Requires Neumann BC on right face and Dirichlet elsewhere.
*     Case 8: Linear velocity field to validate Least Squares method
*/
   const static int problem = 6;
   const static int dim = 2;
   
   /* Various parameters */
   const static bool distort_mesh = false; // Relevant in problems 7
   const static int shocktube = 1; // 1 - Sod, 2 - Lax, 3 - Leblanc
   constexpr static double rotation_angle = 0; // 0 - 1D horizontal velocity

   /* Saltzman Problem parameters */
   const static bool change_CFL = true;
   constexpr static double CFL_first = 0.01;
   constexpr static double CFL_second = 0.25;
   constexpr static double CFL_time_change = 0.01; // From Boscheri's paper
};
}
}

#endif // COMPILE_TIME_VALS