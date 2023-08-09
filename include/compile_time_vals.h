#ifndef COMPILE_TIME_VALS
#define COMPILE_TIME_VALS

namespace mfem
{
namespace hydrodynamics
{
class CompileTimeVals
{
public:
   /* See initial_vals.hpp for test problems */
   // const static int problem = 1;
   const static int dim = 1;
   
   /* Various parameters */
   const static bool distort_mesh = false; // Relevant in problems 7
   const static int shocktube = 1; // 1 - Sod, 2 - Lax, 3 - Leblanc
   constexpr static double rotation_angle = 0; // 0 - 1D horizontal velocity

   /* Saltzman Problem parameters */
   const static bool change_CFL = false;
   constexpr static double CFL_first = 0.01;
   constexpr static double CFL_second = 0.25;
   constexpr static double CFL_time_change = 0.01; // From Boscheri's paper
};
}
}

#endif // COMPILE_TIME_VALS