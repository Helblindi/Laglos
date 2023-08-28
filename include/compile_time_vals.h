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
   const static int dim = 2;
   
   /* Various parameters */
   const static int shocktube = 1; // 1 - Sod, 2 - Lax, 3 - Leblanc
   constexpr static double rotation_angle = 0; // 0 - 1D horizontal velocity   
};
}
}

#endif // COMPILE_TIME_VALS