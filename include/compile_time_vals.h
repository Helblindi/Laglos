#ifndef COMPILE_TIME_VALS
#define COMPILE_TIME_VALS

namespace mfem
{
namespace hydrodynamics
{
namespace CompileTimeVals
{
   /* Problem specific */
   const static int dim = 1;

   /**
    * Directory specific information
    * Note: Can change this to a local path that you would like results output to. 
    * */
   const static std::string results_dir = std::string(LAGLOS_DIR) + "build/results/";
}
}
}

#endif // COMPILE_TIME_VALS