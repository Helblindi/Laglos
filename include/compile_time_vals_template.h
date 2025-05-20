/**
 * DO NOT MODIFY THIS FILE. This is a template file for your systems
 * configuration.  CMake copies this file into include/compile_time_vals.h
 * and edits should be made to that file.  This is a system specific 
 * configuration file to be ignored by git.
*/

#ifndef COMPILE_TIME_VALS
#define COMPILE_TIME_VALS

namespace mfem
{
namespace hydroLO
{
namespace CompileTimeVals
{
   /**
    * Directory specific information
    * Note: Can change this to a local path that you would like results output to. 
    * */
   const static std::string results_dir = std::string(LAGLOS_DIR) + "build/results/";
}
}
}

#endif // COMPILE_TIME_VALS
