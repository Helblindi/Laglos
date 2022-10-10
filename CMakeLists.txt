cmake_minimum_required(VERSION 3.10)

# set the project name
project(Laglos LANGUAGES C CXX Fortran)

# Verify C++ and Fortran can communicate
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

# Load user settings before the defaults - this way the defaults will not
# overwrite the user set options. If the user has not set all options, we still
# have the defaults.
set(USER_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/config/user.cmake" CACHE PATH
  "Path to optional user configuration file.")
message(STATUS "(optional) USER_CONFIG = ${USER_CONFIG}")
include("${USER_CONFIG}" OPTIONAL)
include("${CMAKE_CURRENT_SOURCE_DIR}/config/defaults.cmake")

# Allow overwriting of the compiler by setting CXX/MPICXX on the command line or
# in user.cmake.
# if (NOT CMAKE_CXX_COMPILER)
#   if (CXX)
#     set(CMAKE_CXX_COMPILER ${CXX})
#     # Avoid some issues when CXX is defined
#     unset(CXX)
#     unset(CXX CACHE)
#   endif()
#   if (MPICXX)
#     # In parallel MPICXX takes precedence, if defined.
#     set(CMAKE_CXX_COMPILER ${MPICXX})
#     # Setting the variables below circumvents autodetection, see FindMPI.cmake.
#     set(MPI_CXX_INCLUDE_PATH "")
#     set(MPI_CXX_LIBRARIES "")
#   endif()
# endif()

# Import MFEM. The following variables can be used to help CMake find MFEM:
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
#  * mfem_DIR - absolute path to where MFEMConfig.cmake is.
message(STATUS "Looking for mfem ...")
set(MFEM_DIR "" CACHE PATH "Path to the MFEM build or install prefix.")
if (MFEM_DIR)
   # find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}"
   #              "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)
   find_package(mfem REQUIRED)
else()
   find_package(mfem REQUIRED NAMES MFEM)
endif()
message(STATUS "Found mfem config in: ${mfem_DIR} (version ${MFEM_VERSION})")
# Use the same C++ compiler as MFEM. This is needed when MFEM was built using
# an MPI wrapper and we do not have explicitly the MPI compile and link flags.
if (NOT CMAKE_CXX_COMPILER AND MFEM_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER "${MFEM_CXX_COMPILER}")
endif()

# Include paths and libraries needed by MFEM
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MFEM_CXX_FLAGS}")
list(APPEND LAGLOS_INCLUDE_DIRS "${MFEM_INCLUDE_DIRS}")
list(APPEND LAGLOS_INCLUDE_DIRS "${MFEM_LIBRARIES}")


# add the executable
add_executable(Laglos laglos.cpp)

# Add fortran
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/Fortran")
add_library(comp_lambda STATIC "${CMAKE_CURRENT_SOURCE_DIR}/Fortran/Lagrangian_arbitrary_eos.f90")

# Link MFEM and Fortran libraries
target_link_libraries(Laglos "${MFEM_LIBRARIES}" comp_lambda)

# Specify Include directories for target
target_include_directories(Laglos PUBLIC ${MFEM_DIR} "${MFEM_INCLUDE_DIRS}")
