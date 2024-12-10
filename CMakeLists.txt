cmake_minimum_required(VERSION 3.10)
message(STATUS "CMake version: ${CMAKE_VERSION}")
set(USER_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/config/user.cmake" CACHE PATH
  "Path to optional user configuration file.")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Load user settings before the defaults - this way the defaults will not
# overwrite the user set options. If the user has not set all options, we still
# have the defaults.
include("${USER_CONFIG}" OPTIONAL RESULT_VARIABLE USER_CONFIG_LOADED)
if (USER_CONFIG_LOADED)
  set(USER_CONFIG_LOADED "LOADED")
endif()
message(STATUS "Loading USER_CONFIG = ${USER_CONFIG} (${USER_CONFIG_LOADED})")
include("${CMAKE_CURRENT_SOURCE_DIR}/config/defaults.cmake")

# Prohibit in-source builds
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are prohibited.")
endif ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

project(Laglos LANGUAGES C CXX Fortran)

# Import MFEM. The following variables can be used to help CMake find MFEM:
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
#  * mfem_DIR - absolute path to where MFEMConfig.cmake is.
message(STATUS "Looking for mfem ...")
if (MFEM_DIR)
   find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}"
                "${MFEM_DIR}/build" NO_DEFAULT_PATH)
else()
   find_package(mfem REQUIRED NAMES MFEM)
endif()
message(STATUS "Found mfem config in: ${mfem_DIR} (version ${MFEM_VERSION})")
# Use the same C++ compiler as MFEM. This is needed when MFEM was built using
# an MPI wrapper and we do not have explicitly the MPI compile and link flags.
if (NOT CMAKE_CXX_COMPILER AND MFEM_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER "${MFEM_CXX_COMPILER}")
endif()

# Default options match the Makefile
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()


#
# Start finding everything
#
set(_laglos_compile_defs)
set(_laglos_compile_opts)
set(_laglos_include_dirs)
set(_laglos_libraries)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)


# if (CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
#   list(APPEND _glvis_compile_defs "GLVIS_DEBUG")
#   list(APPEND _glvis_compile_opts "-Wall")
# endif()

# Include paths and libraries needed by MFEM
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MFEM_CXX_FLAGS}")
list(APPEND _laglos_include_dirs "${MFEM_INCLUDE_DIRS}")
list(APPEND _laglos_libraries "${MFEM_LIBRARIES}")
if (LAGLOS_USE_HIOP)
  include(LaglosFindHIOP)
  list(APPEND _laglos_libraries HiOp::hiop_tpl)
endif()

message(STATUS "Laglos build type: CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
# message(STATUS "GLVis opts: ${_glvis_compile_opts}")
# message(STATUS "GLVis include dirs: ${_glvis_include_dirs}")
# message(STATUS "GLVis libraries: ${_glvis_libraries}")
message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")

# Boost required for filesystem modification
# Avoid setting CMP0167 if using a CMake version that predates 3.26
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost REQUIRED COMPONENTS filesystem)
# Include libraries needed by lboost_filesystem
# https://stackoverflow.com/questions/9235679/create-a-directory-if-it-doesnt-exist
list(APPEND _laglos_libraries ${Boost_FILESYSTEM_LIBRARY})

# Add in Fortran libraries
list(APPEND _laglos_libraries comp_lambda_lagrangian_greedy comp_lambda_lagrangian comp_lambda_eulerian vdw_fortran)

# Add in quick fix to sparse stuff
# list(APPEND _laglos_libraries "/opt/homebrew/Cellar/gcc/14.1.0_2/lib/gcc/14/libgfortran.a" "/opt/homebrew/Cellar/gcc/14.1.0_2/lib/gcc/14/libquadmath.a" "/opt/homebrew/Cellar/gcc/14.1.0_2/lib/gcc/14/gcc/aarch64-apple-darwin23/14/libgcc.a")

#
# Setup the GLVis library target
#
include_directories(include problems ${CMAKE_BINARY_DIR})
add_subdirectory(Fortran)
add_subdirectory(src)
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/results/convergence/temp_output")
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/results/state_vectors")

set(LAGLOS_DIR "${PROJECT_SOURCE_DIR}/")
configure_file("${PROJECT_SOURCE_DIR}/config/var-config.h.in" "${PROJECT_BINARY_DIR}/var-config.h")
configure_file("${PROJECT_SOURCE_DIR}/include/compile_time_vals_template.h" "${PROJECT_SOURCE_DIR}/include/compile_time_vals.h" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/include/hiop_fr_template.options" "${PROJECT_BINARY_DIR}/hiop_fr.options" COPYONLY)

# Finally, add the executable
add_executable(laglos-exe src/laglos.cpp)
set_target_properties(laglos-exe PROPERTIES OUTPUT_NAME laglos)
target_link_libraries(laglos-exe PRIVATE laglos)

# Enable testing
# if(ENABLE_TESTS)
enable_testing()
add_subdirectory(tests)
# endif(ENABLE_TESTS)
