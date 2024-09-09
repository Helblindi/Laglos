cmake_minimum_required(VERSION 3.10)
message(STATUS "CMake version: ${CMAKE_VERSION}")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LAGLOS_USE_HIOP_SPARSE "Utilize sparse linear algebra" ON)

# Prohibit in-source builds
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are prohibited.")
endif ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

project(Laglos LANGUAGES C CXX Fortran)

# Import MFEM. The following variables can be used to help CMake find MFEM:
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
#  * mfem_DIR - absolute path to where MFEMConfig.cmake is.
message(STATUS "Looking for mfem ...")
set(MFEM_DIR "" CACHE PATH "Path to the MFEM build or install prefix.")
if (MFEM_DIR)
   find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}"
                "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)
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
  set(CMAKE_BUILD_TYPE "Release")
endif()


#
# Start finding everything
#

set(_laglos_compile_defs)
set(_laglos_compile_opts)
set(_laglos_include_dirs)
set(_laglos_libraries)


# if (CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
#   list(APPEND _glvis_compile_defs "GLVIS_DEBUG")
#   list(APPEND _glvis_compile_opts "-Wall")
# endif()

# Include paths and libraries needed by MFEM
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MFEM_CXX_FLAGS}")
list(APPEND _laglos_include_dirs "${MFEM_INCLUDE_DIRS}")
list(APPEND _laglos_libraries "${MFEM_LIBRARIES}")

message(STATUS "Laglos build type: CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message(STATUS "Laglos defines: ${_glvis_compile_defs}")
# message(STATUS "GLVis opts: ${_glvis_compile_opts}")
# message(STATUS "GLVis include dirs: ${_glvis_include_dirs}")
# message(STATUS "GLVis libraries: ${_glvis_libraries}")
message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")

# Boost required for filesystem modification
cmake_policy(SET CMP0167 NEW)
find_package(Boost REQUIRED COMPONENTS filesystem)
# Include libraries needed by lboost_filesystem
# https://stackoverflow.com/questions/9235679/create-a-directory-if-it-doesnt-exist
list(APPEND _laglos_libraries ${Boost_FILESYSTEM_LIBRARY})

# Add in Fortran libraries
list(APPEND _laglos_libraries comp_lambda_lagrangian comp_lambda_eulerian vdw_fortran)

# Add in quick fix to sparse stuff
list(APPEND _laglos_libraries "${PROJECT_SOURCE_DIR}/../coinhsl/lib/libcoinhsl.a" "/opt/homebrew/Cellar/gcc/14.1.0_2/lib/gcc/14/libgfortran.a" "/opt/homebrew/Cellar/gcc/14.1.0_2/lib/gcc/14/libquadmath.a" "/opt/homebrew/Cellar/gcc/14.1.0_2/lib/gcc/14/gcc/aarch64-apple-darwin23/14/libgcc.a")

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

# Finally, add the executable
add_executable(laglos-exe src/laglos.cpp)
set_target_properties(laglos-exe PROPERTIES OUTPUT_NAME laglos)
target_link_libraries(laglos-exe PRIVATE laglos)

# Enable testing
# if(ENABLE_TESTS)
enable_testing()
add_subdirectory(tests)
# endif(ENABLE_TESTS)
