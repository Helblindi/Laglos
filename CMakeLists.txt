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
set(${PROJECT_NAME}_VERSION 1.0)
string(TOUPPER "${PROJECT_NAME}" PROJECT_NAME_UC)

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

# CMake needs to know where to find things
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)

# Load Laglos CMake utilities
include(LaglosCmakeUtilities)

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
if (LAGLOS_USE_HIOP)
  include(LaglosFindHIOP)
  list(APPEND _laglos_libraries HiOp::hiop_tpl)
endif()

message(STATUS "Laglos build type: CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
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

# List all possible libraries in order of dependencies.
# [METIS < SuiteSparse]:
#    With newer versions of SuiteSparse which include METIS header using 64-bit
#    integers, the METIS header (with 32-bit indices, as used by mfem) needs to
#    be before SuiteSparse.
set(LAGLOS_TPLS mfem BOOST COINHSL HIOP)

# Add all *_FOUND libraries in the variable TPL_LIBRARIES.
set(TPL_LIBRARIES "")
set(TPL_INCLUDE_DIRS "")
foreach(TPL IN LISTS LAGLOS_TPLS)
  if (${TPL}_FOUND)
    message(STATUS "LAGLOS: using package ${TPL}")
    list(APPEND TPL_LIBRARIES ${${TPL}_LIBRARIES})
    list(APPEND TPL_INCLUDE_DIRS ${${TPL}_INCLUDE_DIRS})
  endif()
endforeach(TPL)
list(REVERSE TPL_LIBRARIES)
list(REMOVE_DUPLICATES TPL_LIBRARIES)
list(REVERSE TPL_LIBRARIES)
list(REMOVE_DUPLICATES TPL_INCLUDE_DIRS)
message(STATUS "TPL_INCLUDE_DIRS = ${TPL_INCLUDE_DIRS}")

#-------------------------------------------------------------------------------
# Define and configure the MFEM library
#-------------------------------------------------------------------------------

# Headers and sources
set(SOURCES "")
set(HEADERS "")
set(LAGLOS_SOURCE_DIRS include/ problems/) # appending / will copy the contents of this directory
set(LAGLOS_FORTRAN_LIBS "")
include_directories(include problems ${CMAKE_BINARY_DIR})
add_subdirectory(Fortran)
add_subdirectory(src)
# Add in Fortran libraries
list(APPEND _laglos_libraries ${LAGLOS_FORTRAN_LIBS})

MESSAGE(STATUS "LAGLOS_FORTRAN_LIBS: ${LAGLOS_FORTRAN_LIBS}")

#-------------------------------------------------------------------------------
# Installation
#-------------------------------------------------------------------------------
include(GenerateExportHeader)
set(INSTALL_INCLUDE_DIR include
  CACHE PATH "Relative path for installing header files.")
set(INSTALL_BIN_DIR bin
  CACHE PATH "Relative path for installing the binaries.")
set(INSTALL_LIB_DIR lib
  CACHE PATH "Relative path for installing the library.")
set(INSTALL_CMAKE_DIR lib/cmake
  CACHE PATH "Relative path for installing cmake config files.")
add_library(Laglos SHARED ${SOURCES} ${HEADERS})
set_property(TARGET Laglos PROPERTY VERSION ${Laglos_VERSION})

target_compile_definitions(Laglos PUBLIC "${_laglos_compile_defs}")
target_compile_options(Laglos PUBLIC "${_laglos_compile_opts}")
target_include_directories(Laglos SYSTEM PUBLIC "${_laglos_include_dirs}")
target_link_libraries(Laglos PUBLIC "${_laglos_libraries}")

install(TARGETS ${PROJECT_NAME} comp_lambda_lagrangian_greedy comp_lambda_lagrangian comp_lambda_eulerian vdw_fortran
   EXPORT ${PROJECT_NAME_UC}Targets
   RUNTIME DESTINATION ${INSTALL_BIN_DIR}
   LIBRARY DESTINATION ${INSTALL_LIB_DIR}
   ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
   INCLUDES DESTINATION ${INSTALL_INCLUDE_DIR}
)

install(DIRECTORY ${LAGLOS_SOURCE_DIRS}
  DESTINATION ${INSTALL_INCLUDE_DIR}
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")

# Package the whole thing up nicely
include(CMakePackageConfigHelpers)

export(TARGETS ${PROJECT_NAME} comp_lambda_lagrangian_greedy comp_lambda_lagrangian comp_lambda_eulerian vdw_fortran
   FILE "${PROJECT_BINARY_DIR}/LaglosTargets.cmake")
export(PACKAGE ${PROJECT_NAME})

###
# version to be installed
set(INCLUDE_INSTALL_DIRS ${INSTALL_INCLUDE_DIR}  ${TPL_INCLUDE_DIRS})
set(LIB_INSTALL_DIR ${INSTALL_LIB_DIR})
configure_package_config_file(CMake/LAGLOSConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/LAGLOSConfig.cmake
  INSTALL_DESTINATION ${INSTALL_CMAKE_DIR}
  PATH_VARS INCLUDE_INSTALL_DIRS LIB_INSTALL_DIR)

# Write the version file (same for build and install tree)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/LAGLOSConfigVersion.cmake
  VERSION ${${PROJECT_NAME}_VERSION}
  COMPATIBILITY AnyNewerVersion )

# Install the config files
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/LAGLOSConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/LAGLOSConfigVersion.cmake
  DESTINATION ${INSTALL_CMAKE_DIR})

# Install the export set for use with the install-tree
install(EXPORT ${PROJECT_NAME_UC}Targets
   FILE ${PROJECT_NAME}Targets.cmake
   DESTINATION ${INSTALL_CMAKE_DIR})


###


file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/results/convergence/temp_output")
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/results/state_vectors")

set(LAGLOS_DIR "${PROJECT_SOURCE_DIR}/")
configure_file("${PROJECT_SOURCE_DIR}/config/var-config.h.in" "${PROJECT_BINARY_DIR}/var-config.h")
configure_file("${PROJECT_SOURCE_DIR}/include/compile_time_vals_template.h" "${PROJECT_SOURCE_DIR}/include/compile_time_vals.h" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/include/hiop_fr_template.options" "${PROJECT_BINARY_DIR}/hiop_fr.options" COPYONLY)

# Finally, add the executable
add_executable(laglos-exe src/laglos.cpp)
set_target_properties(laglos-exe PROPERTIES OUTPUT_NAME Laglos)
target_link_libraries(laglos-exe PRIVATE Laglos)

# Enable testing
# if(ENABLE_TESTS)
enable_testing()
add_subdirectory(tests)
# endif(ENABLE_TESTS)
