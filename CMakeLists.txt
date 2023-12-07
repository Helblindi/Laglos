cmake_minimum_required(VERSION 3.10)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE) # Look into this when have time.

# Require C++11 and disable compiler-specific extensions
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(HOMEBREW_PREFIX "/opt/homebrew"
  CACHE PATH "Path to Homebrew installation")

# set(CMAKE_C_COMPILER "${HOMEBREW_PREFIX}/bin/gcc-12")
# set(CMAKE_CXX_COMPILER "${HOMEBREW_PREFIX}/bin/g++-12")

# edits
if (NOT CMAKE_CXX_COMPILER)
  if (CXX)
    set(CMAKE_CXX_COMPILER ${CXX})
    # Avoid some issues when CXX is defined
    unset(CXX)
    unset(CXX CACHE)
  endif()
endif()

# set the project name
project(Laglos LANGUAGES C CXX Fortran)

# Prevent in source builds
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
include(PreventInSourceBuilds)

# Verify C++ and Fortran can communicate
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

# Boost required for filesystem modification
find_package(Boost REQUIRED COMPONENTS filesystem)

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

# Include paths and libraries needed by MFEM
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MFEM_CXX_FLAGS}")
list(APPEND LAGLOS_INCLUDE_DIRS "${MFEM_INCLUDE_DIRS}")
list(APPEND LAGLOS_INCLUDE_DIRS "${MFEM_LIBRARIES}")

# Set LAGLOS_DIR as a directory location that can be referred to in the C++ files using a configure file
# https://cmake.org/cmake/help/latest/command/configure_file.html
set(LAGLOS_DIR "${PROJECT_SOURCE_DIR}/")
configure_file("${PROJECT_SOURCE_DIR}/config/var-config.h.in" "${PROJECT_BINARY_DIR}/var-config.h")

# Create compile_time_vals.h file that can be changed but ignored by git
configure_file("${PROJECT_SOURCE_DIR}/include/compile_time_vals_template.h" "${PROJECT_SOURCE_DIR}/include/compile_time_vals.h" COPYONLY)

### Include Directories
INCLUDE_DIRECTORIES(include problems)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})

# populate matplotlib repository
include(FetchContent)
FetchContent_Declare(
    matplotlib
    GIT_REPOSITORY https://github.com/lava/matplotlib-cpp.git
    GIT_TAG        f23347fca25219d1c42cbb91608b5556814bf572
)
FetchContent_GetProperties(matplotlib)
if(NOT matplotlib_POPULATED)
    FetchContent_Populate(matplotlib)
endif()
include_directories(SYSTEM ${matplotlib_SOURCE_DIR})

add_subdirectory(Fortran)
add_subdirectory(src)
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/results/convergence/temp_output")
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/results/state_vectors")

# Enable testing
# https://cmake.org/cmake/help/latest/command/add_test.html
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

# Add test folder
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()