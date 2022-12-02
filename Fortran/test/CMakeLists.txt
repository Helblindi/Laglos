cmake_minimum_required(VERSION 3.13)  # CMake version check
project(my-project LANGUAGES CXX Fortran)

enable_language(C)
enable_language(CXX)
enable_language(Fortran)
include(FortranCInterface)
# Check if the Fortran and C++ compilers are able to communicate
FortranCInterface_VERIFY(CXX)

set (SOURCE_FILES test.cpp ../Lagrangian_arbitrary_eos.f90)

add_executable(example ${SOURCE_FILES})
