cmake_minimum_required(VERSION 3.10)

# Add Fortran library (and link)
add_library(comp_lambda_lagrangian STATIC "${PROJECT_SOURCE_DIR}/Fortran/lagrangian-mws/lagrangian_arbitrary_eos.f90")
add_library(comp_lambda_eulerian STATIC "${PROJECT_SOURCE_DIR}/Fortran/eulerian-mws/eulerian_arbitrary_eos.f90")
add_library(vdw_fortran STATIC "${PROJECT_SOURCE_DIR}/Fortran/vdw/vdw.f90")
