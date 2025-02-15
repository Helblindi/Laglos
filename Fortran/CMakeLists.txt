cmake_minimum_required(VERSION 3.10)

# Add Fortran library (and link)
add_library(comp_lambda_lagrangian STATIC "lagrangian-mws/lagrangian_arbitrary_eos.f90")
add_library(comp_lambda_lagrangian_greedy STATIC "lagrangian-mws/greedy_viscosity_Lagrangian_for_arbitrary_eos.f90")
add_library(comp_lambda_eulerian STATIC "eulerian-mws/eulerian_arbitrary_eos.f90")
add_library(vdw_fortran STATIC "vdw/vdw.f90")

set(LAGLOS_FORTRAN_LIBS ${LAGLOS_FORTRAN_LIBS}
   comp_lambda_eulerian
   comp_lambda_lagrangian
   comp_lambda_lagrangian_greedy
   vdw_fortran
   PARENT_SCOPE)

install(TARGETS comp_lambda_lagrangian_greedy comp_lambda_lagrangian comp_lambda_eulerian vdw_fortran
        EXPORT LaglosFortranTargets
        LIBRARY DESTINATION ${INSTALL_LIB_DIR}
        ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
        INCLUDES DESTINATION ${INSTALL_INCLUDE_DIR})
