cmake_minimum_required(VERSION 3.10)

include_directories(${TEST_SRC_DIR}/src
                    ${MFEM_INCLUDE_DIRS}
                    ${MFEM_DIR})

# TESTING
# add_executable(TestRTvel EXCLUDE_FROM_ALL test_RT_vel.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
# target_link_libraries(TestRTvel PRIVATE "${MFEM_LIBRARIES}" comp_lambda_lagrangian comp_lambda_eulerian)
# add_test(NAME TestRTvel COMMAND TestRTvel)

# Raviart-Thomas Velocity Construction
# add_executable(TestRaviart EXCLUDE_FROM_ALL test_Raviart_vel.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
# target_link_libraries(TestRaviart PRIVATE "${MFEM_LIBRARIES}" comp_lambda_lagrangian comp_lambda_eulerian)
# add_test(NAME TestRaviart COMMAND TestRaviart)

# add_executable(Test1D EXCLUDE_FROM_ALL test_1D.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
# target_link_libraries(Test1D PRIVATE "${MFEM_LIBRARIES}" comp_lambda_lagrangian comp_lambda_eulerian)
# add_test(NAME Test1D COMMAND Test1D)

# add_executable(Test2D EXCLUDE_FROM_ALL test_2D.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
# target_link_libraries(Test2D PRIVATE "${MFEM_LIBRARIES}" comp_lambda_lagrangian comp_lambda_eulerian)
# add_test(NAME Test2D COMMAND Test2D)

# add_executable(TestMM EXCLUDE_FROM_ALL test_mm.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
# target_link_libraries(TestMM PRIVATE "${MFEM_LIBRARIES}" comp_lambda_lagrangian comp_lambda_eulerian vdw_fortran)
# add_test(NAME TestMM COMMAND TestMM)

# Validate Fortran Max Wave Speed calculations
add_executable(TestEulerianMWS EXCLUDE_FROM_ALL test_fortran_eulerian_mws.cpp "${MFEM_DIR}")
target_link_libraries(TestEulerianMWS PRIVATE "${MFEM_LIBRARIES}" comp_lambda_eulerian)
add_test(NAME TestEulerianMWS
         COMMAND TestEulerianMWS "${PROJECT_SOURCE_DIR}/build/tests/eulerian-mws-output.txt")
add_test(NAME EulerMWSCompare
         COMMAND ${CMAKE_COMMAND} -E compare_files
               "${PROJECT_SOURCE_DIR}/build/tests/eulerian-mws-output.txt"
               "${PROJECT_SOURCE_DIR}/Fortran/eulerian-mws/output_ref_c")

# TODO: Lagrangian Fortran Subroutine still does not return correct values
add_executable(TestLagrangianMWS EXCLUDE_FROM_ALL test_fortran_lagrangian_mws.cpp)
target_link_libraries(TestLagrangianMWS PRIVATE "${MFEM_LIBRARIES}" comp_lambda_lagrangian)
add_test(NAME TestLagrangianMWS
         COMMAND TestLagrangianMWS "${PROJECT_SOURCE_DIR}/build/tests/lagrangian-mws-output.txt")
add_test(NAME LagrangeMWSCompare
         COMMAND ${CMAKE_COMMAND} -E compare_files
               "${PROJECT_SOURCE_DIR}/build/tests/lagrangian-mws-output.txt"
               "${PROJECT_SOURCE_DIR}/Fortran/lagrangian-mws/output_ref_c")

# Test Lagrange Multipliers implementation
# add_executable(TestLM EXCLUDE_FROM_ALL test_lagrange_multipliers.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/geometry.cpp")
# target_link_libraries(TestLM PRIVATE "${MFEM_LIBRARIES}")
# add_test(NAME TestLM COMMAND TestLM)

# Test Equation of State
add_executable(TestEOS EXCLUDE_FROM_ALL test_eos.cpp)
target_include_directories(TestEOS PRIVATE ${PROJECT_SOURCE_DIR}/include)
add_test(NAME TestEOS COMMAND TestEOS)


### Test programs to compare output of various numerical tests
# Sod 1D
file(COPY test_sod_1d.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_sod_1d.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestSod1D
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_sod_1d.sh ${PROJECT_SOURCE_DIR} "sod1d_output.txt"
)

file(COPY compare_ignore_last_lines.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(NAME CompareSod1D
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/sod1d_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/sod1d.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

# Sod 2D
file(COPY test_sod_2d.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_sod_2d.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestSod2D
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_sod_2d.sh ${PROJECT_SOURCE_DIR} "sod2d_output.txt"
)
add_test(NAME CompareSod2D
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/sod2d_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/sod2d.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

# Sod radial

# Noh
file(COPY test_noh.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_noh.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestNoh
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_noh.sh ${PROJECT_SOURCE_DIR} "noh_output.txt"
)
add_test(NAME CompareNoh
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/noh_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/noh.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

# Sedov
file(COPY test_sedov.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_sedov.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestSedov
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_sedov.sh ${PROJECT_SOURCE_DIR} "sedov_output.txt"
)
add_test(NAME CompareSedov
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/sedov_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/sedov.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

# Saltzman
file(COPY test_saltzman.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_saltzman.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestSaltzman
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_saltzman.sh ${PROJECT_SOURCE_DIR} "saltzman_output.txt"
)
add_test(NAME CompareSaltzman
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/saltzman_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/saltzman.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

# Isentropic Vortex


#### Test elastic functions ####
## Elastic Impact
file(COPY test_elastic_impact.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_elastic_impact.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestElasticImpact
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_elastic_impact.sh ${PROJECT_SOURCE_DIR} "elastic_impact_output.txt"
)

## Elastic Shocktube

# Create a new target to build and run tests
# This ensure make all doesn't make every test at each compilation stage
add_custom_target(build_tests)
add_custom_target(run_tests COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(run_tests build_tests)
add_dependencies(build_tests TestEulerianMWS TestLagrangianMWS TestEOS)

# Custom target to clean tmp testing directory
add_custom_target(clean_tests
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_BINARY_DIR}/tests
)