cmake_minimum_required(VERSION 3.10)

include_directories(${TEST_SRC_DIR}/src
                    ${MFEM_INCLUDE_DIRS}
                    ${MFEM_DIR})

# TESTING
add_executable(TestInstantiator EXCLUDE_FROM_ALL testing.cpp "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(TestInstantiator PRIVATE "${MFEM_LIBRARIES}" comp_lambda_lagrangian comp_lambda_eulerian)
add_test(NAME TestInstantiator COMMAND TestInstantiator)

add_executable(TestRTvel EXCLUDE_FROM_ALL test_RT_vel.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(TestRTvel PRIVATE "${MFEM_LIBRARIES}" comp_lambda_lagrangian comp_lambda_eulerian)
add_test(NAME TestRTvel COMMAND TestRTvel)

add_executable(Test1D EXCLUDE_FROM_ALL test_1D.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(Test1D PRIVATE "${MFEM_LIBRARIES}" comp_lambda_lagrangian comp_lambda_eulerian)
add_test(NAME Test1D COMMAND Test1D)

add_executable(Test2D EXCLUDE_FROM_ALL test_2D.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(Test2D PRIVATE "${MFEM_LIBRARIES}" comp_lambda_lagrangian comp_lambda_eulerian)
add_test(NAME Test2D COMMAND Test2D)

# Validate Fortran Max Wave Speed calculations
add_executable(TestEulerianMWS EXCLUDE_FROM_ALL test_fortran_eulerian_mws.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(TestEulerianMWS PRIVATE "${MFEM_LIBRARIES}" comp_lambda_eulerian)
add_test(NAME TestEulerianMWS
         COMMAND TestEulerianMWS "${PROJECT_SOURCE_DIR}/build/tests/eulerian-mws-output.txt")
add_test(NAME EulerMWSCompare
         COMMAND ${CMAKE_COMMAND} -E compare_files
               "${PROJECT_SOURCE_DIR}/build/tests/eulerian-mws-output.txt"
               "${PROJECT_SOURCE_DIR}/Fortran/eulerian-mws/output_ref_c")

# TODO: Lagrangian Fortran Subroutine still does not return correct values
add_executable(TestLagrangianMWS EXCLUDE_FROM_ALL test_fortran_lagrangian_mws.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(TestLagrangianMWS PRIVATE "${MFEM_LIBRARIES}" comp_lambda_lagrangian)
add_test(NAME TestLagrangianMWS
         COMMAND TestLagrangianMWS "${PROJECT_SOURCE_DIR}/build/tests/lagrangian-mws-output.txt")
add_test(NAME LagrangeMWSCompare
         COMMAND ${CMAKE_COMMAND} -E compare_files
               "${PROJECT_SOURCE_DIR}/build/tests/lagrangian-mws-output.txt"
               "${PROJECT_SOURCE_DIR}/Fortran/lagrangian-mws/output_ref_c")


# Create a new target to build and run tests
# This ensure make all doesn't make every test at each compilation stage
add_custom_target(build_tests)
add_custom_target(run_tests COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(run_tests build_tests)
add_dependencies(build_tests TestInstantiator TestRTvel Test1D Test2D TestEulerianMWS TestLagrangianMWS)