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

add_executable(TestCij EXCLUDE_FROM_ALL test_cij_construction.cpp)
target_link_libraries(TestCij PRIVATE Laglos)
add_test(NAME TestCij COMMAND TestCij)

add_executable(TestLimiter EXCLUDE_FROM_ALL test_limiter.cpp)
target_link_libraries(TestLimiter PRIVATE Laglos)
add_test(NAME TestLimiter COMMAND TestLimiter)

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
# Smooth 1D
file(COPY test_smooth_1d.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_smooth_1d.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestSmooth1D
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_smooth_1d.sh ${PROJECT_SOURCE_DIR} "smooth1d_output.txt"
)

file(COPY compare_ignore_last_lines.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(NAME CompareSmooth1D
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/smooth1d_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/smooth1d.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

# Sod 1D
file(COPY test_sod_1d.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_sod_1d.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestSod1D
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_sod_1d.sh ${PROJECT_SOURCE_DIR} "sod1d_output.txt"
)
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
file(COPY test_iv.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_iv.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestIV
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_iv.sh ${PROJECT_SOURCE_DIR} "iv_output.txt"
)
add_test(NAME CompareIV
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/iv_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/iv.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

# Taylor Green Vortex
file(COPY test_tg.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_tg.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestTG
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_tg.sh ${PROJECT_SOURCE_DIR} "tg_output.txt"
)
add_test(NAME CompareTG
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/tg_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/tg.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

# Kidder shell
file(COPY test_kidder.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_kidder.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestKidder
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_kidder.sh ${PROJECT_SOURCE_DIR} "kidder_output.txt"
)
add_test(NAME CompareKidder
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/kidder_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/kidder.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

# VDW
# 2
file(COPY test_vdw2.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_vdw2.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestVdw2
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_vdw2.sh ${PROJECT_SOURCE_DIR} "vdw2_output.txt"
)
add_test(NAME CompareVdw2
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/vdw2_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/vdw2.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)
# 4
file(COPY test_vdw4.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_vdw4.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestVdw4
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_vdw4.sh ${PROJECT_SOURCE_DIR} "vdw4_output.txt"
)
add_test(NAME CompareVdw4
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/vdw4_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/vdw4.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

### High order tests ###
## Sod 1D
file(COPY test_sod_1d_HO.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_sod_1d_HO.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestSod1DHO
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_sod_1d_HO.sh ${PROJECT_SOURCE_DIR} "HOsod1d.txt"
)
add_test(NAME CompareSod1DHO
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/HOsod1d.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/HOsod1d.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

## Isentropic Vortex ot1
file(COPY test_ot1_iv.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_ot1_iv.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestOT1IV
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ot1_iv.sh ${PROJECT_SOURCE_DIR} "ot1_iv_output.txt"
)
add_test(NAME CompareOT1IV
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/ot1_iv_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/ot1_iv.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

## Taylor Green Vortex ot1
file(COPY test_ot1_tg.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_ot1_tg.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestOT1TG
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ot1_tg.sh ${PROJECT_SOURCE_DIR} "ot1_tg_output.txt"
)
add_test(NAME CompareOT1TG
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/ot1_tg_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/ot1_tg.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

#### Test elastic functions ####
## Elastic Impact
file(COPY test_elastic_impact.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_elastic_impact.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestElasticImpactNH
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_elastic_impact.sh ${PROJECT_SOURCE_DIR} "elastic_impact_output.txt"
)
add_test(NAME CompareElasticImpactNH
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/elastic_impact_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/elastic_impact_NH.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

## Elastic Rotation with MR model
file(COPY test_elastic_rotation_MR.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_elastic_rotation_MR.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestElasticRotationMR
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_elastic_rotation_MR.sh ${PROJECT_SOURCE_DIR} "elastic_rotation_MR_output.txt"
)
add_test(NAME CompareElasticRotationMR
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/elastic_rotation_MR_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/elastic_rotation_MR.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

## Elastic projectile plate with Aortic model
file(COPY test_elastic_pp_aortic.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_elastic_pp_aortic.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
add_test(
    NAME TestElasticPPAortic
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_elastic_pp_aortic.sh ${PROJECT_SOURCE_DIR} "elastic_pp_aortic_output.txt"
)
add_test(NAME CompareElasticPPAortic
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/compare_ignore_last_lines.sh
    ${CMAKE_CURRENT_BINARY_DIR}/elastic_pp_aortic_output.txt
    ${PROJECT_SOURCE_DIR}/tests/output_ref/elastic_pp_aortic.txt
    11 # Ignore timing data at the end of the file as this fluctuates
)

## Elastic Shocktube

# Create a new target to build and run tests
# This ensure make all doesn't make every test at each compilation stage
add_custom_target(build_tests)
add_custom_target(run_tests COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(run_tests build_tests)
add_dependencies(build_tests TestEulerianMWS TestLagrangianMWS TestEOS TestCij TestLimiter)

# Custom target to clean tmp testing directory
add_custom_target(clean_tests
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_BINARY_DIR}/tests
)