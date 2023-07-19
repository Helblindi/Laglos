cmake_minimum_required(VERSION 3.10)

include_directories(${TEST_SRC_DIR}/src
                    ${MFEM_INCLUDE_DIRS}
                    ${MFEM_DIR})

# TESTING
add_executable(TestInstantiator EXCLUDE_FROM_ALL testing.cpp "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(TestInstantiator PRIVATE "${MFEM_LIBRARIES}" comp_lambda)
target_include_directories(TestInstantiator PUBLIC "${MFEM_DIR}" "${MFEM_INCLUDE_DIRS}")
add_test(NAME TestInstantiator COMMAND TestInstantiator)

add_executable(TestMM EXCLUDE_FROM_ALL test_mm.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(TestMM PRIVATE "${MFEM_LIBRARIES}" comp_lambda)
target_include_directories(TestMM PUBLIC "${MFEM_DIR}" "${MFEM_INCLUDE_DIRS}")
add_test(NAME TestMM COMMAND TestMM)

# add_executable(TestInitiation EXCLUDE_FROM_ALL test_laglos_creation.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
# target_link_libraries(TestInitiation PRIVATE "${MFEM_LIBRARIES}" comp_lambda)
# target_include_directories(TestInitiation PUBLIC "${MFEM_DIR}" "${MFEM_INCLUDE_DIRS}")
# add_test(NAME TestInitiation COMMAND TestInitiation)

add_executable(TestCR EXCLUDE_FROM_ALL test_CR.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(TestCR PRIVATE "${MFEM_LIBRARIES}" comp_lambda)
target_include_directories(TestCR PUBLIC "${MFEM_DIR}" "${MFEM_INCLUDE_DIRS}")
add_test(NAME TestCR COMMAND TestCR)

add_executable(TestRTvel EXCLUDE_FROM_ALL test_RT_vel.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(TestRTvel PRIVATE "${MFEM_LIBRARIES}" comp_lambda)
target_include_directories(TestRTvel PUBLIC "${MFEM_DIR}" "${MFEM_INCLUDE_DIRS}")
add_test(NAME TestRTvel COMMAND TestRTvel)

add_executable(Test1D EXCLUDE_FROM_ALL test_1D.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(Test1D PRIVATE "${MFEM_LIBRARIES}" comp_lambda)
target_include_directories(Test1D PUBLIC "${MFEM_DIR}" "${MFEM_INCLUDE_DIRS}")
add_test(NAME Test1D COMMAND Test1D)

add_executable(Test2D EXCLUDE_FROM_ALL test_2D.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(Test2D PRIVATE "${MFEM_LIBRARIES}" comp_lambda)
target_include_directories(Test2D PUBLIC "${MFEM_DIR}" "${MFEM_INCLUDE_DIRS}")
add_test(NAME Test2D COMMAND Test2D)

# Create a new target to build and run tests
# This ensure make all doesn't make every test at each compilation stage
add_custom_target(build_tests)
add_custom_target(run_tests COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(run_tests build_tests)
# add_dependencies(build_tests TestInstantiator TestMM TestInitiation TestCR TestRTvel Test1D Test2D)
add_dependencies(build_tests TestInstantiator TestMM TestCR TestRTvel Test1D Test2D)