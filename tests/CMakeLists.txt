cmake_minimum_required(VERSION 3.10)

include_directories(${TEST_SRC_DIR}/src
                    ${MFEM_INCLUDE_DIRS}
                    ${MFEM_DIR})

# TESTING
add_executable(TestInstantiator testing.cpp "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(TestInstantiator PRIVATE "${MFEM_LIBRARIES}" comp_lambda)
target_include_directories(TestInstantiator PUBLIC "${MFEM_DIR}" "${MFEM_INCLUDE_DIRS}")
add_test(NAME TestInstantiator COMMAND TestInstantiator)

# add_executable(TestMM test_mm.cpp "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
# target_include_directories(TestMM PUBLIC "${MFEM_DIR}" "${MFEM_INCLUDE_DIRS}")
# add_test(NAME TestMM COMMAND TestMM)

add_executable(TestMM test_mm.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(TestMM PRIVATE "${MFEM_LIBRARIES}" comp_lambda)
target_include_directories(TestMM PUBLIC "${MFEM_DIR}" "${MFEM_INCLUDE_DIRS}")
add_test(NAME TestMM COMMAND TestMM)