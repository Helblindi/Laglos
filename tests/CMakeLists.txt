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

add_executable(TestInitiation test_laglos_creation.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(TestInitiation PRIVATE "${MFEM_LIBRARIES}" comp_lambda)
target_include_directories(TestInitiation PUBLIC "${MFEM_DIR}" "${MFEM_INCLUDE_DIRS}")
add_test(NAME TestInitiation COMMAND TestInitiation)

add_executable(TestCR test_CR.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(TestCR PRIVATE "${MFEM_LIBRARIES}" comp_lambda)
target_include_directories(TestCR PUBLIC "${MFEM_DIR}" "${MFEM_INCLUDE_DIRS}")
add_test(NAME TestCR COMMAND TestCR)

add_executable(TestVDLFGI test_VDLFGI.cpp "${MFEM_DIR}" "${PROJECT_SOURCE_DIR}/src/laglos_solver.cpp")
target_link_libraries(TestVDLFGI PRIVATE "${MFEM_LIBRARIES}" comp_lambda)
target_include_directories(TestVDLFGI PUBLIC "${MFEM_DIR}" "${MFEM_INCLUDE_DIRS}")
add_test(NAME TestVDLFGI COMMAND TestVDLFGI)