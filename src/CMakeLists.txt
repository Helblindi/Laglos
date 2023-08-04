# add the executable
add_executable(Laglos laglos.cpp laglos_solver.cpp)
set(EXECUTABLE_OUTPUT_PATH "../")

# Add fortran
set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_SOURCE_DIR}/Fortran")
message(STATUS "src message")
message(STATUS "LAGLOS_SOURCE_DIR=${PROJECT_SOURCE_DIR}")
message(STATUS "MFEM_LIBRARIES=${MFEM_DIR}")

# Link MFEM and Fortran libraries
target_link_libraries(Laglos "${MFEM_LIBRARIES}" comp_lambda_lagrangian comp_lambda_eulerian)

# Specify Include directories for target
target_include_directories(Laglos PUBLIC ${MFEM_DIR} "${MFEM_INCLUDE_DIRS}")