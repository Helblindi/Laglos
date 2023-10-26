# add the executable
add_executable(Laglos laglos.cpp laglos_solver.cpp)
set(EXECUTABLE_OUTPUT_PATH "../")

# Add fortran
message(STATUS "src message")
message(STATUS "LAGLOS_SOURCE_DIR=${PROJECT_SOURCE_DIR}")
message(STATUS "MFEM_LIBRARIES=${MFEM_DIR}")

# Link MFEM
target_link_libraries(Laglos "${MFEM_LIBRARIES}")

# Include libraries needed by lboost_filesystem
# https://stackoverflow.com/questions/9235679/create-a-directory-if-it-doesnt-exist

target_link_libraries(Laglos Boost::filesystem)

# set(LBOOST_COMPILE_FLAGS "-lboost_system")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LBOOST_COMPILE_FLAGS}")

# Specify Include directories for target
target_include_directories(Laglos PUBLIC ${MFEM_DIR} "${MFEM_INCLUDE_DIRS}")
target_link_libraries(Laglos comp_lambda_lagrangian comp_lambda_eulerian vdw_fortran)