list(APPEND HEADERS
   "${PROJECT_SOURCE_DIR}/include/compile_time_vals.h"
   "${PROJECT_SOURCE_DIR}/include/initial_vals.hpp"
   "${PROJECT_SOURCE_DIR}/include/laglos_solver.hpp"
   "${PROJECT_SOURCE_DIR}/include/riemann1D.hpp"
   "${PROJECT_SOURCE_DIR}/include/lambda_max_lagrange.h"
   "${PROJECT_SOURCE_DIR}/include/sedov_exact.hpp"
   "${PROJECT_SOURCE_DIR}/include/geometry.hpp"
   "${PROJECT_SOURCE_DIR}/include/lagrange_multiplier.hpp"
   "${PROJECT_SOURCE_DIR}/../mfem-git/build/mfem.hpp")

list(APPEND SOURCES 
   laglos.cpp
   laglos_solver.cpp
   # lagrange_multiplier.cpp
   geometry.cpp)


MESSAGE(STATUS "_laglos_libraries: ${_laglos_libraries}")
MESSAGE(STATUS "_laglos_include_dirs: ${_laglos_include_dirs}")
MESSAGE(STATUS "HEADERS: ${HEADERS}")
MESSAGE(STATUS "SOURCES: ${SOURCES}")

add_library(laglos ${SOURCES} ${HEADERS})

target_compile_definitions(laglos PUBLIC "${_laglos_compile_defs}")
target_compile_options(laglos PUBLIC "${_laglos_compile_opts}")
target_include_directories(laglos SYSTEM PUBLIC "${_laglos_include_dirs}")
target_link_libraries(laglos PUBLIC "${_laglos_libraries}")
