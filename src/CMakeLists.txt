list(APPEND SRCS 
   laglos.cpp
   laglos_solver.cpp
   sedov_exact.cpp
   # lagrange_multiplier.cpp
   geometry.cpp)

convert_filenames_to_full_paths(SRCS)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)