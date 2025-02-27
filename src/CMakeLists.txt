list(APPEND SRCS 
   laglos.cpp
   laglos_solver.cpp
   laglos_assembly.cpp
   # lagrange_multiplier.cpp
   geometry.cpp)

convert_filenames_to_full_paths(SRCS)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)