list(APPEND HDRS 
   geometry.hpp
   compile_time_vals.h 
   lagrange-multiplier.hpp
   riemann1D.hpp
   initial_vals.hpp
   laglos_assembly.hpp
   laglos_solver.hpp
   lagrange_multiplier.hpp
   lambda_max_lagrange.h
   sedov_exact.hpp
   limiter.h
   eos.h
   elastic.hpp
   # "${PROJECT_SOURCE_DIR}/../mfem/build/mfem.hpp"
   )

convert_filenames_to_full_paths(HDRS)

set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)