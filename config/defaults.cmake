set(MFEM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../mfem/)
set(mfem_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../mfem/)

set(LAGLOS_COINHSL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../coinhsl/)
set(LAGLOS_HIOP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../hiop/)

option(LAGLOS_USE_HIOP "Utilize sparse linear algebra" ON)
option(LAGLOS_USE_COINHSL "Use COINHSL for sparse linear algebra" ON)
