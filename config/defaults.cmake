set(MFEM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../mfem/build/install
    CACHE PATH "absolute path to the MFEM build or install prefix")
set(mfem_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../mfem/
    CACHE PATH "absolute path to where MFEMConfig.cmake is")

set(LAGLOS_HIOP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../hiop/)

option(LAGLOS_USE_HIOP "Utilize sparse linear algebra" OFF)
option(LAGLOS_USE_BOOST "Utilize BOOST" ON)
