
#[[

Exports target `COINHSL`

Users may set the following variables:

- HIOP_COINHSL_DIR

]]

MESSAGE(STATUS "***** SEARCHING FOR COINHSL *****")

find_library(COINHSL_LIBRARY
  NAMES
  coinhsl
  PATHS
  ${COINHSL_DIR} $ENV{COINHSL_DIR} ${LAGLOS_COINHSL_DIR}
  PATH_SUFFIXES
  lib64 lib a so lib/x86_64-linux-gnu)

if(COINHSL_LIBRARY)
  get_filename_component(COINHSL_LIBRARY_DIR ${COINHSL_LIBRARY} DIRECTORY)
endif()

find_path(COINHSL_INCLUDE_DIR
  NAMES
  CoinHslConfig.h
  PATHS
  ${COINHSL_DIR} $ENV{COINHSL_DIR} ${LAGLOS_COINHSL_DIR} ${COINHSL_LIBRARY_DIR}/..
  PATH_SUFFIXES
  include/coin/ThirdParty
  include/coinor
  include/coin-or
  include/coin-or/hsl
  include)

if(COINHSL_LIBRARY)
  message(STATUS "Found coinhsl include: ${COINHSL_INCLUDE_DIR}")
  message(STATUS "Found coinhsl library: ${COINHSL_LIBRARY}")
  add_library(COINHSL INTERFACE)
  target_link_libraries(COINHSL INTERFACE ${COINHSL_LIBRARY} ${CMAKE_DL_LIBS})
  target_include_directories(COINHSL INTERFACE ${COINHSL_INCLUDE_DIR})
  message(STATUS "Found COINHSL library: ${COINHSL_LIBRARY}")
  install(TARGETS COINHSL EXPORT hiop-targets)
else()
  message(STATUS "COINHSL was not found.")
endif()

set(COINHSL_INCLUDE_DIR CACHE PATH "Path to coinhsl.h")
set(COINHSL_LIBRARY CACHE PATH "Path to coinhsl library")
