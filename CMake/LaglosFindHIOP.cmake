
#[[

Exports target `HiOp::hiop_tpl`

Users may set the following variables:

]]

MESSAGE(STATUS "***** SEARCHING FOR HIOP *****")

find_package(hiop REQUIRED NAMES hiop HiOp HIOP HINTS LAGLOS_HIOP_DIR "${LAGLOS_HIOP_DIR}/install")

if (hiop_FOUND)
  MESSAGE(STATUS "HIOP FOUND")
endif()

