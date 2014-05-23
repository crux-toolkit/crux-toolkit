project(xdr)
cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

if (WIN32 AND NOT CYGWIN)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  set (
    xdr_lib_files
    xdr.c
    xdr_array.c
    xdr_float.c
    xdr_mem.c
    xdr_private.c
    xdr_rec.c
    xdr_reference.c
    xdr_sizeof.c
    xdr_stdio.c
  )

add_library(
  xdr
  STATIC
  ${xdr_lib_files}
)

set_property(
  DIRECTORY
  APPEND
  PROPERTY
  COMPILE_DEFINITIONS
  _CRT_SECURE_NO_WARNINGS
)

install(TARGETS xdr DESTINATION ${CMAKE_PREFIX_PATH}/lib) # Important to use relative path here (used by CPack)!
endif (WIN32 AND NOT CYGWIN)
