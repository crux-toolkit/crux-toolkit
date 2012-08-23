project(percolator)

cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

include_directories(..)
include_directories(../ArrayLibrary)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_property(
    DIRECTORY
    PROPERTY
    COMPILE_DEFINITIONS
      _CRT_SECURE_NO_WARNINGS
    )
else (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_property(
    DIRECTORY
    PROPERTY
    COMPILE_DEFINITIONS
    _FILE_OFFSET_BITS=64
    _GNU_SOURCE
    _LARGEFILE_SOURCE
    GCC
  )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

add_library(
  percolator
  STATIC
  PercolatorApplication.cpp
)
