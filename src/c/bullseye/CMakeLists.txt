project(bullseye)

cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

include_directories(..)
include_directories(../hardklor)
include_directories(../../external/MSToolkit/include)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_property(
      DIRECTORY
      PROPERTY
      COMPILE_DEFINITIONS
      _CRT_SECURE_NO_WARNINGS
      CRUX
      SEARCH_ENABLED
    )
else (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_property(
    DIRECTORY
    PROPERTY
    COMPILE_DEFINITIONS
    _LARGEFILE_SOURCE
    _FILE_OFFSET_BITS=64
    _GNU_SOURCE
    CRUX
    GCC
    HAVE_EXPAT_CONFIG_H
    SEARCH_ENABLED
  )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

add_library(
  bullseye
  STATIC
  bullseye.cpp
  CKronik2.cpp
  CruxBullseyeApplication.cpp
)
