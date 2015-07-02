project(bullseye)

cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src/app/hardklor)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_1_56_0)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_aux)
include_directories(${CMAKE_BINARY_DIR}/ext/include)
include_directories(${CMAKE_BINARY_DIR}/ext/include/MSToolkit)
if (WIN32 AND NOT Cygwin)
  # Needed to put DLL containing type libraries
  # in include path for Windows
  include_directories(${CMAKE_BINARY_DIR}/ext/lib)
endif (WIN32 AND NOT Cygwin)

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
