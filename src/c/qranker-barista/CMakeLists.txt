project(barista)

cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

include_directories(${CMAKE_SOURCE_DIR}/src/c)
include_directories(${CMAKE_BINARY_DIR}/src/external/build/src/ProteoWizard/libraries/boost_1_54_0)
include_directories(${CMAKE_BINARY_DIR}/src/external/build/src/ProteoWizard/libraries/boost_aux)
include_directories(${CMAKE_BINARY_DIR}/src/external/include)

if (WIN32 AND NOT CYGWIN)
  add_definitions(-DBOOST_ALL_NO_LIB)
endif (WIN32 AND NOT CYGWIN)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_property(
    DIRECTORY
    PROPERTY
    COMPILE_DEFINITIONS
      _CRT_SECURE_NO_WARNINGS
      _SCL_SECURE_NO_WARNINGS
      CRUX
      SEARCH_ENABLED
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
  barista
  STATIC
  Barista.cpp
  BipartiteGraph.cpp
  CruxParser.cpp
  DataSet.cpp
  NeuralNet.cpp
  PepRanker.cpp
  PepScores.cpp
  ProtScores.cpp
  PSMScores.cpp
  QRanker.cpp
  SpecFeatures.cpp
  SQTParser.cpp
  TabDelimParser.cpp
)
