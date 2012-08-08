project(barista)

cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

include_directories(..)

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
