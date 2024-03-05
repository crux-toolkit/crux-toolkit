cmake_minimum_required(VERSION 3.15)
cmake_policy(VERSION 3.15)

project(hardklor)

if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-pie")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
endif (UNIX)

include_directories(${CMAKE_SOURCE_DIR}/src)
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
    NOMINMAX
  )
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
else (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_property(
    DIRECTORY
    PROPERTY
    COMPILE_DEFINITIONS
    _FILE_OFFSET_BITS=64
    _GNU_SOURCE
    _LARGEFILE_SOURCE
    CRUX
    GCC
    HAVE_EXPAT_CONFIG_H
  )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

add_library(
  hardklor
  STATIC
  CAveragine.cpp
  CHardklor.cpp
  CHardklor2.cpp
  CHardklorParser.cpp
  CHardklorProtein.cpp
  CHardklorSetting.cpp
  CHardklorVariant.cpp
  CMercury8.cpp
  CModelLibrary.cpp
  CNoiseReduction.cpp
  CPeriodicTable.cpp
  CSpecAnalyze.cpp
  CSplitSpectrum.cpp
  FFT.cpp
  FFT-HK.cpp
  S2N.cpp
  Smooth.cpp
  SpecAnalyzeSupport.cpp
  CruxHardklorApplication.cpp
)
