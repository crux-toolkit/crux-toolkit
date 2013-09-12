project(hardklor)

cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

include_directories(..)
include_directories(../../external/include)
include_directories(../../external/build/src/ProteoWizard/libraries/boost_1_43_0)
include_directories(../../external/build/src/ProteoWizard/libraries/boost_aux)
include_directories(../../external/build/src/comet/src/comet/MSToolkit/include)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_property(
      DIRECTORY
      PROPERTY
      COMPILE_DEFINITIONS
      _CRT_SECURE_NO_WARNINGS
      CRUX
    )
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
  CHardklorFileReader.cpp
  CHardklorParser.cpp
  CHardklorProtein.cpp
  CHardklorSetting.cpp
  CHardklorVariant.cpp
  CMercury8.cpp
  CModelLibrary.cpp
  CNoiseReduction.cpp
  CPeriodicTable.cpp
  CSettings.cpp
  CSpecAnalyze.cpp
  CSplitSpectrum.cpp
  FFT.cpp
  FFT-HK.cpp
  HardklorApp.cpp
  S2N.cpp
  Smooth.cpp
  SpecAnalyzeSupport.cpp
  CruxHardklorApplication.cpp
)
