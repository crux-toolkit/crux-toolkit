project(xlink)

cmake_minimum_required(VERSION 2.6)
cmake_policy(VERSION 2.8.1)

include_directories(..)
include_directories(../external/MSToolkit/include)

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
  xlink
  STATIC
  XHHC_Peptide.cpp
  LinkedPeptide.cpp
  xhhc.cpp
  LinkedIonSeries.cpp
  xhhc_score_peptide_spectrum.cpp
  xhhc_scorer.cpp
  xhhc_score_peptide_spectrum.cpp
  xhhc_scorer.cpp
  xhhc_search.cpp
  xlink_search.cpp
  SearchForXLinks.cpp
  xlink_compute_qvalues.cpp
  XLinkSite.cpp
  XLinkBondMap.cpp
  XLinkablePeptide.cpp
)
