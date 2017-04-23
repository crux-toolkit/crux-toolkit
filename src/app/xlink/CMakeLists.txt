project(xlink)

cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_1_56_0)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_aux)
include_directories(${CMAKE_BINARY_DIR}/ext/include)

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
  xlink_assign_ions.cpp
  xlink_search.cpp
  SearchForXLinks.cpp
  xlink_compute_qvalues.cpp
  LinearPeptide.cpp
  MonoLinkPeptide.cpp
  SelfLoopPeptide.cpp
  XLinkDatabase.cpp
  XLinkIonSeriesCache.cpp
  XLink.cpp
  XLinkScorer.cpp
  XLinkSite.cpp
  XLinkBondMap.cpp
  XLinkablePeptide.cpp
  XLinkablePeptideIterator.cpp
  XLinkablePeptideIteratorTopN.cpp
  XLinkMatch.cpp
  XLinkMatchCollection.cpp
  XLinkPeptide.cpp
  Weibull.cpp
)
