# -- NeoPepXMLParser ----------------------------------------------------------
cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0091 NEW)
project(neopepxmlparser)
# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${MSTOOLKIT_INCLUDE_PATH} 
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-pie")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

#
# The core NeoPepXMLParser files
#
#
set(
  NEOSRC
  CnpxAffectedChannel.cpp
  CnpxAlternativeProtein.cpp
  CnpxAminoAcidModification.cpp
  CnpxAminoAcidSubstitution.cpp
  CnpxAnalysisResult.cpp
  CnpxAnalysisSummary.cpp
  CnpxAnalysisTimestamp.cpp
  CnpxContributingChannel.cpp
  CnpxDatabaseRefreshTimestamp.cpp
  CnpxDecoyAnalysis.cpp
  CnpxDecoyAnalysisSummary.cpp
  CnpxDistributionPoint.cpp
  CnpxEnzymaticSearchConstraint.cpp
  CnpxErrorPoint.cpp
  CnpxFragmentMasses.cpp
  CnpxInputFile.cpp
  CnpxIntensity.cpp
  CnpxInteractSummary.cpp
  CnpxInterprophetResult.cpp
  CnpxInterprophetSummary.cpp
  CnpxIsotopicContributions.cpp
  CnpxLability.cpp
  CnpxLibraResult.cpp
  CnpxLibraSummary.cpp
  CnpxLinkedPeptide.cpp
  CnpxMSMSPipelineAnalysis.cpp
  CnpxMSMSRunSummary.cpp
  CnpxMixtureModel.cpp
  CnpxMixtureModelDistribution.cpp
  CnpxMixture_Model.cpp
  CnpxModAminoAcidMass.cpp
  CnpxModAminoAcidProbability.cpp
  CnpxModTerminalProbability.cpp
  CnpxModificationInfo.cpp
  CnpxNegModelDistribution.cpp
  CnpxPTMProphetResult.cpp
  CnpxPTMProphetSummary.cpp
  CnpxParameter.cpp
  CnpxPepXMLQuantResult.cpp
  CnpxPeptideProphetResult.cpp
  CnpxPeptideprophetSummary.cpp
  CnpxPosModelDistribution.cpp
  CnpxQuanticResult.cpp
  CnpxQuanticSummary.cpp
  CnpxROCDataPoint.cpp
  CnpxROCErrorData.cpp
  CnpxSampleEnzyme.cpp
  CnpxSearchDatabase.cpp
  CnpxSearchHit.cpp
  CnpxSearchResult.cpp
  CnpxSearchScore.cpp
  CnpxSearchScoreSummary.cpp
  CnpxSearchSummary.cpp
  CnpxSpecificity.cpp
  CnpxSpectrumQuery.cpp
  CnpxTerminalModification.cpp
  CnpxUIPSM.cpp
  CnpxUIPipeline.cpp
  CnpxUIRunSummary.cpp
  CnpxUISpectra.cpp
  CnpxXLink.cpp
  CnpxXLinkScore.cpp
  CnpxXpressLabelFreeResult.cpp
  CnpxXpressLabelFreeSummary.cpp
  NeoPepXMLParser.cpp
)

add_library(
  neopepxmlparser
  STATIC 
  ${NEOSRC}
)

if (UNIX)
  target_compile_options(
    neopepxmlparser
    PUBLIC
    -O3 -static -I. -I$(HEADER_PATH) -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC 
  )
elseif(WIN32 AND NOT CYGWIN)
  target_compile_definitions(
    neopepxmlparser
    PUBLIC
    _CRT_SECURE_NO_WARNINGS
  )
endif(UNIX)
  
################################################################################
								  
# INSTALL NeoPepXMLParser
install(TARGETS neopepxmlparser EXPORT neopepxmlparser DESTINATION lib)

