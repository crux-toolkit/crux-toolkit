cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0091 NEW)
project(mstoolkit)
# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)

#############################################################################
# COMPILE MSToolkit
###############################################################################

set(EXPAT_PATH ./src/expat-2.2.9/lib)
set(MST_PATH ./src/MSToolkit)
set(MZIMLTOOLS_PATH ./src/mzIMLTools)
set(MZPARSER_PATH ./src/mzParser)
set(ZLIB_PATH ./src/zlib-1.2.11)
set(EXPAT_WITH_GETRANDOM OFF)
set(EXPAT_WITH_SYS_GETRANDOM OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

set(
  EXPAT 
  ${EXPAT_PATH}/xmlparse.c 
  ${EXPAT_PATH}/xmlrole.c 
  ${EXPAT_PATH}/xmltok.c
)
if (WIN32 AND NOT Cygwin)
  set(
    MSTOOLKIT 
    ${MST_PATH}/MSObject.cpp 
    ${MST_PATH}/MSReader.cpp 
    ${MST_PATH}/mzMLWriter.cpp 
    ${MST_PATH}/pepXMLWriter.cpp
    ${MST_PATH}/RAWReader.cpp 
    ${MST_PATH}/Spectrum.cpp 
  )
else (WIN32 AND NOT Cygwin)
  set(
    MSTOOLKIT 
    ${MST_PATH}/MSObject.cpp 
    ${MST_PATH}/MSReader.cpp 
    ${MST_PATH}/mzMLWriter.cpp 
    ${MST_PATH}/pepXMLWriter.cpp
    ${MST_PATH}/Spectrum.cpp 
  )
endif (WIN32 AND NOT Cygwin)

set(
  MZIMLTOOLS
  ${MZIMLTOOLS_PATH}/CAdditionalSearchParams.cpp 
  ${MZIMLTOOLS_PATH}/CAffiliation.cpp
  ${MZIMLTOOLS_PATH}/CAmbiguousResidue.cpp
  ${MZIMLTOOLS_PATH}/CAnalysisCollection.cpp 
  ${MZIMLTOOLS_PATH}/CAnalysisData.cpp
  ${MZIMLTOOLS_PATH}/CAnalysisParams.cpp 
  ${MZIMLTOOLS_PATH}/CAnalysisProtocolCollection.cpp
  ${MZIMLTOOLS_PATH}/CAnalysisSampleCollection.cpp 
  ${MZIMLTOOLS_PATH}/CAnalysisSoftware.cpp 
  ${MZIMLTOOLS_PATH}/CAnalysisSoftwareList.cpp 
  ${MZIMLTOOLS_PATH}/CAuditCollection.cpp 
  ${MZIMLTOOLS_PATH}/CBibliographicReference.cpp 
  ${MZIMLTOOLS_PATH}/CContactRole.cpp 
  ${MZIMLTOOLS_PATH}/CCvList.cpp 
  ${MZIMLTOOLS_PATH}/CDBSequence.cpp
  ${MZIMLTOOLS_PATH}/CDataCollection.cpp 
  ${MZIMLTOOLS_PATH}/CDatabaseFilters.cpp 
  ${MZIMLTOOLS_PATH}/CDatabaseName.cpp 
  ${MZIMLTOOLS_PATH}/CEnzyme.cpp 
  ${MZIMLTOOLS_PATH}/CEnzymeName.cpp 
  ${MZIMLTOOLS_PATH}/CEnzymes.cpp
  ${MZIMLTOOLS_PATH}/CFileFormat.cpp 
  ${MZIMLTOOLS_PATH}/CFilter.cpp 
  ${MZIMLTOOLS_PATH}/CFilterType.cpp 
  ${MZIMLTOOLS_PATH}/CFragmentArray.cpp 
  ${MZIMLTOOLS_PATH}/CFragmentation.cpp 
  ${MZIMLTOOLS_PATH}/CFragmentationTable.cpp 
  ${MZIMLTOOLS_PATH}/CFragmentTolerance.cpp 
  ${MZIMLTOOLS_PATH}/CExclude.cpp 
  ${MZIMLTOOLS_PATH}/CInclude.cpp 
  ${MZIMLTOOLS_PATH}/CInputs.cpp 
  ${MZIMLTOOLS_PATH}/CIonType.cpp
  ${MZIMLTOOLS_PATH}/CMassTable.cpp 
  ${MZIMLTOOLS_PATH}/CMeasure.cpp
  ${MZIMLTOOLS_PATH}/CModification.cpp 
  ${MZIMLTOOLS_PATH}/CModificationParams.cpp 
  ${MZIMLTOOLS_PATH}/CMzIdentML.cpp 
  ${MZIMLTOOLS_PATH}/COrganization.cpp 
  ${MZIMLTOOLS_PATH}/CParent.cpp 
  ${MZIMLTOOLS_PATH}/CParentTolerance.cpp 
  ${MZIMLTOOLS_PATH}/CPeptide.cpp 
  ${MZIMLTOOLS_PATH}/CPeptideEvidence.cpp 
  ${MZIMLTOOLS_PATH}/CPeptideHypothesis.cpp
  ${MZIMLTOOLS_PATH}/CPerson.cpp 
  ${MZIMLTOOLS_PATH}/CProteinAmbiguityGroup.cpp 
  ${MZIMLTOOLS_PATH}/CProteinDetection.cpp 
  ${MZIMLTOOLS_PATH}/CProteinDetectionHypothesis.cpp  
  ${MZIMLTOOLS_PATH}/CProteinDetectionList.cpp 
  ${MZIMLTOOLS_PATH}/CProteinDetectionProtocol.cpp
  ${MZIMLTOOLS_PATH}/CProvider.cpp 
  ${MZIMLTOOLS_PATH}/CPSM.cpp 
  ${MZIMLTOOLS_PATH}/CResidue.cpp 
  ${MZIMLTOOLS_PATH}/CRole.cpp
  ${MZIMLTOOLS_PATH}/CSample.cpp
	${MZIMLTOOLS_PATH}/CSearchDatabase.cpp 
  ${MZIMLTOOLS_PATH}/CSearchModification.cpp 
  ${MZIMLTOOLS_PATH}/CSearchType.cpp 
  ${MZIMLTOOLS_PATH}/CSequenceCollection.cpp 
  ${MZIMLTOOLS_PATH}/CSiteRegexp.cpp 
  ${MZIMLTOOLS_PATH}/CSoftwareName.cpp
  ${MZIMLTOOLS_PATH}/CSourceFile.cpp 
  ${MZIMLTOOLS_PATH}/CSpecificityRules.cpp 
  ${MZIMLTOOLS_PATH}/CSpectraData.cpp 
  ${MZIMLTOOLS_PATH}/CSpectrumIDFormat.cpp
  ${MZIMLTOOLS_PATH}/CSpectrumIdentification.cpp 
  ${MZIMLTOOLS_PATH}/CSpectrumIdentificationItem.cpp
  ${MZIMLTOOLS_PATH}/CSpectrumIdentificationList.cpp 
  ${MZIMLTOOLS_PATH}/CSpectrumIdentificationProtocol.cpp 
  ${MZIMLTOOLS_PATH}/CSpectrumIdentificationResult.cpp
  ${MZIMLTOOLS_PATH}/CSubSample.cpp 
  ${MZIMLTOOLS_PATH}/CThreshold.cpp
)
set(
  MZPARSERLITE
  ${MZPARSER_PATH}/MSNumpress.cpp 
  ${MZPARSER_PATH}/mzp_base64.cpp 
  ${MZPARSER_PATH}/BasicSpectrum.cpp 
  ${MZPARSER_PATH}/mzParser.cpp 
  ${MZPARSER_PATH}/RAMPface.cpp 
  ${MZPARSER_PATH}/saxhandler.cpp 
  ${MZPARSER_PATH}/saxmzmlhandler.cpp
  ${MZPARSER_PATH}/saxmzxmlhandler.cpp 
  ${MZPARSER_PATH}/Czran.cpp 
  ${MZPARSER_PATH}/mz5handler.cpp 
  ${MZPARSER_PATH}/mzpMz5Config.cpp 
  ${MZPARSER_PATH}/mzpMz5Structs.cpp 
  ${MZPARSER_PATH}/BasicChromatogram.cpp 
  ${MZPARSER_PATH}/PWIZface.cpp
)
set(
  ZLIB 
  ${ZLIB_PATH}/adler32.c 
  ${ZLIB_PATH}/compress.c 
  ${ZLIB_PATH}/crc32.c
  ${ZLIB_PATH}/deflate.c 
  ${ZLIB_PATH}/inffast.c 
  ${ZLIB_PATH}/inflate.c
  ${ZLIB_PATH}/infback.c 
  ${ZLIB_PATH}/inftrees.c 
  ${ZLIB_PATH}/trees.c
  ${ZLIB_PATH}/uncompr.c 
  ${ZLIB_PATH}/zutil.c
)
add_library(
  mstoolkitlite 
  STATIC 
  ${ZLIB} 
  ${EXPAT} 
  ${MZIMLTOOLS} 
  ${MZPARSERLITE}
  ${MSTOOLKIT}
)
								  
if (UNIX)
  target_compile_options(
    mstoolkitlite
    PUBLIC
    -O3 -static -I. -I$(HEADER_PATH) -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC
    -D_NOSQLITE -DXML_POOR_ENTROPY
  )
elseif(WIN32 AND NOT CYGWIN)
  target_compile_definitions(
    mstoolkitlite
    PUBLIC
    _CRT_SECURE_NO_WARNINGS
  )
endif(UNIX)
  

################################################################################
								  
# INSTALL MSTOOLKIT
install(TARGETS mstoolkitlite EXPORT mstoolkitlite DESTINATION lib)
install(DIRECTORY include/ DESTINATION include/MSToolkit FILES_MATCHING PATTERN "*") 


