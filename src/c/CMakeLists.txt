cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

configure_file (
 "CruxConfig.h.in"
 "${CMAKE_CURRENT_BINARY_DIR}/CruxConfig.h"
 )

add_subdirectory(bullseye)
add_subdirectory(hardklor)
add_subdirectory(posterior-error)
add_subdirectory(qranker-barista)
add_subdirectory(xlink)
add_subdirectory(tide)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(.)
include_directories(bullseye)
include_directories(hardklor)
include_directories(qranker-barista)
include_directories(posterior-error)
include_directories(xlink)
include_directories(${CMAKE_SOURCE_DIR}/src/external)
include_directories(${CMAKE_BINARY_DIR}/src/external)
include_directories(${CMAKE_BINARY_DIR}/src/external/include)
include_directories(${CMAKE_BINARY_DIR}/src/external/include/CometSearch)
include_directories(${CMAKE_BINARY_DIR}/src/external/include/MSToolkit)
include_directories(${CMAKE_BINARY_DIR}/src/external/build/src/ProteoWizard/libraries/boost_1_54_0)
include_directories(${CMAKE_BINARY_DIR}/src/external/build/src/ProteoWizard/libraries/boost_aux)
include_directories(${CMAKE_BINARY_DIR}/src/external/build/src/percolator/src/)
include_directories(${CMAKE_BINARY_DIR}/src/external/build/src/percolator/src/fido)
include_directories(${CMAKE_BINARY_DIR}/src/external/build/src/percolator-build/src/)
include_directories(${CMAKE_BINARY_DIR}/src/external/build/src/percolator-build/src/fido)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/tide/protoobj)


link_directories(.)
link_directories(./bullseye)
link_directories(./hardklor)
link_directories(./posterior-error)
link_directories(./qranker-barista)
link_directories(./xlink)
link_directories(${CMAKE_BINARY_DIR}/src/external/lib)
link_directories(${CMAKE_BINARY_DIR}/src/external/build/src/percolator-build/src) 
link_directories(${CMAKE_BINARY_DIR}/src/external/build/src/percolator-build/src/fido)

set (
  crux_lib_files
  CruxConfig.h
  Alphabet.cpp
  analyze_psms.cpp
  carp.cpp
  CarpStreamBuf.cpp
  CometApplication.cpp
  ComputeQValues.cpp
  CreateIndex.cpp
  crux-file-utils.cpp
  crux-utils.cpp
  Database.cpp
  DatabaseProteinIterator.cpp
  DatabasePeptideIterator.cpp
  CruxApplication.cpp
  CruxApplicationList.cpp
  DelimitedFile.cpp
  DelimitedFileReader.cpp
  DelimitedFileWriter.cpp
  ExtractColumns.cpp
  ExtractRows.cpp
  FilteredSpectrumChargeIterator.cpp
  GenerateDecoys.cpp
  GeneratePeptidesIterator.cpp
  GetMs2Spectrum.cpp
  hash.cpp
  Index.cpp
  IndexFile.cpp
  IndexMap.cpp
  IndexPeptideIterator.cpp
  IndexFilteredPeptideIterator.cpp
  Ion.cpp
  IonConstraint.cpp
  IonFilteredIterator.cpp
  IonSeries.cpp
  linked_list.cpp
  LineFileReader.cpp
  MakePinApplication.cpp
  mass.cpp
  Match.cpp
  MatchColumns.cpp
  MatchFileReader.cpp
  MatchFileWriter.cpp
  MatchSearch.cpp
  MatchCollection.cpp
  MatchCollectionIterator.cpp
  MatchCollectionParser.cpp
  MatchIterator.cpp
  modifications.cpp
  ModifiedPeptidesIterator.cpp
  MSToolkitSpectrumCollection.cpp
  MzIdentMLReader.cpp
  MzIdentMLWriter.cpp
  OutputFiles.cpp
  parameter.cpp
  parse_arguments.cpp
  Peak.cpp
  Peptide.cpp
  PeptideConstraint.cpp
  peptide_modifications.cpp
  PeptideSrc.cpp
  PepXMLReader.cpp
  PepXMLWriter.cpp
  PercolatorAdapter.cpp
  PercolatorApplication.cpp
  GeneratePeptides.cpp
  PinXMLWriter.cpp
  PostProcessProtein.cpp
  PredictPeptideIons.cpp
  PrintProcessedSpectra.cpp
  PrintVersion.cpp
  Protein.cpp
  ProteinPeptideIterator.cpp
  ProteinIndex.cpp
  ProteinIndexIterator.cpp
  ProteinMatchCollection.cpp
  AbstractMatch.cpp
  ProteinMatch.cpp
  PeptideMatch.cpp
  PMCDelimitedFileWriter.cpp
  PMCPepXMLWriter.cpp
  PMCSQTWriter.cpp
  SpectrumMatch.cpp
  PWIZSpectrumCollection.cpp
  q-value.cpp
  ReadSpectrumRecordsApplication.cpp
  SequestSearch.cpp
  SortColumn.cpp
  Scorer.cpp
  sorter.cpp
  SpectralCounts.cpp
  SpectrumCollection.cpp
  SpectrumCollectionFactory.cpp
  Spectrum.cpp
  SpectrumRecordWriter.cpp
  SpectrumZState.cpp
  StatColumn.cpp
  SQTReader.cpp
  SQTWriter.cpp
  TideIndexApplication.cpp
  TideMatchSet.cpp
  TideSearchApplication.cpp
  utils.cpp
)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(
    crux_lib_files
    ${crux_lib_files}
    WinCrux.cpp
  )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

add_library(
  crux-support
  STATIC
  ${crux_lib_files}
)
add_dependencies(crux-support tide-support)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_property(
      DIRECTORY
      APPEND
      PROPERTY
      COMPILE_DEFINITIONS
      _CRT_SECURE_NO_WARNINGS
    )
else (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_property(
	  DIRECTORY
    APPEND
	  PROPERTY
	  COMPILE_DEFINITIONS
	  _LARGEFILE_SOURCE
	  _LARGEFILE64_SOURCE
	  _BSD_SOURCE
	  _POSIX_C_SOURCE=199506L
	  NDEBUG
	  _FILE_OFFSET_BITS=64
	  _GNU_SOURCE
	  GCC
	  CRUX
	  HAVE_EXPAT_CONFIG_H
	)
	if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	  set_property(
      DIRECTORY
      APPEND
      PROPERTY
      COMPILE_DEFINITIONS
      DARWIN
      _DARWIN_C_SOURCE
	  )
	  if (${CMAKE_SYSTEM_VERSION} VERSION_LESS "12")
      set_property(
        DIRECTORY
        APPEND
        PROPERTY
        COMPILE_DEFINITIONS
        PRE_MOUNTAIN_LION
      )
	  endif (${CMAKE_SYSTEM_VERSION} VERSION_LESS "12")
	  if (${CMAKE_SYSTEM_VERSION} VERSION_LESS "11")
      set_property(
        DIRECTORY
        APPEND
        PROPERTY
        COMPILE_DEFINITIONS
        PRE_LION
      )
    endif(${CMAKE_SYSTEM_VERSION} VERSION_LESS "11")
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

add_executable(crux crux-main.cpp)
set_property(
  TARGET
  crux
  PROPERTY
  COMPILE_DEFINITIONS
  MAIN
)

add_executable(xlink-predict-peptide-ions xlink/xhhc_predict_peptide_ions.cpp)
set_property(
  TARGET
  xlink-predict-peptide-ions
  PROPERTY
  COMPILE_DEFINITIONS
  MAIN
)

add_executable(xlink-score-spectrum xlink/xhhc_score_peptide_spectrum.cpp)
set_property(
  TARGET
  xlink-score-spectrum
  PROPERTY
  COMPILE_DEFINITIONS
  MAIN
)

add_executable(xlink-assign-ions xlink/xlink_assign_ions.cpp)
set_property(
  TARGET
  xlink-assign-ions
  PROPERTY
  COMPILE_DEFINITIONS
  MAIN
)

# The rt library is not needed for OS X
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(RT_LIBRARY "")
else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(RT_LIBRARY "rt")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  target_link_libraries(
    crux
    barista
    bullseye
    crux-support
    hardklor
    cometsearch
    mstoolkitlite
    perclibrary
    posterior_error
    xlink
    ws2_32
  )
  target_link_libraries(
    xlink-predict-peptide-ions
    xlink
    crux-support
    cometsearch
    mstoolkitlite
    ws2_32
  )
  target_link_libraries(
    xlink-score-spectrum
    xlink
    crux-support
    mstoolkitlite
    ws2_32
  )
  target_link_libraries(
    xlink-assign-ions
    xlink
    crux-support
    mstoolkitlite
    ws2_32
  )
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  target_link_libraries(
    crux
    barista
    bullseye
    crux-support
    tide-support
    protobuf
    gflags
    cometsearch
    dl
    hardklor
    m
    ${RT_LIBRARY}
    mstoolkitlite
    posterior_error
    perclibrary
    pthread
    xlink
    xerces-c
    boost_filesystem
    boost_system
    fido
    pwiz_data_identdata
    pwiz_data_identdata_version
    pwiz_data_misc
    pwiz_data_msdata
    pwiz_data_msdata_version
    pwiz_data_proteome
    pwiz_data_common
    pwiz_utility_chemistry
    pwiz_utility_minimxml
    pwiz_utility_misc
    pwiz_version
    SHA1
    boost_system-mt-s
    boost_iostreams-mt-s
    z-mt-s-1_2
    boost_filesystem-mt-s
    boost_regex-mt-s
    boost_thread-mt-s

  )
    
  target_link_libraries(
    xlink-predict-peptide-ions
    dl
    m
    ${RT_LIBRARY}
    pthread
    xlink
    crux-support
    mstoolkitlite
    pwiz_data_msdata
    pwiz_data_common
    pwiz_data_identdata
    pwiz_data_msdata_version
    pwiz_data_misc
    pwiz_utility_chemistry
    pwiz_utility_misc
    pwiz_utility_minimxml
    pwiz_data_common
    pwiz_data_identdata
    pwiz_data_identdata_version
    pwiz_data_proteome
    SHA1
    boost_system-mt-s
    boost_iostreams-mt-s
    z-mt-s-1_2
    boost_filesystem-mt-s
    boost_regex-mt-s
    boost_thread-mt-s
  )
  target_link_libraries(
    xlink-score-spectrum
    dl
    m
    ${RT_LIBRARY}
    pthread
    xlink
    crux-support
    mstoolkitlite
    pwiz_data_msdata
    pwiz_data_common
    pwiz_data_identdata
    pwiz_data_msdata_version
    pwiz_data_misc
    pwiz_utility_chemistry
    pwiz_utility_misc
    pwiz_utility_minimxml
    pwiz_data_common
    pwiz_data_identdata
    pwiz_data_identdata_version
    pwiz_data_proteome
    SHA1
    boost_system-mt-s
    boost_iostreams-mt-s
    z-mt-s-1_2
    boost_filesystem-mt-s
    boost_regex-mt-s
    boost_thread-mt-s
  )
  target_link_libraries(
    xlink-assign-ions
    dl
    m
    ${RT_LIBRARY}
    pthread
    xlink
    crux-support
    mstoolkitlite
    pwiz_data_msdata
    pwiz_data_common
    pwiz_data_identdata
    pwiz_data_msdata_version
    pwiz_data_misc
    pwiz_utility_chemistry
    pwiz_utility_misc
    pwiz_utility_minimxml
    pwiz_data_common
    pwiz_data_identdata
    pwiz_data_identdata_version
    pwiz_data_proteome
    SHA1
    boost_system-mt-s
    boost_iostreams-mt-s
    z-mt-s-1_2
    boost_filesystem-mt-s
    boost_regex-mt-s
    boost_thread-mt-s
  )
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install (
  TARGETS
  crux
  xlink-predict-peptide-ions
  xlink-score-spectrum
  xlink-assign-ions
  DESTINATION
  bin
)
