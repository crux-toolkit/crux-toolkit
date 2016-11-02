cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

add_subdirectory(app/bullseye)
add_subdirectory(app/hardklor)
add_subdirectory(app/qranker-barista)
add_subdirectory(app/xlink)
add_subdirectory(app/tide)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(.)
include_directories(app/bullseye)
include_directories(app/hardklor)
include_directories(app/qranker-barista)
include_directories(app/xlink)
include_directories(io)
include_directories(model)
include_directories(${CMAKE_SOURCE_DIR}/ext)
include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_BINARY_DIR}/ext)
include_directories(${CMAKE_BINARY_DIR}/ext/include)
include_directories(${CMAKE_BINARY_DIR}/ext/include/CometSearch)
include_directories(${CMAKE_BINARY_DIR}/ext/include/MSToolkit)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_1_56_0)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_aux)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator/src/)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator/src/fido)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator/src/fisher)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src/)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src/fido)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/app/tide/protoobj)
if (WIN32 AND NOT Cygwin)
  # Needed to put DLL containing type libraries
  # in include path for Windows
  include_directories(${CMAKE_BINARY_DIR}/ext/lib)
  # Needed to include gflags/gflags.h on Windows
  include_directories(${CMAKE_BINARY_DIR}/ext/build/src/google-gflags/src/windows)
  # Needed to include google/protobuf on Windows
  include_directories(${CMAKE_BINARY_DIR}/ext/build/src/protobuf-2.5.0/src)
  SET(
    CMAKE_EXE_LINKER_FLAGS 
    "${CMAKE_EXE_LINKER_FLAGS} /VERBOSE:LIB /NODEFAULTLIB:libboost_date_time-vc120-mt-1_56 /NODEFAULTLIB:libboost_system-vc120-mt-1_56 /NODEFAULTLIB:libboost_thread-vc120-mt-1_56"
  )
  SET(
  CMAKE_EXE_LINKER_FLAGS_DEBUG 
  "${CMAKE_EXE_LINKER_FLAGS}
  /VERBOSE:LIB /NODEFAULTLIB:libboost_chrono-vc120-mt-gd-1_56 /NODEFAULTLIB:libboost_date_time-vc120-mt-gd-1_56 /NODEFAULTLIB:libboost_system-vc120-mt-gd-1_56 /NODEFAULTLIB:libboost_thread-vc120-mt-gd-1_56"
)
  #  Need Comet base directory for Windows
  include_directories(${CMAKE_BINARY_DIR}/ext/build/src/comet/)

  # On Windows we should build vendor libraries by default
  if (NOT DEFINED INCLUDE_VENDOR_LIBRARIES)
    set(INCLUDE_VENDOR_LIBRARIES TRUE)
  endif (NOT DEFINED INCLUDE_VENDOR_LIBRARIES)

endif (WIN32 AND NOT Cygwin)

link_directories(.)
link_directories(./app/bullseye)
link_directories(./app/hardklor)
link_directories(./app/qranker-barista)
link_directories(./app/xlink)
link_directories(${CMAKE_BINARY_DIR}/ext/lib)
link_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src) 
link_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src/fido)
if (WIN32 AND NOT Cygwin)
 # Needed for windows build of vendor libraries
  if (CMAKE_GENERATOR MATCHES "^.*Win64$")
    link_directories(
      ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/pwiz_aux/msrc/utility/vendor_api/Waters/vc12_x64
      ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/pwiz_aux/msrc/utility/vendor_api/Bruker/x64
    )
  else (CMAKE_GENERATOR MATCHES "^.*Win64$")
    link_directories(
      ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/pwiz_aux/msrc/utility/vendor_api/Waters/vc12_x86
      ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/pwiz_aux/msrc/utility/vendor_api/Bruker/x86
    )
  endif (CMAKE_GENERATOR MATCHES "^.*Win64$")
endif (WIN32 AND NOT Cygwin)

set (
  crux_lib_files
  app/SubtractIndexApplication.cpp
  app/CascadeSearchApplication.cpp
  app/AssignConfidenceApplication.cpp
  util/Alphabet.cpp
  util/AminoAcidUtil.cpp
  util/ArgParser.cpp
  util/CacheableMass.cpp
  util/GlobalParams.cpp
  io/carp.cpp
  util/CarpStreamBuf.cpp
  app/CometApplication.cpp
  app/ComputeQValues.cpp
  app/CreateDocs.cpp
  util/crux-utils.cpp
  model/Database.cpp
  model/DatabaseProteinIterator.cpp
  model/DatabasePeptideIterator.cpp
  app/CruxApplication.cpp
  app/CruxApplicationList.cpp
  io/DelimitedFile.cpp
  io/DelimitedFileReader.cpp
  io/DelimitedFileWriter.cpp
  app/ExtractColumns.cpp
  app/ExtractRows.cpp
  util/FileUtils.cpp
  model/FilteredSpectrumChargeIterator.cpp
  model/GeneratePeptidesIterator.cpp
  app/GetMs2Spectrum.cpp
  io/HTMLWriter.cpp
  model/Ion.cpp
  model/IonConstraint.cpp
  model/IonFilteredIterator.cpp
  model/IonSeries.cpp
  util/linked_list.cpp
  io/LineFileReader.cpp
  util/mass.cpp
  app/MakePinApplication.cpp
  model/Match.cpp
  io/MatchColumns.cpp
  io/MatchFileReader.cpp
  io/MatchFileWriter.cpp
  model/MatchCollection.cpp
  model/MatchCollectionIterator.cpp
  io/MatchCollectionParser.cpp
  model/MatchIterator.cpp
  util/MathUtil.cpp
  model/Modification.cpp
  util/modifications.cpp
  model/ModifiedPeptidesIterator.cpp
  io/MSToolkitSpectrumCollection.cpp
  io/MzIdentMLReader.cpp
  io/MzIdentMLWriter.cpp
  io/OutputFiles.cpp
  parameter.cpp
  app/ParamMedicApplication.cpp
  util/Params.cpp
  model/Peak.cpp
  model/Peptide.cpp
  model/PeptideConstraint.cpp
  util/peptide_modifications.cpp
  model/PeptideSrc.cpp
  io/PepXMLReader.cpp
  io/PepXMLWriter.cpp
  app/PercolatorAdapter.cpp
  app/PercolatorApplication.cpp
  io/PinWriter.cpp
  app/Pipeline.cpp
  app/GeneratePeptides.cpp
  model/PostProcessProtein.cpp
  app/PredictPeptideIons.cpp
  app/PrintProcessedSpectra.cpp
  app/PrintVersion.cpp
  model/Protein.cpp
  model/ProteinPeptideIterator.cpp
  model/ProteinIndex.cpp
  model/ProteinIndexIterator.cpp
  model/ProteinMatchCollection.cpp
  app/PSMConvertApplication.cpp
  io/PSMReader.cpp
  io/PSMWriter.cpp
  model/AbstractMatch.cpp
  model/ProteinMatch.cpp
  model/PeptideMatch.cpp
  io/PMCDelimitedFileWriter.cpp
  io/PMCPepXMLWriter.cpp
  io/PMCSQTWriter.cpp
  model/SpectrumMatch.cpp
  io/PWIZSpectrumCollection.cpp
  app/ReadSpectrumRecordsApplication.cpp
  app/ReadTideIndex.cpp
  app/SortColumn.cpp
  model/Scorer.cpp
  app/SpectralCounts.cpp
  io/SpectrumCollection.cpp
  io/SpectrumCollectionFactory.cpp
  model/Spectrum.cpp
  io/SpectrumRecordWriter.cpp
  model/SpectrumZState.cpp
  app/StatColumn.cpp
  util/StringUtils.cpp
  io/SQTReader.cpp
  io/SQTWriter.cpp
  app/TideIndexApplication.cpp
  app/TideMatchSet.cpp
  app/TideSearchApplication.cpp
  util/utils.cpp
)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(
    crux_lib_files
    ${crux_lib_files}
    util/WinCrux.cpp
  )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

add_library(
  crux-support
  STATIC
  ${crux_lib_files}
)
add_dependencies(crux-support tide-support)

find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG")

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_property(
      DIRECTORY
      APPEND
      PROPERTY
      COMPILE_DEFINITIONS
      _CRT_SECURE_NO_WARNINGS
     NOMINMAX
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
if (WIN32 AND NOT CYGWIN)
  set_property(
    TARGET crux 
    PROPERTY 
      COMPILE_DEFINITIONS 
      MAIN
      GFLAGS_DLL_DECL= 
      GFLAGS_DLL_DECLARE_FLAG= 
      GFLAGS_DLL_DEFINE_FLAG=
  )
else (WIN32 AND NOT CYGWIN)
  set_property(
    TARGET
    crux
    PROPERTY
    COMPILE_DEFINITIONS
    MAIN
  )
endif (WIN32 AND NOT CYGWIN)

if (WIN32 AND NOT CYGWIN)
  # These definititions are needed to set the
  # the proper calling conventions (static lib)
  # for Boost and expat. Expat is used in the
  # MSToolkit library.
  add_definitions(-DBOOST_ALL_NO_LIB)
  add_definitions(-DXML_STATIC)
  if (INCLUDE_VENDOR_LIBRARIES)
    add_definitions(-DINCLUDE_VENDOR_LIBRARIES)
  endif (INCLUDE_VENDOR_LIBRARIES)
endif (WIN32 AND NOT CYGWIN)

# The rt library is not needed for OS X
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(RT_LIBRARY "")
else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(RT_LIBRARY "rt")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if (WIN32 AND NOT CYGWIN)
  if (INCLUDE_VENDOR_LIBRARIES)
    set (
      PWIZ_LIBRARIES
      libpwiz_data_identdata
      libpwiz_data_identdata_version
      libpwiz_data_misc
      libpwiz_data_msdata
      libpwiz_data_msdata_version
      libpwiz_data_proteome
      libpwiz_data_common
      libpwiz_data_vendor_readers
      libpwiz_reader_abi
      libpwiz_reader_abi_t2d
      libpwiz_reader_agilent
      libpwiz_reader_shimadzu
      libpwiz_reader_thermo
      libpwiz_reader_waters
      libpwiz_vendor_api_abi
      libpwiz_vendor_api_agilent
      libpwiz_vendor_api_shimadzu
      libpwiz_vendor_api_thermo
      libpwiz_utility_chemistry
      libpwiz_utility_minimxml
      libpwiz_utility_misc
      libpwiz_version
      libSHA1
      libzlib
      MassLynxRaw
      cdt
    )
  else (INCLUDE_VENDOR_LIBRARIES)
    set (
      PWIZ_LIBRARIES
      libpwiz_data_identdata
      libpwiz_data_identdata_version
      libpwiz_data_misc
      libpwiz_data_msdata
      libpwiz_data_msdata_version
      libpwiz_data_proteome
      libpwiz_data_common
      libpwiz_utility_chemistry
      libpwiz_utility_minimxml
      libpwiz_utility_misc
      libpwiz_version
      libSHA1
      libzlib
    )
  endif (INCLUDE_VENDOR_LIBRARIES)
  target_link_libraries(
    crux
    barista
    bullseye
    hardklor
    cometsearch
    perclibrary
    fido
    xlink
    gflags
    crux-support
    tide-support
    mstoolkitlite
    protobuf
    ${PWIZ_LIBRARIES}
    ws2_32
    shlwapi
    optimized libboost_system-vc120-mt
    optimized libboost_iostreams-vc120-mt
    optimized libboost_nowide-vc120-mt
    optimized libboost_filesystem-vc120-mt
    optimized libboost_thread-vc120-mt
    optimized libboost_chrono-vc120-mt
    debug libboost_system-vc120-mt-gd
    debug libboost_iostreams-vc120-mt-gd
    debug libboost_nowide-vc120-mt-gd
    debug libboost_filesystem-vc120-mt-gd
    debug libboost_thread-vc120-mt-gd
    debug libboost_chrono-vc120-mt-gd
  )
else()
  target_link_libraries(
    crux
    xlink
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
    perclibrary
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
    boost_thread-mt-s
    boost_chrono
    boost_filesystem
    boost_system
    pthread
  )
endif(WIN32 AND NOT CYGWIN)

install (
  TARGETS
  crux
  DESTINATION
  bin
)

