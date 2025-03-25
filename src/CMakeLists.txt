cmake_minimum_required(VERSION 3.15)
cmake_policy(VERSION 3.15)


if (APPLE)
    enable_language(OBJC)
    find_library(FOUNDATION Foundation)
    set (ASYNCPOST util/AsyncPost.mm)
endif(APPLE)
if (NOT APPLE AND NOT WIN32)
  set(
    SSL_SUPPORT
    crypto
    ssl
  )
endif (NOT APPLE AND NOT WIN32)

add_subdirectory(app/bullseye)
add_subdirectory(app/hardklor)
add_subdirectory(app/tide)

if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-pie")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
elseif (WIN32)
  # disable autolinking in boost
  add_definitions( -DBOOST_ALL_NO_LIB )
endif(UNIX)

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src/app/bullseye)
include_directories(${CMAKE_SOURCE_DIR}/src/app/hardklor)
include_directories(${CMAKE_SOURCE_DIR}/src/app/tide/protoobj)
include_directories(${CMAKE_SOURCE_DIR}/src/io)
include_directories(${CMAKE_SOURCE_DIR}/src/model)
include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_BINARY_DIR}/ext/include)
include_directories(${CMAKE_BINARY_DIR}/ext/include/CometSearch)
include_directories(${CMAKE_BINARY_DIR}/ext/include/MSToolkit)
include_directories(${CMAKE_BINARY_DIR}/ext/include/kojak)
include_directories(${CMAKE_BINARY_DIR}/ext/include/gflags)
include_directories(${CMAKE_BINARY_DIR}/ext/include/google/protobuf)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_aux)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/pwiz)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator/src/)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator/src/blas)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator/src/picked_protein)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src/)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/app/tide/protoobj)
if (WIN32 AND NOT Cygwin)
  # Needed to put DLL containing type libraries
  # in include path for Windows
  include_directories(${CMAKE_BINARY_DIR}/ext/lib)
  SET(
    CMAKE_EXE_LINKER_FLAGS 
    "${CMAKE_EXE_LINKER_FLAGS} /VERBOSE:LIB "
  )
  SET(
    CMAKE_EXE_LINKER_FLAGS_DEBUG 
    "${CMAKE_EXE_LINKER_FLAGS} /VERBOSE:LIB "
  )
  #  Need Comet base directory for Windows
  include_directories(${CMAKE_BINARY_DIR}/ext/build/src/comet/)

  # On Windows we should build vendor libraries by default
  if (NOT DEFINED INCLUDE_VENDOR_LIBRARIES)
    set(INCLUDE_VENDOR_LIBRARIES TRUE)
  endif (NOT DEFINED INCLUDE_VENDOR_LIBRARIES)

  # On windows set static linking of runtime libraries
  set(
    CMAKE_MSVC_RUNTIME_LIBRARY 
    "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL"
  )
endif (WIN32 AND NOT Cygwin)

link_directories(${CMAKE_BINARY_DIR}/app/bullseye)
link_directories(${CMAKE_BINARY_DIR}/app/hardklor)
link_directories(${CMAKE_BINARY_DIR}/app/tide-support)
link_directories(${CMAKE_BINARY_DIR}/ext/lib)
link_directories(${CMAKE_BINARY_DIR}/ext/lib64)
link_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src) 
link_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src/blas)
link_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src/picked_protein)
if (WIN32 AND NOT Cygwin)
 # Needed for windows build of vendor libraries
  link_directories(
    ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/pwiz_aux/msrc/utility/vendor_api/Waters/vc12_x64
    ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/pwiz_aux/msrc/utility/vendor_api/Bruker/x64
  )
  SET(
    CMAKE_SHARED_LINKER_FLAGS 
    ${CMAKE_SHARED_LINKER_FLAGS} 
    "/NODEFAULTLIBS:MSVCRT" 
  )
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
  app/CometIndexApplication.cpp
  app/ComputeQValues.cpp
  app/CreateDocs.cpp
  ${ASYNCPOST}
  util/crux-utils.cpp
  model/Database.cpp
  model/DatabaseProteinIterator.cpp
  model/DatabasePeptideIterator.cpp
  app/CruxApplication.cpp
  app/CruxApplicationList.cpp
  io/DelimitedFile.cpp
  io/DelimitedFileReader.cpp
  io/DelimitedFileWriter.cpp
  util/FileUtils.cpp
  model/FilteredSpectrumChargeIterator.cpp
  model/GeneratePeptidesIterator.cpp
  app/GetMs2Spectrum.cpp
  io/HTMLWriter.cpp
  model/Ion.cpp
  model/IonConstraint.cpp
  model/IonFilteredIterator.cpp
  model/IonSeries.cpp
  app/KojakApplication.cpp
  util/linked_list.cpp
  io/LineFileReader.cpp
  app/LocalizeModification.cpp
  util/mass.cpp
  app/MakePinApplication.cpp
  model/Match.cpp
  io/MatchColumns.cpp
  io/MatchFileReader.cpp
  io/MatchFileWriter.cpp
  model/MatchCollection.cpp
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
  app/ParamMedicInference.cpp
  util/Params.cpp
  model/Peak.cpp
  model/Peptide.cpp
  model/PeptideConstraint.cpp
  util/peptide_modifications.cpp
  model/PeptideSrc.cpp
  io/PepXMLReader.cpp
  io/PepXMLWriter.cpp
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
  model/Scorer.cpp
  app/SpectralCounts.cpp
  io/SpectrumCollection.cpp
  io/SpectrumCollectionFactory.cpp
  model/Spectrum.cpp
  io/SpectrumRecordSpectrumCollection.cpp
  io/SpectrumRecordWriter.cpp
  model/SpectrumZState.cpp
  util/StringUtils.cpp
  io/SQTReader.cpp
  io/SQTWriter.cpp
  app/TideIndexApplication.cpp
  app/TideMatchSet.cpp
  app/SpectrumConvertApplication.cpp
  app/TideSearchApplication.cpp
  io/DIAmeterFeatureScaler.cpp
  io/DIAmeterPSMFilter.cpp
  io/DIAmeterCVSelector.cpp
  app/DIAmeterApplication.cpp
  util/utils.cpp
)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(
    crux_lib_files
    ${crux_lib_files}
    util/WinCrux.cpp
  )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

add_definitions(-DCRUX)
add_definitions(-DBOOST_ALL_NO_LIB)

add_library(
  crux-support
  STATIC
  ${crux_lib_files}
)

if (APPLE)
  set(OPENMP_FOUND FALSE)
else(APPLE)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif(OPENMP_FOUND)
endif(APPLE)

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
	  _DEFAULT_SOURCE
	  _POSIX_C_SOURCE=199506L
	  _FILE_OFFSET_BITS=64
	  _GNU_SOURCE
	  GCC
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

if (WIN32 AND NOT CYGWIN)
  if (INCLUDE_VENDOR_LIBRARIES)
    set (
      PWIZ_LIBRARIES
      libhdf5.lib
      libhdf5pp.lib
      libpwiz_data_common.lib
      libpwiz_data_identdata.lib
      libpwiz_data_identdata_examples.lib
      libpwiz_data_identdata_version.lib
      libpwiz_data_misc.lib
      libpwiz_data_msdata.lib
      libpwiz_data_msdata_core.lib
      libpwiz_data_msdata_examples.lib
      libpwiz_data_msdata_mz5.lib
      libpwiz_data_msdata_mzmlb.lib
      libpwiz_data_msdata_version.lib
      libpwiz_data_obo.lib
      libpwiz_data_proteome.lib
      libpwiz_data_proteome_examples.lib
      libpwiz_data_proteome_version.lib
      libpwiz_data_vendor_readers.lib
      libpwiz_reader_abi.lib
      libpwiz_reader_abi_t2d.lib
      libpwiz_reader_agilent.lib
      libpwiz_reader_shimadzu.lib
      libpwiz_reader_thermo.lib
      libpwiz_reader_uimf.lib
      libpwiz_reader_unifi.lib
      libpwiz_reader_waters.lib
      libpwiz_utility_chemistry.lib
      libpwiz_utility_math.lib
      libpwiz_utility_minimxml.lib
      libpwiz_utility_misc.lib
      libpwiz_utility_vendor_reader_test_harness.lib
      libpwiz_vendor_api_agilent.lib
      libpwiz_vendor_api_shimadzu.lib
      libpwiz_vendor_api_thermo.lib
      libpwiz_vendor_api_uimf.lib
      libpwiz_vendor_api_unifi.lib
      libpwiz_version.lib
      libSHA1
      libsqlite3.lib
      libsqlite3pp.lib
      libsvm
      libzlib
      baf2sql_c
      timsdata
    )
  else (INCLUDE_VENDOR_LIBRARIES)
    set (
      PWIZ_LIBRARIES
      libpwiz_data_identdata_examples
      libpwiz_data_identdata_version
      libpwiz_data_misc
      libpwiz_data_msdata
      libpwiz_data_msdata_core
      libpwiz_data_msdata_examples
      libpwiz_data_msdata_version
      libpwiz_data_proteome_examples
      libpwiz_data_proteome_version
      libpwiz_data_common
      libpwiz_utility_chemistry
      libpwiz_utility_math
      libpwiz_utility_minimxml
      libpwiz_version
      libSHA1
      libzlib
    )
  endif (INCLUDE_VENDOR_LIBRARIES)
  SET(
    BOOST_LIBRARIES
    optimized libboost_chrono-mt
    debug libboost_chrono-mt-gd
    optimized libboost_filesystem-mt
    debug libboost_filesystem-mt-gd
    optimized libboost_iostreams-mt
    debug libboost_iostreams-mt-gd
    optimized libboost_nowide-mt
    debug libboost_nowide-mt-gd
    optimized libboost_serialization-mt
    debug libboost_serialization-mt-gd
    optimized libboost_system-mt
    debug libboost_system-mt-gd
    optimized libboost_thread-mt.lib
    debug libboost_thread-vc142-mt-gd
  )
  target_link_libraries(
    crux
    bullseye
    hardklor
    cometsearch
    kojaksearch
    neopepxmlparser.lib
    perclibrary
    blas
    picked_protein
    optimized gflags_static
    debug gflags_static_debug
    crux-support
    tide-support
    mstoolkitlite
    optimized libprotobuf
    debug libprotobufd
    ${PWIZ_LIBRARIES}
    ws2_32
    shlwapi
    ${BOOST_LIBRARIES}
  )
else()
  # UNIX SYSTEMS
  target_link_libraries(
    crux
    crux-support
    tide-support
    cometsearch
    kojaksearch
    libneopepxmlparser.a
    gflags
    bullseye
    hardklor
    mstoolkitlite
    perclibrary
    blas
    picked_protein
    pwiz_data_msdata
    pwiz_data_msdata_mz5
    pwiz_data_msdata_mzmlb
    pwiz_data_msdata_version
    pwiz_data_msdata_core
    pwiz_data_identdata
    pwiz_data_identdata_version
    pwiz_data_proteome
    pwiz_data_proteome_version
    pwiz_data_misc
    pwiz_data_common
    pwiz_utility_chemistry
    pwiz_utility_minimxml
    pwiz_utility_misc
    pwiz_version
    boost_chrono
    boost_filesystem
    boost_iostreams
    boost_serialization
    boost_system
    boost_thread
    protobuf
    protoc
    dl
    m
    ${SSL_SUPPORT}
    pthread
    hdf5pp
    hdf5
    SHA1
    zlib
    ${FOUNDATION}
  )
endif(WIN32 AND NOT CYGWIN)

install (
  TARGETS
  crux
  DESTINATION
  bin
)

