# CMakeList.txt : CMake project for crux-lfq, include source and define
# project specific logic here.
#
cmake_minimum_required(VERSION 3.15)

project("crux")

set(TEST_DATA_SOURCE_DIR ${CMAKE_SOURCE_DIR}/test/unit-tests/test-data)
set(TEST_DATA_DESTINATION_DIR "${CMAKE_BINARY_DIR}/test/unit-tests")
file(COPY ${TEST_DATA_SOURCE_DIR} DESTINATION ${TEST_DATA_DESTINATION_DIR})
add_definitions(-DTEST_DATA_PATH="${TEST_DATA_DESTINATION_DIR}")

if(UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-pie")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
elseif(WIN32)
  # disable autolinking in boost
  add_definitions(-DBOOST_ALL_NO_LIB)
endif(UNIX)

enable_testing()

include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
)

include_directories(${CMAKE_SOURCE_DIR}/ext)
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src/app/bullseye)
include_directories(${CMAKE_SOURCE_DIR}/src/app/hardklor)
include_directories(${CMAKE_SOURCE_DIR}/src/app/crux-lfq)
include_directories(${CMAKE_SOURCE_DIR}/src/app/tide)
include_directories(${CMAKE_SOURCE_DIR}/src/app)

include_directories(${CMAKE_SOURCE_DIR}/src/io)
include_directories(${CMAKE_SOURCE_DIR}/src/model)
include_directories(${CMAKE_SOURCE_DIR}/src/util)

include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_BINARY_DIR}/ext)
include_directories(${CMAKE_BINARY_DIR}/ext)
include_directories(${CMAKE_BINARY_DIR}/ext/include)
include_directories(${CMAKE_BINARY_DIR}/ext/include/CometSearch)
include_directories(BEFORE ${CMAKE_BINARY_DIR}/ext/include/MSToolkit)
include_directories(${CMAKE_BINARY_DIR}/ext/include/kojak)
include_directories(${CMAKE_BINARY_DIR}/ext/include/gflags)
include_directories(${CMAKE_BINARY_DIR}/ext/include/google/protobuf)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_aux)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator/src/)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator/src/blas)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator/src/picked_protein)
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src/)
include_directories(${CMAKE_BINARY_DIR}/src/app/tide/protoobj)

link_directories(${CMAKE_BINARY_DIR}/src)
link_directories(${CMAKE_BINARY_DIR}/src/app/bullseye)
link_directories(${CMAKE_BINARY_DIR}/src/app/hardklor)
link_directories(${CMAKE_BINARY_DIR}/src/app/tide/protoobj)
link_directories(${CMAKE_BINARY_DIR}/ext/lib)
link_directories(${CMAKE_BINARY_DIR}/ext/lib64)
link_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src)
link_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src/blas)
link_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src/picked_protein)

if(WIN32 AND NOT Cygwin)
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

  # Need Comet base directory for Windows
  include_directories(${CMAKE_BINARY_DIR}/ext/build/src/comet/)

  # On Windows we should build vendor libraries by default
  if(NOT DEFINED INCLUDE_VENDOR_LIBRARIES)
    set(INCLUDE_VENDOR_LIBRARIES TRUE)
  endif(NOT DEFINED INCLUDE_VENDOR_LIBRARIES)

  # On windows set static linking of runtime libraries
  set(
    CMAKE_MSVC_RUNTIME_LIBRARY
    "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL"
  )
endif(WIN32 AND NOT Cygwin)

if(WIN32 AND NOT Cygwin)
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
endif(WIN32 AND NOT Cygwin)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(
    crux_lib_files
    ${crux_lib_files}
    util/WinCrux.cpp
  )
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(crux_lfq_test crux_lfq_test.cpp)

if(WIN32 AND NOT CYGWIN)
  if(INCLUDE_VENDOR_LIBRARIES)
    set(
      PWIZ_LIBRARIES
      libhdf5.lib
      libhdf5pp.lib
      libpwiz_analysis_chromatogram_processing.lib
      libpwiz_analysis_common.lib
      libpwiz_analysis_demux.lib
      libpwiz_analysis_frequency.lib
      libpwiz_analysis_passive.lib
      libpwiz_analysis_passive_core.lib
      libpwiz_analysis_peakdetect.lib
      libpwiz_analysis_peptideid.lib
      libpwiz_analysis_proteome_processing.lib
      libpwiz_analysis_spectrum_processing.lib
      libpwiz_analysis_spectrum_processing_basic.lib
      libpwiz_analysis_version.lib
      libpwiz_data_common.lib
      libpwiz_data_identdata.lib
      libpwiz_data_identdata_examples.lib
      libpwiz_data_identdata_version.lib
      libpwiz_data_misc.lib
      libpwiz_data_msdata.lib
      libpwiz_data_msdata_core.lib
      libpwiz_data_msdata_examples.lib
      libpwiz_data_msdata_mz5.lib
      libpwiz_data_msdata_version.lib
      libpwiz_data_obo.lib
      libpwiz_data_proteome.lib
      libpwiz_data_proteome_examples.lib
      libpwiz_data_proteome_version.lib
      libpwiz_data_vendor_readers.lib
      libpwiz_reader_abi.lib
      libpwiz_reader_abi_t2d.lib
      libpwiz_reader_agilent.lib
      libpwiz_reader_bruker.lib
      libpwiz_reader_shimadzu.lib
      libpwiz_reader_thermo.lib
      libpwiz_reader_uimf.lib
      libpwiz_reader_unifi.lib
      libpwiz_reader_waters.lib
      libpwiz_utility_chemistry.lib
      libpwiz_utility_image.lib
      libpwiz_utility_math.lib
      libpwiz_utility_minimxml.lib
      libpwiz_utility_misc.lib
      libpwiz_utility_vendor_reader_test_harness.lib
      libpwiz_vendor_api_agilent.lib
      libpwiz_vendor_api_bruker.lib
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
  else(INCLUDE_VENDOR_LIBRARIES)
    set(
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
      libpwiz_utility_image
      libpwiz_utility_math
      libpwiz_utility_minimxml
      libpwiz_version
      libSHA1
      libzlib
    )
  endif(INCLUDE_VENDOR_LIBRARIES)

  SET(
    BOOST_LIBRARIES
    optimized libboost_chrono-vc142-mt
    debug libboost_chrono-vc142-mt-gd
    optimized libboost_filesystem-vc142-mt
    debug libboost_filesystem-vc142-mt-gd
    optimized libboost_iostreams-vc142-mt
    debug libboost_iostreams-vc142-mt-gd
    optimized libboost_nowide-vc142-mt
    debug libboost_nowide-vc142-mt-gd
    optimized libboost_serialization-vc142-mt
    debug libboost_serialization-vc142-mt-gd
    optimized libboost_system-vc142-mt
    debug libboost_system-vc142-mt-gd
    optimized libboost_thread-vc142-mt
    debug libboost_thread-vc142-mt-gd
  )
  target_link_libraries(
    crux_lfq_test
    GTest::gtest
    GTest::gtest_main
    bullseye
    hardklor
    cometsearch
    kojaksearch
    crux-lfq-support
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
    optimized libfreetype-vc142-mt-2_4
    debug libfreetype-vc142-mt-gd-2_4
    optimized libgd-vc142-mt-2_1
    debug libgd-vc142-mt-gd-2_1
    libpng
  )
else()
  # UNIX SYSTEMS
  target_link_libraries(
    crux_lfq_test
    GTest::gtest
    GTest::gtest_main
    crux-support
    tide-support
    cometsearch
    kojaksearch
    crux-lfq-support
    libneopepxmlparser.a
    gflags
    bullseye
    hardklor
    mstoolkitlite
    perclibrary
    blas
    picked_protein
    protobuf
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
    dl
    m
    pthread
    hdf5pp
    hdf5
    ${RT_LIBRARY}
    SHA1
    zlib
    ssl
    crypto
  )
endif(WIN32 AND NOT CYGWIN)

include(GoogleTest)
gtest_discover_tests(crux_lfq_test)
