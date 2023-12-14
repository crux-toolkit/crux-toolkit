# CMakeList.txt : CMake project for crux-lfq, include source and define
# project specific logic here.
#
cmake_minimum_required(VERSION 3.15)

project ("crux")

set(TEST_DATA_SOURCE_DIR ${CMAKE_SOURCE_DIR}/test/unit-tests/test-data)
set(TEST_DATA_DESTINATION_DIR "${CMAKE_BINARY_DIR}/test/unit-tests")
file(COPY ${TEST_DATA_SOURCE_DIR} DESTINATION ${TEST_DATA_DESTINATION_DIR})
add_definitions(-DTEST_DATA_PATH="${TEST_DATA_DESTINATION_DIR}")


macro(use_cxx11)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  else (CMAKE_VERSION VERSION_LESS "3.1")
    set(CMAKE_CXX_STANDARD 11)
  endif (CMAKE_VERSION VERSION_LESS "3.1")
endmacro(use_cxx11)

use_cxx11()


enable_testing()

include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1
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
include_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator/src/fido)
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
link_directories(${CMAKE_BINARY_DIR}/ext/build/src/percolator-build/src/fido)
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

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(crux_lfq_test crux_lfq_test.cpp)

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
  fido
  picked_protein
  protobuf
  pwiz_data_msdata
  pwiz_data_msdata_mz5
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
)
include(GoogleTest)
gtest_discover_tests(crux_lfq_test)