cmake_minimum_required(VERSION 3.15)
cmake_policy(VERSION 3.15)

###### Set up include and link directories  ######
include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/ext/include
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj
  ${CMAKE_BINARY_DIR}/ext/build/src/ProtocolBuffers/src
  ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/
  ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_aux
)
if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-pie")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
elseif (WIN32 AND NOT CYGWIN)
  # On windows set static linking of runtime libraries
  set(
    CMAKE_MSVC_RUNTIME_LIBRARY 
    "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL"
  )
  # disable autolinking in boost
  add_definitions( -DBOOST_ALL_NO_LIB )
endif(UNIX)

link_directories(
  ${CMAKE_BINARY_DIR}/ext/lib
)

set_directory_properties(
  PROPERTIES
  ADDITIONAL_MAKE_CLEAN_FILES
  protoobj
)

########### For compiling proto files  ###########
set(
  proto_files_compiled
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj/header.pb.cc
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj/peptides.pb.cc
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj/raw_proteins.pb.cc
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj/results.pb.cc
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj/spectrum.pb.cc
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj/residue_stats.pb.cc
)
set(
  proto_files
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/header.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/peptides.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/raw_proteins.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/results.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/spectrum.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/residue_stats.proto
)
add_custom_command(
  OUTPUT ${proto_files_compiled}
  DEPENDS ${proto_files}
  COMMAND ${CMAKE_COMMAND} -E make_directory protoobj
  COMMAND
    ${CMAKE_BINARY_DIR}/ext/bin/protoc 
      --proto_path=${CMAKE_CURRENT_SOURCE_DIR}/protobuf/ 
      --cpp_out=${CMAKE_CURRENT_BINARY_DIR}/protoobj
      --python_out=${CMAKE_CURRENT_BINARY_DIR}/protoobj 
      ${proto_files}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

########## Set up tide-support library  ##########
set(
  tide_lib_files
  ${proto_files_compiled}
  abspath.cc
  ActivePeptideQueue.cc  
  crux_sp_spectrum.cc
  fifo_alloc.cc
  index_settings.cc
  make_peptides.cc
  mass_constants.cc
  max_mz.cc
  peptide.cc
  peptide_mods3.cc
  peptide_peaks.cc
  spectrum_collection.cc
  spectrum_preprocess2.cc
)
if (WIN32 AND NOT CYGWIN)
  set(
    tide_lib_files
    ${tide_lib_files}
    mman.c
  )
endif (WIN32 AND NOT CYGWIN)
add_library(tide-support STATIC ${tide_lib_files})

if (WIN32 AND NOT CYGWIN)
  set_property(
    TARGET tide-support 
    PROPERTY 
      COMPILE_DEFINITIONS 
      GFLAGS_DLL_DECL= 
      GFLAGS_DLL_DECLARE_FLAG= 
      GFLAGS_DLL_DEFINE_FLAG=
  )
elseif (APPLE)
  set_property(
    TARGET tide-support 
    PROPERTY 
      COMPILE_DEFINITIONS 
      DARWIN 
  )
endif (WIN32 AND NOT CYGWIN)

