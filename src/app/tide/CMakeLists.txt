cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

###### Set up include and link directories  ######
if (WIN32 AND NOT CYGWIN)
include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/ext/include
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj
  ${CMAKE_BINARY_DIR}/ext/build/src/protobuf-2.5.0/src
  ${CMAKE_BINARY_DIR}/ext/build/src/google-gflags/src/windows
  ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_1_56_0
  ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_aux
)
else (WIN32 AND NOT CYGWIN)
include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/ext/include
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj
  ${CMAKE_BINARY_DIR}/ext/build/src/ProtocolBuffers/src
  ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_1_56_0
  ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_aux
)
endif (WIN32 AND NOT CYGWIN)

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
)
set(
  proto_files
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/header.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/peptides.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/raw_proteins.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/results.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/spectrum.proto
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
if (WIN32 AND NOT CYGWIN)
  set(
    tide_lib_files
    ${proto_files_compiled}
    abspath.cc
    active_peptide_queue.cc
    crux_sp_spectrum.cc
    fifo_alloc.cc
    index_settings.cc
    make_peptides.cc
    mass_constants.cc
    max_mz.cc
    mman.c
    peptide.cc
    peptide_mods3.cc
    peptide_peaks.cc
    sp_scorer.cc
    spectrum_collection.cc
    spectrum_preprocess2.cc
  )
else (WIN32 AND NOT CYGWIN)
  set(
    tide_lib_files
    ${proto_files_compiled}
    abspath.cc
    active_peptide_queue.cc
    crux_sp_spectrum.cc
    fifo_alloc.cc
    index_settings.cc
    make_peptides.cc
    mass_constants.cc
    max_mz.cc
    peptide.cc
    peptide_mods3.cc
    peptide_peaks.cc
    sp_scorer.cc
    spectrum_collection.cc
    spectrum_preprocess2.cc
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

