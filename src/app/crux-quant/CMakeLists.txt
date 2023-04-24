cmake_minimum_required(VERSION 3.15)
cmake_policy(VERSION 3.15)

###### Set up include and link directories  ######
if (WIN32 AND NOT CYGWIN)
  include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/ext/include
    ${CMAKE_CURRENT_BINARY_DIR}/protoobj
    ${CMAKE_BINARY_DIR}/ext/build/src/protobuf-2.5.0/src
    ${CMAKE_BINARY_DIR}/ext/build/src/google-gflags/src/windows
    ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_1_76_0
    ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_aux
  )
  # On windows set static linking of runtime libraries
  set(
    CMAKE_MSVC_RUNTIME_LIBRARY 
    "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL"
  )
else (WIN32 AND NOT CYGWIN)
include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/ext/include
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj
  ${CMAKE_BINARY_DIR}/ext/build/src/ProtocolBuffers/src
  ${CMAKE_BINARY_DIR}/ext/build/src/ProteoWizard/libraries/boost_1_76_0
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
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj/indexed_mass_spectral_peak.pb.cc
)
set(
  proto_files
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/indexed_mass_spectral_peak.proto
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

########## Set up crux-quant-support library  ##########
if (WIN32 AND NOT CYGWIN)
  set(
    crux_quant_lib_files
    ${proto_files_compiled}
    IndexedMassSpectralPeak.cpp
    PeakIndexingEngine.cpp
  )
else (WIN32 AND NOT CYGWIN)
  set(
    crux_quant_lib_files
    ${proto_files_compiled}
    IndexedMassSpectralPeak.cpp
    PeakIndexingEngine.cpp
  )
endif (WIN32 AND NOT CYGWIN)
add_library(crux-quant-support STATIC ${crux_quant_lib_files})

if (WIN32 AND NOT CYGWIN)
  set_property(
    TARGET crux-quant-support 
    PROPERTY 
      COMPILE_DEFINITIONS 
      GFLAGS_DLL_DECL= 
      GFLAGS_DLL_DECLARE_FLAG= 
      GFLAGS_DLL_DEFINE_FLAG=
  )
elseif (APPLE)
  set_property(
    TARGET crux-quant-support 
    PROPERTY 
      COMPILE_DEFINITIONS 
      DARWIN 
  )
endif (WIN32 AND NOT CYGWIN)
