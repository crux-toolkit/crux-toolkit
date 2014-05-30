cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

###### Set up include and link directories  ######
if (WIN32 AND NOT CYGWIN)
include_directories(
  ${CMAKE_SOURCE_DIR}/src/c
  ${CMAKE_BINARY_DIR}/src/external/include
  ${CMAKE_BINARY_DIR}/src/external/build/src/ProtocolBuffers/src
  ${CMAKE_BINARY_DIR}/src/external/build/src/google-gflags/src/windows
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj
)
else (WIN32 AND NOT CYGWIN)
include_directories(
  ${CMAKE_SOURCE_DIR}/src/c
  ${CMAKE_BINARY_DIR}/src/external/include
  ${CMAKE_BINARY_DIR}/src/external/build/src/ProtocolBuffers/src
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj
)
endif (WIN32 AND NOT CYGWIN)
link_directories(
  ${CMAKE_BINARY_DIR}/src/external/lib
)

set_directory_properties(PROPERTIES
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
  ${CMAKE_CURRENT_BINARY_DIR}/protoobj/view_ms1_cmds.pb.cc
)
set(
  proto_files
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/header.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/peptides.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/raw_proteins.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/results.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/spectrum.proto
  ${CMAKE_CURRENT_SOURCE_DIR}/protobuf/view_ms1_cmds.proto
)
add_custom_command(
  OUTPUT ${proto_files_compiled}
  DEPENDS ${proto_files}
  COMMAND ${CMAKE_COMMAND} -E make_directory protoobj
  COMMAND
    ${CMAKE_BINARY_DIR}/src/external/bin/protoc 
      --proto_path=${CMAKE_CURRENT_SOURCE_DIR}/protobuf/ 
      --cpp_out=${CMAKE_CURRENT_BINARY_DIR}/protoobj
      --python_out=${CMAKE_CURRENT_BINARY_DIR}/protoobj 
      ${proto_files}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

########## Set up tide-support library  ##########
if (WIN32 AND NOT CYGWIN)
  set(tide_lib_files
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
    read_fasta.cc
    report.cc
    results_writer.cc
    sp_scorer.cc
    spectrum_collection.cc
    spectrum_preprocess2.cc
    weibull.cc
  )
else (WIN32 AND NOT CYGWIN)
  set(tide_lib_files
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
    read_fasta.cc
    report.cc
    results_writer.cc
    sp_scorer.cc
    spectrum_collection.cc
    spectrum_preprocess2.cc
    weibull.cc
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
endif (WIN32 AND NOT CYGWIN)

###### Build tide-test-spectrum-preprocess  ######
#add_executable(tide-test-spectrum-preprocess test_spectrum_preprocess.cc)
#set_property(TARGET tide-test-spectrum-preprocess PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-test-spectrum-preprocess tide-support gflags protobuf)

########### Build tide-peptide-sorter  ###########
#add_executable(tide-peptide-sorter peptide_sorter.cc)
#set_property(TARGET tide-peptide-sorter PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-peptide-sorter tide-support gflags protobuf)

######### Build tide-peptide-downselect  #########
#add_executable(tide-peptide-downselect peptide_downselect.cc)
#set_property(TARGET tide-peptide-downselect PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-peptide-downselect tide-support gflags protobuf)

############ Build tide-make-peptides ############
#add_executable(tide-make-peptides make_peptides_cmdline.cc)
#set_property(TARGET tide-make-peptides PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-make-peptides tide-support gflags protobuf)

########### Build tide-filter-peptides ###########
#add_executable(tide-filter-peptides filter_peptides.cc)
#set_property(TARGET tide-filter-peptides PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-filter-peptides tide-support gflags protobuf)

########### Build tide-peptide-counter ###########
#add_executable(tide-peptide-counter peptide_counter.cc)
#set_property(TARGET tide-peptide-counter PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-peptide-counter tide-support gflags protobuf)

############### Build tide-search  ###############
add_executable(tide-search search.cc)
if (WIN32 AND NOT CYGWIN)
  set_property(
    TARGET tide-search 
    PROPERTY 
      COMPILE_DEFINITIONS 
      MAIN
      GFLAGS_DLL_DECL= 
      GFLAGS_DLL_DECLARE_FLAG= 
      GFLAGS_DLL_DEFINE_FLAG=
  )
else (WIN32 AND NOT CYGWIN)
  set_property(TARGET tide-search PROPERTY COMPILE_DEFINITIONS MAIN)
endif (WIN32 AND NOT CYGWIN)
target_link_libraries(tide-search tide-support gflags protobuf)

############ Build tide-test-read-ms2 ############
#add_executable(tide-test-read-ms2 test_read_ms2.cc)
#set_property(TARGET tide-test-read-ms2 PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-test-read-ms2 tide-support gflags protobuf)

######## Build tide-read-spectrumrecords  ########
#add_executable(tide-read-spectrumrecords read_spectrumrecords.cc)
#set_property(TARGET tide-read-spectrumrecords PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-read-spectrumrecords tide-support gflags protobuf)

########### Build tide-filter-spectra  ###########
#add_executable(tide-filter-spectra filter_spectra.cc)
#set_property(TARGET tide-filter-spectra PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-filter-spectra tide-support gflags protobuf)

########### Build tide-combine-spectra ###########
#add_executable(tide-combine-spectra combine_spectra.cc)
#set_property(TARGET tide-combine-spectra PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-combine-spectra tide-support gflags protobuf)

############ Build tide-test-records  ############
#add_executable(tide-test-records test_records.cc)
#set_property(TARGET tide-test-records PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-test-records tide-support gflags protobuf)

############# Build tide-read-fasta  #############
#add_executable(tide-read-fasta read_fasta_cmdline.cc)
#set_property(TARGET tide-read-fasta PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-read-fasta tide-support gflags protobuf)

####### Build tide-test-theoretical-peaks  #######
#add_executable(tide-test-theoretical-peaks test_theoretical_peaks.cc)
#set_property(TARGET tide-test-theoretical-peaks PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-test-theoretical-peaks tide-support gflags protobuf)

########### Build tide-test-fifo-alloc ###########
#add_executable(tide-test-fifo-alloc test_fifo_alloc.cc)
#set_property(TARGET tide-test-fifo-alloc PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-test-fifo-alloc tide-support)

############# Build tide-test-order  #############
#add_executable(tide-test-order test_order.cc)
#set_property(TARGET tide-test-order PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-test-order tide-support gflags protobuf)

######### Build tide-read-peptide-diffs  #########
#add_executable(tide-read-peptide-diffs read_peptide_diffs.cc)
#set_property(TARGET tide-read-peptide-diffs PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-read-peptide-diffs tide-support gflags protobuf)

############ Build tide-peptide-mods  ############
#add_executable(tide-peptide-mods peptide_mods.cc)
#set_property(TARGET tide-peptide-mods PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-peptide-mods tide-support gflags protobuf)

############ Build tide-peptide-mods2 ############
#add_executable(tide-peptide-mods2 peptide_mods2.cc)
#set_property(TARGET tide-peptide-mods2 PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-peptide-mods2 tide-support gflags protobuf)

########### Build tide-test-mod-table  ###########
#add_executable(tide-test-mod-table test_mod_table.cc)
#set_property(TARGET tide-test-mod-table PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-test-mod-table tide-support gflags protobuf)

########## Build tide-test-peptide-mods ##########
#add_executable(tide-test-peptide-mods test_peptide_mods.cc)
#set_property(TARGET tide-test-peptide-mods PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-test-peptide-mods tide-support gflags protobuf)

############ Build tide-peptide-peaks ############
#add_executable(tide-peptide-peaks peptide_peaks_cmdline.cc)
#set_property(TARGET tide-peptide-peaks PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-peptide-peaks tide-support gflags protobuf)

################ Build tide-index ################
add_executable(tide-index index.cc)
set_property(TARGET tide-index PROPERTY COMPILE_DEFINITIONS MAIN)
target_link_libraries(tide-index tide-support gflags protobuf)

############## Build tide-view-ms1  ##############
#add_executable(tide-view-ms1 view_ms1.cc)
#set_property(TARGET tide-view-ms1 PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-view-ms1 tide-support gflags protobuf)

############## Build tide-view-ms12 ##############
#add_executable(tide-view-ms12 view_ms12.cc)
#set_property(TARGET tide-view-ms12 PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-view-ms12 tide-support gflags protobuf)

############### Build tide-results ###############
#add_executable(tide-results results.cc)
#set_property(TARGET tide-results PROPERTY COMPILE_DEFINITIONS MAIN)
#target_link_libraries(tide-results tide-support gflags protobuf)

