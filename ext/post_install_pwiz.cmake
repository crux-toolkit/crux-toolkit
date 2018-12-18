# This is a CMake script for installing ProteoWizard for use in building crux.
# We have to rename one of the boost libraries and install the Thermo RAW
# reader DLL

# This script expects two variables to be passed on the command line:
#   BUILD_TYPE, the building configuration: Release or Debug
#   PREFIX, the location to store binaries from the build.

# This macro checks download status codes for errors
macro (check_status status_code)
  if (${status_code} EQUAL 0)
    message(STATUS "Fixed ProteoWizard install.")
  else ()
    message(
      FATAL_ERROR 
      "ProteoWizard install failed"
    )
  endif (${status_code} EQUAL 0)
endmacro (check_status)

if (WIN32 AND NOT CYGWIN)
  if (${BUILD_TYPE} MATCHES "Debug")
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_chrono-vc141-mt-gd.lib
        ${PREFIX}/lib/libboost_chrono-vc141-mt-gd-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_date_time-vc141-mt-gd.lib
        ${PREFIX}/lib/libboost_date_time-vc141-mt-gd-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_filesystem-vc141-mt-gd.lib
        ${PREFIX}/lib/libboost_filesystem-vc141-mt-gd-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_iostreams-vc141-mt-gd.lib
        ${PREFIX}/lib/libboost_iostreams-vc141-mt-gd-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_nowide-vc141-mt-gd.lib
        ${PREFIX}/lib/libboost_nowide-vc141-mt-gd-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_program_options-vc141-mt.lib
        ${PREFIX}/lib/libboost_program_options-vc141-mt-gd-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_system-vc141-mt-gd.lib
        ${PREFIX}/lib/libboost_system-vc141-mt-gd-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_thread-vc141-mt-gd.lib
        ${PREFIX}/lib/libboost_thread-vc141-mt-gd-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/build/src/ProteoWizard/pwiz_aux/msrc/utility/vendor_api/Bruker/x64/timsdata.lib
        ${PREFIX}/lib/timsdata.lib
      RESULT_VARIABLE status
    )
  else()
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_chrono-vc141-mt.lib
        ${PREFIX}/lib/libboost_chrono-vc141-mt-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_date_time-vc141-mt.lib
        ${PREFIX}/lib/libboost_date_time-vc141-mt-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_filesystem-vc141-mt.lib
        ${PREFIX}/lib/libboost_filesystem-vc141-mt-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_iostreams-vc141-mt.lib
        ${PREFIX}/lib/libboost_iostreams-vc141-mt-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_nowide-vc141-mt.lib
        ${PREFIX}/lib/libboost_nowide-vc141-mt-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_program_options-vc141-mt.lib
        ${PREFIX}/lib/libboost_program_options-vc141-mt-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_system-vc141-mt.lib
        ${PREFIX}/lib/libboost_system-vc141-mt-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_thread-vc141-mt.lib
        ${PREFIX}/lib/libboost_thread-vc141-mt-x64-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/build/src/ProteoWizard/pwiz_aux/msrc/utility/vendor_api/Bruker/x64/timsdata.lib
        ${PREFIX}/lib/timsdata.lib
      RESULT_VARIABLE status
    )
  endif (${BUILD_TYPE} MATCHES "Debug")
#  check_status(status)
endif (WIN32 AND NOT CYGWIN)
