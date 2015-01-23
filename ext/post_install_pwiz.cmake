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
        ${PREFIX}/lib/libboost_filesystem-vc100-mt-gd.lib
        ${PREFIX}/lib/libboost_filesystem-vc100-mt-gd-1_54.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_nowide-vc100-mt-gd.lib
        ${PREFIX}/lib/libboost_nowide-vc100-mt-gd-1_54.lib
      RESULT_VARIABLE status
    )
  else()
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_filesystem-vc100-mt.lib
        ${PREFIX}/lib/libboost_filesystem-vc100-mt-1_54.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_nowide-vc100-mt.lib
        ${PREFIX}/lib/libboost_nowide-vc100-mt-1_54.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_chrono-vc100-mt.lib
        ${PREFIX}/lib/libboost_chrono-vc100-mt-1_54.lib
      RESULT_VARIABLE status
    )
  endif (${BUILD_TYPE} MATCHES "Debug")
  check_status(status)
  execute_process(
    COMMAND regsvr32 
      /s
      ${PREFIX}/lib/MSFileReader.XRawfile2.dll
  )
  # Don't care about status of regsvr32 command
endif (WIN32 AND NOT CYGWIN)
