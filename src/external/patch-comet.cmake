# This is a CMake script for patching Comet on Windows for use in the Crux
# build.

# This script expects two variables to be passed on the command line:
#   SOURCE_DIR, the path to the directory containing the patch source.
#   BINARY_DIR, the path to the directory containing the patch destination.

# This macro checks download status codes for errors
macro (check_status status_code)
  if (${status_code} EQUAL 0)
    message(STATUS "Patched Comet.")
  else ()
    message(
      FATAL_ERROR 
      "Comet patch failed"
    )
  endif (${status_code} EQUAL 0)
endmacro (check_status)

# This macro checks download status codes for errors
if (WIN32 AND NOT CYGWIN)
    message(${SOURCE_DIR}/patches/comet/CometSearch/CometSearch.vcxproj)
    message(${BINARY_DIR}/build/src/comet/CometSearch/CometSearch.vcxproj)
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${SOURCE_DIR}/patches/comet/CometSearch/CometSearch.vcxproj
        ${BINARY_DIR}/build/src/comet/CometSearch/CometSearch.vcxproj
        RESULT_VARIABLE status
    )
    check_status(status)
    message(${SOURCE_DIR}/patches/comet/Comet.sln)
    message(${BINARY_DIR}/build/src/comet/Comet.sln)
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${SOURCE_DIR}/patches/comet/Comet.sln
        ${BINARY_DIR}/build/src/comet/Comet.sln
        RESULT_VARIABLE status
    )
    check_status(status)
    message(${SOURCE_DIR}/patches/comet/MSToolkit/src/MSToolkit/MSReader.cpp)
    message(${BINARY_DIR}/build/src/comet/MSToolkit/src/MSToolkit/MSReader.cpp)
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${SOURCE_DIR}/patches/comet/MSToolkit/src/MSToolkit/MSReader.cpp
        ${BINARY_DIR}/build/src/comet/MSToolkit/src/MSToolkit/MSReader.cpp
        RESULT_VARIABLE status
    )
    check_status(status)
endif (WIN32 AND NOT CYGWIN)
