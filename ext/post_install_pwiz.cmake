# This is a CMake script for installing ProteoWizard for use in building crux.
# We have to rename one of the boost libraries because Windows is too clever
# and tries to link using the fully decorated BOOST library names,
# even though ProteoWizard builds the libraries using the undecorated names. 


# This script expects three variables to be passed on the command line:
#   BUILD_TYPE, the building configuration: Release or Debug
#   PREFIX, the location to store binaries from the build.
#   MSBUILD_PLATFORM, Win32 if for 32-bit build, x64 for 64-bit build

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
  # If windows build we need to select 32/64 bit platform
  if (MSBUILD_PLATFORM STREQUAL "x64")
    SET(BOOST_ARCH "x64")
    SET(BOOST_ARCH_DIR "x64")
  elseif (MSBUILD_PLATFORM STREQUAL "Win32")
    SET(BOOST_ARCH "x32")
    SET(BOOST_ARCH_DIR "x86")
  else ()
    message(FATAL_ERROR "${MSBUILD_PLATFORM}: Windows build architecture unknown.")
  endif (MSBUILD_PLATFORM STREQUAL "x64")
  # Select between debug/release builds
  if (${BUILD_TYPE} MATCHES "Debug")
    SET(TYPE "-gd")
  else()
    SET(TYPE "")
  endif (${BUILD_TYPE} MATCHES "Debug")

    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_chrono-vc141-mt${TYPE}.lib
        ${PREFIX}/lib/libboost_chrono-vc141-mt${TYPE}-${BOOST_ARCH}-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_date_time-vc141-mt${TYPE}.lib
        ${PREFIX}/lib/libboost_date_time-vc141-mt${TYPE}-${BOOST_ARCH}-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_filesystem-vc141-mt${TYPE}.lib
        ${PREFIX}/lib/libboost_filesystem-vc141-mt${TYPE}-${BOOST_ARCH}-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_iostreams-vc141-mt${TYPE}.lib
        ${PREFIX}/lib/libboost_iostreams-vc141-mt${TYPE}-${BOOST_ARCH}-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_nowide-vc141-mt${TYPE}.lib
        ${PREFIX}/lib/libboost_nowide-vc141-mt${TYPE}-${BOOST_ARCH}-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_program_options-vc141-mt${TYPE}.lib
        ${PREFIX}/lib/libboost_program_options-vc141-mt${TYPE}-${BOOST_ARCH}-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_system-vc141-mt${TYPE}.lib
        ${PREFIX}/lib/libboost_system-vc141-mt${TYPE}-${BOOST_ARCH}-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/lib/libboost_thread-vc141-mt${TYPE}.lib
        ${PREFIX}/lib/libboost_thread-vc141-mt${TYPE}-${BOOST_ARCH}-1_67.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/build/src/ProteoWizard/pwiz_aux/msrc/utility/vendor_api/Bruker/${BOOST_ARCH_DIR}/timsdata.lib
        ${PREFIX}/lib/timsdata.lib
      RESULT_VARIABLE status
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${PREFIX}/build/src/ProteoWizard/libraries/boost_1_67_0/boost/system/detail/error_code.ipp
        ${PREFIX}/include/boost/system/detail/error_code.ipp
      RESULT_VARIABLE status
    )
#  check_status(status)
endif (WIN32 AND NOT CYGWIN)
