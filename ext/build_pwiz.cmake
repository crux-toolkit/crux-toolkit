# This is a CMake script for building ProteoWizard.
# We construct the command line for building ProteoWizard 
# based upond the platform (UNIX or Windows) and the build
# type (Release or Debug)

# This script expects two variables to be passed on the command line:
#   BUILD_TYPE, the building configuration: Release or Debug
#   PREFIX, the location to store binaries from the build.
#   WORKING_DIR, the location where the expanded source files can be found

# This macro checks download status codes for errors
macro (check_status status_code)
  if (${status_code} EQUAL 0)
    message(STATUS "Built ProteoWizard.")
  else ()
    message(
      FATAL_ERROR 
      "ProteoWizard build failed"
    )
  endif (${status_code} EQUAL 0)
endmacro (check_status)

# Set up arguments needed for ProteoWizard build.
# quickbuild.sh is the script provided by ProteoWizard on UNIX
# and quickbuild.bat is the corresponding script on Windows
set(pwiz_build_args ${pwiz_build_args} --prefix=${PREFIX})

# Get number of processors
if (EXISTS "/proc/cpuinfo")
  file(STRINGS "/proc/cpuinfo" procs REGEX "^processor.: [0-9]+$")
  list(LENGTH procs PROCESSOR_COUNT)
elseif (APPLE)
  find_program(cmd_sys_pro "sysctl")
  if (cmd_sys_pro)
    execute_process(COMMAND ${cmd_sys_pro} -n hw.ncpu OUTPUT_VARIABLE info)
    set(PROCESSOR_COUNT "${info}")
  endif()
elseif (WIN32)
  set(PROCESSOR_COUNT "$ENV{NUMBER_OF_PROCESSORS}")
endif()
if (PROCESSOR_COUNT)
  message("ProteoWizard build will use ${PROCESSOR_COUNT} threads")
  set(pwiz_build_args ${pwiz_build_args} -j${PROCESSOR_COUNT})
endif()

if (WIN32 AND NOT CYGWIN)
  set(pwiz_build "pwiz_tools\\build-apps.bat")
  if (MSBUILD_PLATFORM STREQUAL "x64")
    set(pwiz_build_args ${pwiz_build_args} 64)
  elseif (MSBUILD_PLATFORM STREQUAL "Win32")
    set(pwiz_build_args ${pwiz_build_args} 32)
  else (MSBUILD_PLATFORM STREQUAL "Win32")
    message(FATAL_ERROR "${MSBUILD_PLATFORM}: Windows build architecture unknown.")
  endif (MSBUILD_PLATFORM STREQUAL "x64")
  set(pwiz_build_args ${pwiz_build_args} --layout=versioned)
  set(pwiz_build_args ${pwiz_build_args} --link=static)
  set(pwiz_build_args ${pwiz_build_args} toolset=msvc-12.0)
  set(pwiz_build_args ${pwiz_build_args} --i-agree-to-the-vendor-licenses)
  set(pwiz_build_args ${pwiz_build_args} --without-mz5)
else()
  set(pwiz_build ./quickbuild.sh)
  set(pwiz_build_args ${pwiz_build_args} --without-binary-msdata)
  set(pwiz_build_args ${pwiz_build_args} --layout=system)
  set(pwiz_build_args ${pwiz_build_args} runtime-link=shared)
endif (WIN32 AND NOT CYGWIN)

if (${BUILD_TYPE} MATCHES "Debug")
  set(pwiz_build_args ${pwiz_build_args} variant=debug)
endif (${BUILD_TYPE} MATCHES "Debug")

set(pwiz_build_args ${pwiz_build_args} libraries)

execute_process(
  COMMAND ${pwiz_build} ${pwiz_build_args}
  RESULT_VARIABLE status
)
check_status(status)
