# This is a CMake script for downloading the partial source distribution of 
# Google gflags from its GitHub repository.
#
# We have to use 'wget' for downloading because the GitHub site uses
# HTTPS rather than HTTP. HTTPS is not directly supported by the standard 
# build of CMake.
#
# This script expects two variables to be passed on the command line:
#   DOWNLOAD_DIR, the location to store the downloaded tar ball.
#   WORKING_DIR, the location where the expanded source files

# This macro checks download status codes for errors
macro (check_status download_name status_code error_message)
  if (${status_code} EQUAL 0)
    message(STATUS "Downloaded ${download_name}.")
  else (${status_code} EQUAL 0)
    message(
      FATAL_ERROR 
      "Unable to download ${download_name} "
      "\nError message: ${error_message}." 
      "\nDownload of Google gflags failed."
    )
  endif (${status_code} EQUAL 0)
endmacro (check_status)

# 'wget' has to be available
find_program(wget "wget")
if (${wget} STREQUAL "wget-NOTFOUND")
  message(
    FATAL_ERROR
    "The program 'wget' was not found.\n"
    "'wget' is required in order to download Google gflags."
  )
endif (${wget} STREQUAL "wget-NOTFOUND")
# 'tar' has to be available
find_program(tar "tar")
if (${tar} STREQUAL "tar-NOTFOUND")
  message(
    FATAL_ERROR
    "The program 'tar' was not found.\n"
    "'tar' is required in order to download Google gflags."
  )
endif (${tar} STREQUAL "tar-NOTFOUND")

set(download_url "https://codeload.github.com/gflags/gflags/tar.gz/v2.1.2")
set(filename "gflags-2.1.2.tar.gz")
execute_process(
  COMMAND ${wget} --no-check-certificate -nv -O ${DOWNLOAD_DIR}/${filename} "${download_url}"
  RESULT_VARIABLE status
  ERROR_VARIABLE error_message
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
check_status("Google gflags distribution" status ${error_message})

# Expand the tar file into the working directory
#  COMMAND ${tar} --strip-components=1 -xf ${DOWNLOAD_DIR}/${filename} 
execute_process(
  COMMAND ${CMAKE_COMMAND} -E tar xf ${DOWNLOAD_DIR}/${filename}
  RESULT_VARIABLE status
  WORKING_DIRECTORY ${WORKING_DIR}
)
if (${status} EQUAL 0)
  message(STATUS "Generated Google gflags source files.")
else (${status} EQUAL 0)
  message(FATAL_ERROR "Unable to generate Google gflags source files.")
endif (${status} EQUAL 0)
