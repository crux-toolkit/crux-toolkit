# This is a CMake script for downloading the partial source distribution of 
# Proteowizard from its Team City repository.
#
# We have to use 'wget' for downloading because the Proteowizard site uses
# HTTPS rather than HTTP. HTTPS is not directly supported by the standard 
# build of CMake.
#
# To get the URL for the distribution tarball we have to find the build id, 
# then use the build id to look up the Proteowizard version number.
# Once we have the version number we can build the URL for downloading the distribution.

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
      "\nDownload of Proteowizard failed."
    )
  endif (${status_code} EQUAL 0)
endmacro (check_status)

# 'wget' has to be available
find_program(wget "wget")
if (${wget} STREQUAL "wget-NOTFOUND")
  message(
    FATAL_ERROR
    "The program 'wget' was not found.\n"
    "'wget' is required in order to download ProteoWizard."
  )
endif (${wget} STREQUAL "wget-NOTFOUND")

# Team City performs several different builds. We want the one
# for Pwiz source without Skyline but with vendor library suppport.
set(build_type "bt81")

set(
  build_info_url 
  https://teamcity.labkey.org:/app/rest/buildTypes/id:${build_type}/builds?status=SUCCESS&count=1&guest=1
)
execute_process(
  COMMAND ${wget} --no-check-certificate -nv -O build.info.txt "${build_info_url}"
  RESULT_VARIABLE status
  OUTPUT_VARIABLE build_info 
  ERROR_VARIABLE error_message
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
check_status("build info" status ${error_message})
file (STRINGS "build.info.txt" build_info)

# Using the build id download the version string
string(REGEX REPLACE "^.*build id=\"([0-9]+)\".*$" "\\1" build_id "${build_info}")
set(
  build_info_url https://teamcity.labkey.org/repository/download/${build_type}/${build_id}:id/VERSION?guest=1
)
execute_process(
  COMMAND ${wget} --no-check-certificate -nv -O version.info.txt "${build_info_url}"
  RESULT_VARIABLE status
  OUTPUT_VARIABLE version_info 
  ERROR_VARIABLE error_message
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
check_status("version string" status ${error_message})
file (STRINGS "version.info.txt" version_info)

# Use the version string to build the URL
set(
  base_download_url
  "https://teamcity.labkey.org:/guestAuth/repository/download/${build_type}/.lastSuccessful/"
)
string(REGEX REPLACE "^([0-9]+)\\.([0-9]+)\\.([0-9]+).*$" "\\1_\\2_\\3" version_id ${version_info})
set(filename pwiz-src-${version_id}.tar.bz2)
set(
  download_url 
  "${base_download_url}${filename}"
)
# Using the version string download the file
execute_process(
  COMMAND ${wget} --no-check-certificate -nv -O ${DOWNLOAD_DIR}/${filename} "${download_url}"
  RESULT_VARIABLE status
  ERROR_VARIABLE error_message
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
check_status("Proteowizard distribution" status ${error_message})

# Expand the tar file into the working directory
execute_process(
  COMMAND ${CMAKE_COMMAND} -E tar xf ${DOWNLOAD_DIR}/${filename}
  RESULT_VARIABLE status
  WORKING_DIRECTORY ${WORKING_DIR}
)
if (${status} EQUAL 0)
  message(STATUS "Generated ProteoWizard source files.")
else (${status} EQUAL 0)
  message(FATAL_ERROR "Unable to generate ProteoWizard source files.")
endif (${status} EQUAL 0)
