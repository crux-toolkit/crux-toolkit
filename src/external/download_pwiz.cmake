# This is a CMake script for downloading the partial source distribution of 
# Proteowizard from its Team City repository. To get the URL for this
# distribution we have to find the build id, then use the build id to look up
# the version number. Once we have the version number we can build the URL for 
# downloading the distribution.

# Team City performs several different builds. We want the one
# for Pwiz source without Skyline but with vendor library suppport.

# This script expects to variables to be passed on the command line:
#   DOWNLOAD_DIR, the location to store the downloaded tar ball.
#   WORKING_DIR, the location where the expanded source files

set(build_type "bt81")

# This macro checks download status codes for errors
macro (check_status download_name status)
  list(GET status 0 status_code)
  list(GET status 1 status_message)
  if (${status_code} EQUAL 0)
    message(STATUS "Downloaded ${download_name}.")
  else (${status_code} EQUAL 0)
    message(
      FATAL_ERROR 
      "Unable to download ${download_name}: ${status_message} "
      "\nError code: ${status_code}." 
      "\nDownload of Proteowizard failed."
    )
  endif (${status_code} EQUAL 0)
endmacro (check_status)

# Get the build id
file(
  DOWNLOAD 
  "http://teamcity.labkey.org:8080/app/rest/buildTypes/id:${build_type}/builds?status=SUCCESS&count=1&guest=1"
  ${DOWNLOAD_DIR}/pwiz.build-info.txt
  STATUS status
)
check_status("build ID" status)

# Using the build id download the version string
file(STRINGS "${DOWNLOAD_DIR}/pwiz.build-info.txt" BUILD_INFO)
string(REGEX REPLACE "^.*build id=\"([0-9]+)\".*$" "\\1" BUILD_ID ${BUILD_INFO})
file(
  DOWNLOAD 
  "http://teamcity.labkey.org:8080/repository/download/${build_type}/${BUILD_ID}:id/VERSION?guest=1"
  ${DOWNLOAD_DIR}/pwiz.version-info.txt
  STATUS status
)
check_status("version string" status)

# Use the version string to build the URL
set(
  BASE_DOWNLOAD_URL
  "http://teamcity.labkey.org:8080/guestAuth/repository/download/${build_type}/.lastSuccessful/"
)
file(STRINGS "${DOWNLOAD_DIR}/pwiz.version-info.txt" VERSION_INFO)
string(REGEX REPLACE "^([0-9]+)\\.([0-9]+)\\.([0-9]+).*$" "\\1_\\2_\\3" VERSION_ID ${VERSION_INFO})
set(filename pwiz-src-${VERSION_ID}.tar.bz2)
set(
  DOWNLOAD_URL 
  "${BASE_DOWNLOAD_URL}${filename}"
)
# Using the version string download the file
file(
  DOWNLOAD 
  ${DOWNLOAD_URL}
  "${DOWNLOAD_DIR}/${filename}"
  STATUS status
)
check_status("Proteowizard distribution" status)

# Expand the tar file into the working directory
execute_process(
  COMMAND ${CMAKE_COMMAND} -E tar xf ${DOWNLOAD_DIR}/${filename}
  WORKING_DIRECTORY ${WORKING_DIR}
  RESULT_VARIABLE rv
)
if (${rv} EQUAL 0)
  message(STATUS "Generated ProteoWizard source files.")
else (${rv} EQUAL 0)
  message(FATAL_ERROR "Unable to generate ProteoWizard source files.")
endif (${rv} EQUAL 0)
