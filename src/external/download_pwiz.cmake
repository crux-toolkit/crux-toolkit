# Download the build id
file(
  DOWNLOAD 
  "http://teamcity.labkey.org:8080/app/rest/buildTypes/id:bt81/builds?status=SUCCESS&count=1&guest=1"
  ./build-info.txt
  STATUS status
  LOG log.txt
)
message(${status})
file(STRINGS "build-info.txt" BUILD_INFO)
string(REGEX REPLACE "^.*build id=\"([0-9]+)\".*$" "\\1" BUILD_ID ${BUILD_INFO})
# Using the build number download the version string
file(
  DOWNLOAD 
  "http://teamcity.labkey.org:8080/repository/download/bt81/${BUILD_ID}:id/VERSION?guest=1"
  ./version-info.txt
  STATUS status
  LOG log.txt
)
message(${status})
set(
  BASE_DOWNLOAD_URL
  "http://teamcity.labkey.org:8080/guestAuth/repository/download/bt81/.lastSuccessful/"
)
file(STRINGS "version-info.txt" VERSION_INFO)
string(REGEX REPLACE "^([0-9]+)\\.([0-9]+)\\.([0-9]+).*$" "\\1_\\2_\\3" VERSION_ID ${VERSION_INFO})
set(
  DOWNLOAD_URL 
  "${BASE_DOWNLOAD_URL}pwiz-src-${VERSION_ID}.tar.bz2"
)
message(${DOWNLOAD_URL})
# Using the version string download the file
file(
  DOWNLOAD 
  ${DOWNLOAD_URL}
  ./pwiz-src-3_0_5471.tar.bz2
  STATUS status
  LOG log.txt
)
message(${status})
