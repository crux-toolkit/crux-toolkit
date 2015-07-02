#!/bin/bash

# This script downloads the source tarball for the most recent successful build of
# ProteoWizard including dependencies (build type = bt81)

set -e
set -o pipefail

build_id=$(wget --no-check-certificate -q -O -\
  https://teamcity.labkey.org:/app/rest/buildTypes/id:bt81/builds\?status=SUCCESS\&count=1\&guest=1)
build_id=$(echo $build_id| sed 's/^.*build id=.\([0-9]*\).*/\1/')
echo "build_id=" $build_id
version=$(wget --no-check-certificate -q -O -\
  https://teamcity.labkey.org/repository/download/bt81/$build_id:id/VERSION?guest=1)
version=$(echo $version|sed 's/\./_/g')
echo "version = " $version
wget --no-check-certificate \
  https://teamcity.labkey.org:/guestAuth/repository/download/bt81/.lastSuccessful/pwiz-src-$version.tar.bz2
