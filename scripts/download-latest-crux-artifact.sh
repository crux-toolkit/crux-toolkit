#!/usr/bin/bash
set -eux pipefail
# Authenticate to GitHub
gh auth login --with-token < access-repo.txt
# Set up temp working directory
temp_dir=$(mktemp -d  gh-crux-download.XXX)
cd $temp_dir
# Download short git commit hash for latest sucessful build
wget 'https://charlesegrant.github.io/crux-toolkit.github.io/latest-build.txt' -O 'latest-build.txt'
new_version=$(cat latest-build.txt)
old_version=$(cat /noble/www/htdocs/crux-downloads/daily/latest-build.txt)
echo "old version is " $old_version
echo "new version is " $new_version
if [ "$new_version" != "$old_version" ];
then
  echo versions were different
fi
# Get the id of the latest run.
id_latest_run=$(gh run list -R github.com/CharlesEGrant/crux-toolkit -b master --workflow main.yml | \
  grep -oh "success.*"| head -n 1|awk '{print $7}')
echo $id_latest_run
# Download the artifacts for the latest run
gh run download $id_latest_run --pattern "crux*.*" -R github.com/CharlesEGrant/crux-toolkit 
