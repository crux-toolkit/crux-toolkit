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
# Get the id of the latest run.
id_latest_run=$(gh run list -R github.com/CharlesEGrant/crux-toolkit -b master --workflow main.yml | \
  grep -oh "success.*"| head -n 1|awk '{print $7}')
echo $id_latest_run
# Download the artifacts for the latest run
gh run download $id_latest_run --pattern "crux*.*" -R github.com/CharlesEGrant/crux-toolkit 
# Copy artifacts to web site
cp crux-4.1.$new_version.Source/crux-4.1.Source.tar.gz /noble/www/htdocs/crux-downloads/daily/crux-4.1.$new_version.Source.tar.gz
cp crux-4.1.$new_version.windows/crux-4.1.Windows.i386.zip /noble/www/htdocs/crux-downloads/daily/crux-4.1.$new_version.Windows.AMD64.zip 
cp crux-4.1.$new_version.macos/crux-4.1.Darwin.x86_64.zip /noble/www/htdocs/crux-downloads/daily/crux-4.1.$new_version.Darwin.x86_64.zip
cp crux-4.1.$new_version.centos7/crux-4.1.Linux.x86_64.zip /noble/www/htdocs/crux-downloads/daily/crux-4.1.$new_version.Linux.x86_64.zip
cp latest-build.txt /noble/www/htdocs/crux-downloads/daily/latest-build.txt
rm -rf $temp_dir
else
  echo "No new build available."
fi
