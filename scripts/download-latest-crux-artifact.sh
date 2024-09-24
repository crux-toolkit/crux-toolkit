#!/usr/bin/bash
set -eux pipefail
# Authenticate to GitHub
gh auth login --with-token 
# Set up temp working directory
temp_dir=$(mktemp -d  gh-crux-download.XXX)
cd $temp_dir
# Download short git commit hash for latest sucessful build
wget 'https://crux-toolkit.github.io/latest-build.txt' -O 'latest-build.txt'
new_version=$(cat latest-build.txt)
old_version=$(cat /noble/www/htdocs/crux-downloads/daily/latest-build.txt)
echo "old version is " $old_version
echo "new version is " $new_version
if [ "$new_version" != "$old_version" ];
then
echo "$(date): New version detected $new_version" >> /noble/www/htdocs/crux-downloads/daily/crux-update.log
# Get the id of the latest run.
id_latest_run=$(gh run list -R github.com/crux-toolkit/crux-toolkit --workflow main.yml | \
  grep -oh "success.*"| head -n 1|awk -F '\t' '{print $6}')
echo $id_latest_run
# Download the artifacts for the latest run
gh run download $id_latest_run  -R github.com/crux-toolkit/crux-toolkit
# Copy artifacts to web site
cp crux-4.2.$new_version.Source/crux-4.2.Source.tar.gz /noble/www/htdocs/crux-downloads/daily/crux-4.2.$new_version.Source.tar.gz
echo "$(date): Updated crux-4.2.$new_version.Source.tar.gz" >> /noble/www/htdocs/crux-downloads/daily/crux-update.log
cp crux-4.2.$new_version.windows/crux-4.2.Windows.AMD64.zip /noble/www/htdocs/crux-downloads/daily/crux-4.2.$new_version.Windows.AMD64.zip 
echo "$(date): Updated crux-4.2.$new_version.Windows.AMD64.zip" >> /noble/www/htdocs/crux-downloads/daily/crux-update.log
cp crux-4.2.$new_version.macos/crux-4.2.Darwin.x86_64.zip /noble/www/htdocs/crux-downloads/daily/crux-4.2.$new_version.Darwin.x86_64.zip
echo "$(date): Updated crux-4.2.$new_version.Darwin.x86_64.zip" >> /noble/www/htdocs/crux-downloads/daily/crux-update.log
cp crux-4.2.$new_version.ubuntu.20.04/crux-4.2.Linux.x86_64.zip /noble/www/htdocs/crux-downloads/daily/crux-4.2.$new_version.Linux.x86_64.zip
echo "$(date): Updated crux-4.2.$new_version.Linux.x86_64.zip" >> /noble/www/htdocs/crux-downloads/daily/crux-update.log
cp latest-build.txt /noble/www/htdocs/crux-downloads/daily/latest-build.txt
echo "$(date): Updated latest-build.txt" >> /noble/www/htdocs/crux-downloads/daily/crux-update.log
else
  echo "$(date): No new build available." >> /noble/www/htdocs/crux-downloads/daily/crux-update.log

fi
cd -
rm -r $temp_dir
