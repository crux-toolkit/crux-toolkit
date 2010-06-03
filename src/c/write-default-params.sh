#!/bin/bash

# A script to create the default.params file to distribute in the docs
# and to store in src/c.  Writes a param file with default values.
# Adds a header.  Copies to documents directory.  This script assumes
# that it is being run from crux/src/c.  It is up to the developer to
# run this script after any changes are made to the parameters.

# Create param file
rm -rf crux-output 
./crux search-for-matches smoke/test.ms2 smoke/test.fasta 1>/dev/null 2>&1

if [ ! -e crux-output/search.params.txt ]
then
  echo "Failed to write the param file.";
  exit;
fi

# Write a header
echo "####################################################################
# Sample parameter file
#
# Lines starting with '#' will be ignored.  Don't leave any space
# before or after a parameter setting.  The format is
#
# <parameter-name>=<value>
#
#####################################################################
" > default.params
cat crux-output/search.params.txt >> default.params

# Rename the value for parameter-file, and copy to doc directory.
# N.B. Don't use the "-i" option to 'sed' because it has different
# syntax under Linux and Darwin.  Instead, combine the sed and the copy,
# then copy the sed output back to the local directory.
sed 's/parameter-file=T/parameter-file=F/' default.params \
  > ../../doc/user/default.params
cp ../../doc/user/default.params default.params

# Remove output files
rm -rf crux-output
