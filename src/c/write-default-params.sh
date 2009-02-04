#!/bin/bash

# A script to create the default.params file to distribute in the docs
# and to store in src/c.
# Writes a param file with default values.  Re-orders the options and
# adds a header.  Copies to documents directory.

# Assumes it is being run from crux/src/c

# 1. Create param file
rm -f unordered-default-params test*csm
./crux search-for-matches \
  --write-parameter-file unordered-default-params \
  smoke/test.ms2 smoke/test.fasta

if [ ! -e unordered-default-params ]
then
  echo "Failed to write the param file.";
  exit;
fi

# 2. Write a header
echo "####################################################################
# Sample parameter file
#
# Lines starting with '#' will be ignored
# Don't leave any space before or after a parameter setting
# format: <parameter-name>=<value>
#
#####################################################################
" > default.params

# 3. Reorder the options so they are grouped by function
for op in verbosity version parameter-file write-parameter-file \
          overwrite use-index min-length max-length isotopic-mass \
          fragment-mass mass-window ion-tolerance cleavages \
          missed-cleavages unique-peptides mod cmod nmod max-mods \
          max-aas-modified \
          max-rank-preliminary min-mass max-mass spectrum-min-mass \
          spectrum-max-mass spectrum-charge output-mode \
          match-output-folder sqt-output-file decoy-sqt-output-file \
          number-decoy-set top-match top-match-sqt algorithm \
          feature-file output-trypticity output-sequence sort stats \
          isotope primary-ions neutral-losses flanking \
          precursor-ions nh3 max-ion-charge h2o A C D E F G H I K L M \
          N P Q R S T V W Y ; 
do
  grep -B2 "^$op=" unordered-default-params >> default.params
  echo "" >> default.params
done

# 4. Rename the value for write-parameter-file
sed -i 's/unordered-default-params/__NULL_STR/' default.params

# 5. Copy to crux/doc/user
cp -f default.params ../../doc/user/default.params
