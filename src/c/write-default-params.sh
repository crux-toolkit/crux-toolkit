#!/bin/bash

# A script to create the default.params file to distribute in the docs
# and to store in src/c.
# Writes a param file with default values.  Re-orders the options and
# adds a header.  Copies to documents directory.

# Assumes it is being run from crux/src/c

# 1. Create param file
rm -rf crux-output 
./crux search-for-matches smoke/test.ms2 smoke/test.fasta 1>/dev/null 2>&1

if [ ! -e crux-output/search.params.txt ]
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
for op in verbosity version parameter-file overwrite \
          output-dir fileroot print-search-progress \
          decoy-location num-decoys-per-target reverse-sequence \
          top-match precision search-decoy-pvalue-file \
          min-length max-length isotopic-mass fragment-mass \
          mass-window use-mz-window ion-tolerance \
          min-mass max-mass spectrum-min-mass spectrum-max-mass \
          spectrum-charge max-rank-preliminary \
          enzyme custom-enzyme digestion missed-cleavages \
          mod cmod nmod max-mods max-aas-modified compute-p-values \
          feature-file output-sequence sort stats unique-peptides \
          isotope primary-ions neutral-losses flanking \
          precursor-ions nh3 max-ion-charge h2o A C D E F G H I K L M \
          N P Q R S T V W Y ; 



do
  grep -B2 "^$op=" crux-output/search.params.txt >> default.params
  echo "" >> default.params
done

# 4. Rename the value for write-parameter-file
sed -i 's/parameter-file=T/parameter-file=F/' default.params

# 5. Copy to crux/doc/user
cp -f default.params ../../doc/user/default.params
