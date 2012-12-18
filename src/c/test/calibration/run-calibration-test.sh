#!/bin/bash -efx
# FILE: run-calibration-test.sh
# AUTHOR: William Stafford Noble

# This script runs a calibration test on the crux toolkit.  Its usage
# and output are described in ./calibration.html.  The MS2 file comes
# from the PAnDA paper by Hoopman et al. (JPR 2009).  Extensive
# documentation of this particular data set is available from
# crux-projects/panda-data.

# The database is a shuffled version of
# ../performance-tests/worm+contaminants.fa
# 
# It was created with this command line:
# fasta-shuffle-letters < ../performance-tests/worm+contaminants.fa > shuffle.fa
#
db=shuffle

# The location of the crux binary.
CRUX=../../crux

# Increase the file limit for Crux. (Necessary on MacOS.)
ulimit -n 1024

# Create a parameter file.
params=parameters.txt
echo num-decoys-per-target=1 > $params
echo top-match=1 >> $params
echo output-dir=search >> $params
echo compute-p-values=T >> $params
echo decoy-p-values=T >> $params  # Write raw p-values to a separate file.
echo decoys=peptide-shuffle >> $params
echo precursor-window=3 >> $params

if [[ -e $db ]]; then
  echo Skipping create-index.
else
  $CRUX create-index --parameter-file $params $db.fa $db
fi

ms2=../performance-tests/051708-worm-ASMS-10.ms2


# Run the search.
if [[ -e search/search.target.txt ]]; then
  echo Skipping search-for-matches.
else  
  $CRUX search-for-matches --parameter-file $params $ms2 $db
fi

for peptide in target decoy; do

  # Extract just the top-ranked p-values from the p-value column.
  awk 'NR == 1 || $8 == 1' search/search.$peptide.txt \
    | $CRUX extract-columns - p-value \
    | awk 'NR > 1' \
    > $peptide.pvalues.txt

  # Make a histogram.
  histogram -bar-height distribution -minvalue 0 -binsize 0.01 100 \
      $peptide.pvalues.txt \
    | plot-histogram -xlabel "$peptide p-value" -format png - \
    > $peptide.histogram.png

  # Make a qq plot.
  ./make-qq-plot.py $peptide.pvalues.txt $peptide.qq
done

# Make a histogram of uncorrected decoy p-values.
awk '$1 != "#"' search/search.decoy.p.txt \
  | histogram -bar-height distribution -minvalue 0 -binsize 0.01 100 - \
  | plot-histogram -format png - > decoy.png

