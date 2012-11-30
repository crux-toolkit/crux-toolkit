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
echo num-decoys-per-target=0 >> $params
echo top-match=1 >> $params
echo output-dir=search >> $params
echo compute-p-values=T >> $params
#echo decoy-p-values=T >> $params  # Write raw p-values to a separate file.
echo decoys=peptide-shuffle >> $params

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

# Extract just the p-value column.
$CRUX extract-columns search/search.target.txt p-value \
  | awk 'NR > 1' \
  > pvalues.txt

# Make a histogram.
histogram -bar-height distribution -minvalue 0 -binsize 0.01 100 pvalues.txt \
  | plot-histogram -format png - > histogram.png

# Make a qq plot.
./make-qq-plot.py pvalues.txt qq
