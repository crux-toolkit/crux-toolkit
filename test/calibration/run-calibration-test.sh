#!/bin/bash -efx
# FILE: run-calibration-test.sh
# AUTHOR: William Stafford Noble

# This script runs a calibration test on the crux toolkit.  Its usage
# and output are described in ./calibration.html.  The MS2 file comes
# from the PAnDA paper by Hoopman et al. (JPR 2009).  Extensive
# documentation of this particular data set is available from
# crux-projects/panda-data.

# The location of the crux binary.
CRUX=../../src/crux

# Run the search.
$CRUX tide-search --top-match 1000000 \
       --overwrite T \
       --output-dir . \
       --exact-p-value T \
       ../performance-tests/051708-worm-ASMS-10.ms2 \
       ../performance-tests/worm+contaminants.fa

for peptide in target decoy; do

  # Extract all the p-values.
  $CRUX extract-columns --header F tide-search.$peptide.txt "exact p-value" \
      > $peptide.pvalues.txt

  # Make a histogram.
  histogram -bar-height distribution -minvalue 0 -binsize 0.01 100 \
      $peptide.pvalues.txt \
    | plot-histogram -xlabel "$peptide p-value" -format png - \
    > $peptide.histogram.png

  # Make a qq plot.
  ./make-qq-plot.py --title $peptide $peptide.pvalues.txt $peptide.qq

done
