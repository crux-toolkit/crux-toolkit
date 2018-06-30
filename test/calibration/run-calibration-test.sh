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

# Note that the scripts make-qq-plot.py, histogram and plot-histogram
# are taken from the 2002_wnoble_utilities bitbucket repository.

index=worm+contaminants
if [[ ! -e $index ]]; then
  $CRUX tide-index \
      --output-dir $index \
      ../performance-tests/worm+contaminants.fa \
      $index
fi

for score in xcorr residue-evidence both; do

  # Run the search.
  results=$score/tide-search.target.txt
  if [[ ! -e $results ]]; then
    $CRUX tide-search --top-match 1000000 \
       --use-neutral-loss-peaks F \
       --score-function $score \
       --output-dir $score \
       --exact-p-value T \
       --num-threads 8 \
       ../performance-tests/051708-worm-ASMS-10.ms2 \
       $index
  fi

  for peptide in target decoy; do

    if [[ $score == "xcorr" ]]; then
      column="exact p-value"
    elif [[ $score == "residue-evidence" ]]; then
      column="res-ev p-value"
    else
      column="combined p-value"
    fi
      
    # Extract all the p-values.
    pvalues=$score/$peptide.pvalues.txt
    if [[ ! -e $pvalues ]]; then
      $CRUX extract-columns --header F \
            $score/tide-search.$peptide.txt "$column" \
        > $pvalues
    fi

    # Make a histogram.
    ./histogram -bar-height distribution -minvalue 0 -binsize 0.01 100 $pvalues \
      | ./plot-histogram -xlabel "$peptide p-value" -format png - \
      > $score/$peptide.histogram.png

    # Make a qq plot.
    ./make-qq-plot.py --title "$score $peptide" $score/$peptide.pvalues.txt $score/$peptide.qq

  done
done
