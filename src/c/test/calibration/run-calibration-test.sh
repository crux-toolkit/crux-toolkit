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

if [[ -e $db ]]; then
  echo Skipping create-index.
else
  $CRUX create-index $db.fa $db
fi

ms2=../performance-tests/051708-worm-ASMS-10.ms2

# Do the whole test twice, once for each search tool.
searchtool=search-for-matches

# Run the search.
if [[ -e search/search.target.txt ]]; then
  echo Skipping search-for-matches.
else  
  $CRUX search-for-matches \
    --compute-p-values T \
    --num-decoys-per-target 1 \
    --output-dir search \
    $ms2 $db
fi

# Run compute-q-values.
if [[ -e search/qvalues.target.txt ]]; then
  echo Skipping compute-q-values.
else
  $CRUX compute-q-values \
    --output-dir search \
    $db search
fi

# Run Crux percolator
if [[ -e search/percolator.target.txt ]]; then
  echo Skipping crux percolator.
else
  $CRUX percolator \
    --output-dir search \
    --feature-file T \
    $db search 
fi

# Run q-ranker.
if [[ -e search/qranker.target.txt ]]; then
  echo Skipping q-ranker.
else
  $CRUX q-ranker \
    --output-dir search \
    --feature-file T \
    $db search
fi

