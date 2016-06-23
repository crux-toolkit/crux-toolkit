#!/bin/bash -x

# The main executable for uge memory jobs. Calls setup.sh for commone code and
# then analyzes a run of barista on the search results from setup.sh

if [[ -e /etc/profile.d/modules.sh ]]; then
  source /etc/profile.d/modules.sh
  module load modules modules-init modules-gs modules-noble
fi

module load valgrind/latest

set -o nounset
set -o errexit 
set -o xtrace

FASTA_ZIP=$1
MS2_ZIP=$2
NUM_FILES=$3
CONCAT="concat.txt"

./setup.sh $FASTA_ZIP $MS2_ZIP $NUM_MS2 $TMPDIR $CONCAT

$SEARCH_RESULT="$TMPDIR/exact-p/tide-search.target.txt"
CRUX=./crux

outdir="$HOME/barista_benchmark/"
if [[ ! -e $outdir ]]; then
    mkdir $outdir
fi
root="$outdir/$NUM_FILES.memory"
usage=$root.usage.txt
# TODO: massif may not be the most appropriate tool
# We only need the max memory usage and the output is difficult to parse
valgrind --tool=massif --depth=1 --massif-out-file=$root.massif \
	 $CRUX barista --output-dir $TMPDIR/exact-p $CONCAT \
	 $TMPDIR/MS2 $SEARCH_RESULT
ms_print $root.massif > $root.massif.txt
~wnoble/bin/get-max-memory-usage.py $root.massif.txt > $usage
~wnoble/bin/getsize-ms2.py $MS2/* > $root.size.txt
