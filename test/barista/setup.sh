#!/bin/bash -x
# This script moves all of the necessary files to $4
# unzips them, runs tide index and tide search
# It is the common code between *_job.sh

if [[ -e /etc/profile.d/modules.sh ]]; then
  source /etc/profile.d/modules.sh
  module load modules modules-init modules-gs modules-noble
fi

set -o nounset
set -o errexit 
set -o xtrace

# $1 is fasta.gz file
# $2 dir of ms2.gz
# $3 number of ms2 to copy
# $4 the directory to which to copy the files and unzip them
# $5 the desired name of the concatenated database
FASTA_ZIP=$1
MS2ZIP=$2
NUM_MS2=$3
TO_DIR=$4
CONCAT=$5

# so long as we are consistent with the ordering, just keeping track of the
# number of files to compute is sufficient
MS2="${TO_DIR}/MS2"
mkdir $MS2
find $MS2ZIP -maxdepth 1 -type f |head -$NUM_MS2|xargs cp -t $MS2

cp $FASTA_ZIP -t $TO_DIR

FASTA_ZIP=$(basename $FASTA_ZIP)
gunzip "${TO_DIR}/$FASTA_ZIP"
gunzip $MS2/*

# N.B. Be sure to copy everything here before running this script.
# cp ~/proj/crux/trunk/release/src/crux .

CRUX=./crux
# TODO: verify that this (running tide search / index) works correctly

if [[ ! -e Linfeng-index ]]; then
  $CRUX tide-index \
    --clip-nterm-methionine T \
    --missed-cleavages 2 \
    --peptide-list T \
    --mods-spec C+57.0214,K+229.16293 \
    --nterm-peptide-mods-spec X+229.16293 \
    --max-mods 0 \
    --enzyme trypsin/p \
    --output-dir Linfeng-index \
    --overwrite T \
    $FASTA Linfeng-index
fi
# Search using exact p-values and low-res.
results=exact-p/tide-search.target.txt
if [[ ! -e $results ]]; then
  $CRUX tide-search \
    --precursor-window 10 \
    --precursor-window-type ppm \
    --top-match 1 \
    --exact-p-value T \
    --output-dir exact-p \
    --overwrite T \
    "$MS2/*" Linfeng-index
fi

$DECOY_INDEX="Linfeng-index/tide-index.decoy.fasta"
$(cat $DECOY_INDEX $FASTA > $CONCAT)
