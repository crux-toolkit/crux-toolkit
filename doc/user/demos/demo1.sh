#!/bin/bash -x
# AUTHOR: WSN
# CREATE DATE: 18 December 2013, Tokyo Narita airport

# This file implements the commands described in ../search-tutorial.html.

CRUX=../../../src/c/crux

rm -r yeast-index crux-output
cp ../data/demo.ms2 .
cp ../data/small-yeast.fasta .
$CRUX tide-index small-yeast.fasta yeast-index
$CRUX tide-search --compute-sp T demo.ms2 yeast-index
$CRUX percolator crux-output/tide-search.target.txt 

rdb2html -noformatline crux-output/percolator.target.psms.txt \
  > crux-output/percolator.target.psms.html
