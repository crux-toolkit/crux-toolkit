#!/bin/bash -x
# AUTHOR: WSN
# CREATE DATE: 18 December 2013, Tokyo Narita airport

# This file implements the commands described in ../search-tutorial.html.

if [[ -e ../../../src/crux ]]; then
  CRUX=../../../src/crux
else
  CRUX=../../../release/src/crux
fi    
    
rm -rf yeast-index crux-output
chmod +x ./tsv2html.py
cp ../example-files/demo.ms2 .
cp ../example-files/small-yeast.fasta .

$CRUX tide-index small-yeast.fasta yeast-index

$CRUX tide-search --compute-sp T demo.ms2 yeast-index
$CRUX sort-by-column --column-type real --ascending F \
      crux-output/tide-search.target.txt "xcorr score" \
      > crux-output/tide-search.target.sort.txt
head -11 crux-output/tide-search.target.sort.txt \
     | ./tsv2html.py - \
     > crux-output/tide-search.target.sort.html

$CRUX percolator crux-output/tide-search.target.txt 
$CRUX sort-by-column --column-type real --ascending F \
      crux-output/percolator.target.psms.txt "percolator score" \
      > crux-output/percolator.target.psms.sort.txt
head -11 crux-output/percolator.target.psms.sort.txt \
     | ./tsv2html.py - \
     > crux-output/percolator.target.psms.sort.html
