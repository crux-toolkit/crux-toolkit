#!/bin/bash
if [ ! -e worm.fasta.pepix -o ! -e worm.fasta.protix ]
 then
  echo "Cannot find files created in indexing step. Please run worm-demo-index.sh first." 1>&2
  exit 1
fi
time ./tide-search --peptides=worm.fasta.pepix --proteins=worm.fasta.protix \
              --spectra=worm-06-10000.spectrumrecords  > worm.results
