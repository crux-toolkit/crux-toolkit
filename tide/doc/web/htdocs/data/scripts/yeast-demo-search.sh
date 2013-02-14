#!/bin/bash
if [ ! -e yeast.fasta.pepix -o ! -e yeast.fasta.protix ]
 then
  echo "Cannot find files created in indexing step. Please run yeast-demo-index.sh first." 1>&2
  exit 1
fi
time ./tide-search --peptides=yeast.fasta.pepix --proteins=yeast.fasta.protix \
              --spectra=yeast-02-10000.spectrumrecords  > yeast.results
