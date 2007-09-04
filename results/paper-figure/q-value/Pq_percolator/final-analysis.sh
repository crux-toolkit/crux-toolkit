#!/bin/bash

#------------
. ../SOURCEME

#------------
export OUTPUT=`basename ${PREFIX}`.anal

#------------
../match_analysis $FASTA --match-output-folder `pwd` --parameter-file $PARAM_FILE --algorithm percolator  > $OUTPUT
cut -f1,5 $OUTPUT | cut -f2 -d' ' > pq.xy.tmp
cut -f2 pq.xy.tmp > 1.tmp
cut -f1 pq.xy.tmp > 2.tmp
paste 1.tmp 2.tmp > pq.xy
rm -f *.tmp

