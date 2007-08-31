#!/bin/bash
#------------
. ../SOURCEME
#------------

./match_analysis $FASTA --match-output-folder `pwd` --parameter-file $PARAM_FILE --algorithm percolator  > $OUTPUT
