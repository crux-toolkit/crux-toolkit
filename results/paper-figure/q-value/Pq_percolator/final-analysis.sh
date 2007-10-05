#!/bin/bash

#------------
. ../SOURCEME

#------------
export OUTPUT=`basename ${PREFIX}`.anal

#------------
# (../match_analysis `pwd` $FASTA --feature-file match_analysis.features --parameter-file $PARAM_FILE --algorithm percolator --verbosity 100 > $OUTPUT  )>& error
cut -f1,5 $OUTPUT | cut -f2 -d' ' > pq.xy.tmp
cut -f2 pq.xy.tmp > 1.tmp
cut -f1 pq.xy.tmp > 2.tmp
paste 1.tmp 2.tmp > pq.xy

cut -f1 match_analysis.features > labels.tmp

cut -f2- match_analysis.features > features.tmp
cut -f1 $OUTPUT | sed 's/ /_/' | sed 's/P/T/' > strings.tmp
cut -f1 $OUTPUT | sed 's/ /_/' | sed 's/P/D/' >> strings.tmp
head -2000 features.tmp > features-head.tmp
paste strings.tmp features-head.tmp > a.tmp
paste strings.tmp labels.tmp | head -2000 | sed 's/	0/	-1/' > b.tmp
cat gist-header a.tmp > match_analysis.mtx
cat gist-header b.tmp | cut -f1,2 > match_analysis.labels

rm -f *.tmp

