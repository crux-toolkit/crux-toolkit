# run on gist
 (~/bin/percolator.n001  -g match_analysis.mtx match_analysis.labels > perc-gist.out ) >& perc-gist.err
# run on sqts
# (~/bin/percolator.n001  all-t.sqt all-d.sqt > perc-sqt.out) >& perc-sqt.err
# head -2000 match_analysis.features > 2000.features
#for i in `cat gist-header | cut -f2-`;{
#  echo $i
#  roc -column 0 -score $i 2000.features
#}
#roc -column 0 -score 1 -file xcorr.xy 2000.features

