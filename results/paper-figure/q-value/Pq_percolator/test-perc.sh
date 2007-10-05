# (~/bin/percolator.n001  -g match_analysis.mtx match_analysis.labels > perc-gist.out ) >& perc-gist.err
# (~/bin/percolator.n001  all-t.sqt all-d.sqt > perc-sqt.out) >& perc-sqt.err
for i in `cat gist-header | cut -f2-`;{
  echo $i
  roc -column 0 -score $i match_analysis.features
}

