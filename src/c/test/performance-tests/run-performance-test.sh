#!/bin/bash -efx
# FILE: run-performance-test.sh
# AUTHOR: William Stafford Noble

# This script runs a performance test on the crux toolkit.  Its usage
# and output are described in ./performance.html.  The MS2 file comes
# from the PAnDA paper by Hoopman et al. (JPR 2009).  Extensive
# documentation of this particular data set is available from
# crux-projects/panda-data.

# Increase the file limit for Crux. (Necessary on MacOS.)
ulimit -n 1024

# Initialize the gnuplot script.
gnuplot=performance.gnuplot
echo set output \"/dev/null\" > $gnuplot
echo set terminal png >> $gnuplot
echo set xlabel \"q-value threshold\" >> $gnuplot
echo set ylabel \"Number of accepted PSMs\" >> $gnuplot
echo set xrange \[0:0.1\] >> $gnuplot
echo set key outside >> $gnuplot
echo set size 1.5, 1
# Insert dummy plot so we can use "replot" consistently below.
echo plot 0 notitle with dots >> $gnuplot 

# Create the index.
db=worm+contaminants
if [[ -e $db ]]; then
  echo Skipping create-index.
else
  crux create-index $db.fa $db
fi

ms2=051708-worm-ASMS-10.ms2

# Do the whole test twice, once for each search tool.
for searchtool in sequest-search search-for-matches; do

  if [[ $searchtool == "sequest-search" ]]; then
     shortname=sequest
     search_parameter=""
  else
     shortname=search
     search_parameter="--compute-p-values T"
  fi

  # Run the search.
  if [[ -e $shortname/$shortname.target.txt ]]; then
    echo Skipping $searchtool.
  else  
    crux $searchtool \
      $search_parameter \
      --num-decoys-per-target 1 \
      --output-dir $shortname \
      $ms2 $db
  fi

  # Run Lukas's percolator
  if [[ $searchtool == "search-for-matches" ]]; then
    echo Stand-alone Percolator does not work with crux search-for-matches.
  elif [[ `which percolator` == "" ]]; then
    echo Skipping stand-alone Percolator -- not installed.
  else
    if [[ -e $shortname/l-percolator.features.tsv ]]; then
      echo Skipping stand-alone Percolator.
    else
      percolator \
        --tab-out $shortname/l-percolator.features.tsv \
        --gist-out $shortname/l-percolator.gist.txt \
        --weights $shortname/l-percolator.weights.txt \
        --results $shortname/l-percolator.tsv \
        --sqt-out $shortname/l-percolator.sqt \
        --xml-output $shortname/l-percolator.xml \
        $shortname/sequest.target.sqt \
        $shortname/sequest.decoy-1.sqt
    fi
    echo replot \"$shortname/l-percolator.tsv\" using 3:0 title \"Stand-alone percolator\" with lines >> $gnuplot
  fi
  
  # Run compute-q-values.
  if [[ -e $shortname/qvalues.target.txt ]]; then
    echo Skipping compute-q-values.
  else
    crux compute-q-values \
      --output-dir $shortname \
      $db
  fi
  echo replot \"$shortname/qvalues.target.txt\" using 9:0 title \"$shortname XCorr \(decoy\)\" with lines >> $gnuplot
  if [[ $searchtool == "search-for-matches" ]]; then
    echo replot \"$shortname/qvalues.target.txt\" using 9:0 title \"$shortname XCorr \(Weibull\)\" with lines >> $gnuplot
    echo replot \"$shortname/qvalues.target.txt\" using 9:0 title \"$shortname XCorr \(decoy p-value\)\" with lines >> $gnuplot
  fi
  
  # Run Crux percolator
  if [[ -e $shortname/percolator.target.txt ]]; then
    echo Skipping crux percolator.
  else
    crux percolator \
      --output-dir $shortname \
      --feature-file T \
      $db
  fi
  echo replot \"$shortname/percolator.target.txt\" using 13:0 title \"$shortname crux percolator\" with lines >> $gnuplot
  
  # Run q-ranker.
  if [[ -e $shortname/qranker.target.txt ]]; then
    echo Skipping q-ranker.
  else
    crux q-ranker \
      --output-dir $shortname \
      --feature-file T \
      $db
  fi
  echo plot \"$shortname/qranker.target.txt\" using 12:0 title \"$shortname q-ranker\" with lines >> $gnuplot
  
done

# Finalize the gnuplot script.
echo set output >> $gnuplot
echo replot >> $gnuplot

# Make the plot.
gnuplot $gnuplot > performance.png
