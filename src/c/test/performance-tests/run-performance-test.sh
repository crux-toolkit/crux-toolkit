#!/bin/bash -efx
# FILE: run-performance-test.sh
# AUTHOR: William Stafford Noble

# This script runs a performance test on the crux toolkit.  Its usage
# and output are described in ./performance.html.  The MS2 file comes
# from the PAnDA paper by Hoopman et al. (JPR 2009).  Extensive
# documentation of this particular data set is available from
# crux-projects/panda-data.

# The location of the crux binary.
CRUX=../../crux

# Increase the file limit for Crux. (Necessary on MacOS.)
ulimit -n 1024

# Initialize the gnuplot script.
gnuplot=performance.gnuplot
echo set output \"/dev/null\" > $gnuplot
echo set terminal png >> $gnuplot
echo set xlabel \"q-value threshold\" >> $gnuplot
echo set ylabel \"Number of accepted PSMs\" >> $gnuplot
echo set xrange \[0:0.1\] >> $gnuplot
echo set key center right >> $gnuplot
# Insert dummy plot so we can use "replot" consistently below.
echo plot 0 notitle with dots >> $gnuplot 

# Create the index.
db=worm+contaminants
if [[ -e $db ]]; then
  echo Skipping create-index.
else
  $CRUX create-index --parameter-file crux.param $db.fa $db
fi

ms2=051708-worm-ASMS-10.ms2

# Do the whole test twice, once for each search tool.
for searchtool in sequest-search search-for-matches; do

  if [[ $searchtool == "sequest-search" ]]; then
     shortname=sequest
     search_parameter="--parameter-file crux.param"
  else
     shortname=search
     search_parameter="--parameter-file crux.param --compute-p-values T --compute-sp T"
  fi

  # Run the search.
  if [[ -e $shortname/$shortname.target.txt ]]; then
    echo Skipping $searchtool.
  else  
    $CRUX $searchtool \
      $search_parameter \
      --num-decoys-per-target 1 \
      --output-dir $shortname \
      $ms2 $db
  fi

  # Run compute-q-values.
  if [[ -e $shortname/qvalues.target.txt ]]; then
    echo Skipping compute-q-values.
  else
    $CRUX compute-q-values \
      --output-dir $shortname \
      $db $shortname
  fi
  if [[ $searchtool == "sequest-search" ]]; then
    $CRUX extract-columns $shortname/qvalues.target.txt "decoy q-value (xcorr)" > $shortname/qvalues.xcorr.decoy.txt
    echo replot \"$shortname/qvalues.xcorr.decoy.txt\" using 1:0 title \"$shortname XCorr \(decoy\)\" with lines >> $gnuplot
  else  
    $CRUX extract-columns $shortname/qvalues.target.txt "Weibull est. q-value" > $shortname/qvalues.weibull.txt
    echo replot \"$shortname/qvalues.weibull.txt\" using 1:0 title \"$shortname XCorr \(Weibull\)\" with lines >> $gnuplot
  fi
  
  # Run Crux percolator
  if [[ -e $shortname/percolator.target.txt ]]; then
    echo Skipping crux percolator.
  else
    $CRUX percolator \
      --output-dir $shortname \
      --feature-file T \
      $db $shortname 
  fi

  $CRUX extract-columns $shortname/percolator.target.txt "percolator q-value" > $shortname/qvalues.percolator.txt
  echo replot \"$shortname/qvalues.percolator.txt\" using 1:0 title \"$shortname crux percolator\" with lines >> $gnuplot

  # Run q-ranker.
  if [[ -e $shortname/qranker.target.txt ]]; then
    echo Skipping q-ranker.
  else
    $CRUX q-ranker \
      --output-dir $shortname \
      --feature-file T \
      $db $shortname
  fi
  $CRUX extract-columns $shortname/qranker.target.txt "q-ranker q-value" > $shortname/qvalues.qranker.txt

  
  echo replot \"$shortname/qvalues.qranker.txt\" using 1:0 title \"$shortname q-ranker\" with lines >> $gnuplot
  
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
        $shortname/sequest.decoy.sqt
    fi
    echo replot \"$shortname/l-percolator.tsv\" using 3:0 title \"Stand-alone percolator\" with lines >> $gnuplot
  fi
  
done

# Finalize the gnuplot script.
echo set output >> $gnuplot
echo replot >> $gnuplot

# Make the plot.
gnuplot $gnuplot > performance.png
