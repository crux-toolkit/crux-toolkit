#!/bin/bash -efx
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
echo set yrange \[0:3000\] >> $gnuplot
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

# Do the whole test twice, once for each search tool.  (N.B. This loop
# is temporarily not doing anything useful.  Once we get jimmy-search
# integrated into Crux, then this loop will be useful again.)
for searchtool in search-for-matches; do

  # Run the search.
  if [[ -e $searchtool/search.target.txt ]]; then
    echo Skipping search-for-matches.
  else  
    $CRUX $searchtool \
      --parameter-file crux.param \
      --num-decoys-per-target 1 \
      --output-dir $searchtool \
      $ms2 $db
  fi

  # Calibrate the scores.
  if [[ -e $searchtool/qvalues.target.txt ]]; then
    echo Skipping calibrate-scores \(Weibull\).
  else
    $CRUX calibrate-scores \
      --output-dir $searchtool \
      $db $searchtool
  fi

  $CRUX extract-columns $searchtool/qvalues.target.txt "Weibull est. q-value" > $searchtool/qvalues.weibull.txt
  echo replot \"$searchtool/qvalues.weibull.txt\" using 1:0 title \"XCorr \(Weibull\)\" with lines >> $gnuplot

  # Re-do the calibration, but without the Weibull p-values.
  if [[ -e decoys/qvalues.target.txt ]]; then
    echo Skipping calibrate-scores \(decoy\).
  else
    # Make copies of the search results without the p-value column.
    mkdir -p decoys
    cut -f 1-10,12-22 search-for-matches/search.target.txt \
      > decoys/search.target.txt
    cut -f 1-10,12-22 search-for-matches/search.decoy.txt \
      > decoys/search.decoy.txt
    $CRUX calibrate-scores \
      --output-dir decoys \
      $db decoys
  fi

  $CRUX extract-columns decoys/qvalues.target.txt "decoy q-value (xcorr)" > decoys/qvalues.xcorr.decoy.txt
  echo replot \"decoys/qvalues.xcorr.decoy.txt\" using 1:0 title \"XCorr \(decoy\)\" with lines >> $gnuplot
  
  # Run Crux percolator
  if [[ -e $searchtool/percolator.target.txt ]]; then
    echo Skipping crux percolator.
  else
    $CRUX percolator \
      --output-dir $searchtool \
      --feature-file T \
      $db $searchtool 
  fi

  $CRUX extract-columns $searchtool/percolator.target.txt "percolator q-value" > $searchtool/qvalues.percolator.txt
  echo replot \"$searchtool/qvalues.percolator.txt\" using 1:0 title \"percolator\" with lines >> $gnuplot

  # Run q-ranker.
  if [[ -e $searchtool/q-ranker.target.psms.txt ]]; then
    echo Skipping q-ranker.
  else
    $CRUX q-ranker \
      --output-dir $searchtool \
      --feature-file T \
      --separate-searches $searchtool/search.decoy.txt \
      $ms2 $searchtool/search.target.txt
  fi

  $CRUX extract-columns $searchtool/q-ranker.target.psms.txt "q-ranker q-value" > $searchtool/qvalues.qranker.txt
  echo replot \"$searchtool/qvalues.qranker.txt\" using 1:0 title \"q-ranker\" with lines >> $gnuplot

  # Run Barista.
  if [[ -e $searchtool/barista.target.psms.txt ]]; then
    echo Skipping barista.
  else
    $CRUX barista \
      --output-dir $searchtool \
      --feature-file T \
      --overwrite T \
      --separate-searches $searchtool/search.decoy.txt \
      $db $ms2 $searchtool/search.target.txt
  fi

  $CRUX extract-columns $searchtool/barista.target.psms.txt "q-value" > $searchtool/qvalues.barista.txt
  echo replot \"$searchtool/qvalues.barista.txt\" using 1:0 title \"barista\" with lines >> $gnuplot
  
done

# Finalize the gnuplot script.
echo set output >> $gnuplot
echo replot >> $gnuplot

# Make the plot.
gnuplot $gnuplot > performance.png
