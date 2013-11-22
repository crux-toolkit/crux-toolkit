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
fasta=$db.fa
if [[ -e $db ]]; then
  echo Skipping create-index.
else
  $CRUX tide-index --output-dir tide-index --decoy-format reverse $fasta $db
fi

ms2=051708-worm-ASMS-10.ms2

# Do the whole test twice, once for each search tool.
for searchtool in comet tide-search; do

  # Do we use an index or the fasta?
  if [[ $searchtool == "comet" ]]; then
    params="--parameter-file crux.param"
    proteins=$fasta
  else
    params=""
    proteins=$db
  fi

  # Run the search.
  if [[ -e $searchtool/$searchtool.target.txt ]]; then
    echo Skipping search-for-matches.
  else  
    $CRUX $searchtool \
      $params --output-dir $searchtool \
      $ms2 $proteins
  fi

  # Run calibrate-scores
  if [[ -e $searchtool/qvalues.target.txt ]]; then
    echo Skipping crux calibrate-scores.
  else
    $CRUX calibrate-scores \
      --output-dir $searchtool \
      $searchtool/$searchtool.target.txt
  fi

  $CRUX extract-columns $searchtool/qvalues.target.txt "decoy q-value (xcorr)" > $searchtool/qvalues.xcorr.txt
  echo replot \"$searchtool/qvalues.xcorr.txt\" using 1:0 title \"$searchtool decoy \(xcorr\)\" with lines >> $gnuplot

  # Run Crux percolator
  if [[ -e $searchtool/percolator.target.psms.txt ]]; then
    echo Skipping crux percolator.
  else
    $CRUX percolator \
      --output-dir $searchtool \
      --feature-file T \
      $searchtool/$searchtool.target.txt
  fi

  $CRUX extract-columns $searchtool/percolator.target.psms.txt "percolator q-value" > $searchtool/qvalues.percolator.txt
  echo replot \"$searchtool/qvalues.percolator.txt\" using 1:0 title \"$searchtool percolator\" with lines >> $gnuplot

  # Run q-ranker.
  if [[ -e $searchtool/q-ranker.target.psms.txt ]]; then
    echo Skipping q-ranker.
  else
    $CRUX q-ranker \
      --output-dir $searchtool \
      --feature-file T \
      --separate-searches $searchtool/$searchtool.decoy.txt \
      $ms2 $searchtool/$searchtool.target.txt
  fi

  $CRUX extract-columns $searchtool/q-ranker.target.psms.txt "q-ranker q-value" > $searchtool/qvalues.qranker.txt
  echo replot \"$searchtool/qvalues.qranker.txt\" using 1:0 title \"$searchtool q-ranker\" with lines >> $gnuplot

  # Run Barista.
  # if [[ -e $searchtool/barista.target.psms.txt ]]; then
  #   echo Skipping barista.
  # else
  #   $CRUX barista \
  #     --output-dir $searchtool \
  #     --feature-file T \
  #     --overwrite T \
  #     --separate-searches $searchtool/$searchtool.decoy.txt \
  #     $fasta $ms2 $searchtool/$searchtool.target.txt
  # fi

  # $CRUX extract-columns $searchtool/barista.target.psms.txt "barista q-value" > $searchtool/qvalues.barista.txt
  # echo replot \"$searchtool/qvalues.barista.txt\" using 1:0 title \"$searchtool barista\" with lines >> $gnuplot
  
done

# Finalize the gnuplot script.
echo set output >> $gnuplot
echo replot >> $gnuplot

# Make the plot.
gnuplot $gnuplot > performance.png

# Extract a subset of columns for use in comparing XCorr scores.
for searchtool in tide-search comet; do
  searchFile=$searchtool/$searchtool.target.txt
  reducedFile=$searchtool/$searchtool.target.reduced.txt
  $CRUX extract-columns $searchFile \
     "scan,charge,sequence,xcorr score" \
     | awk 'NR > 1' \
     | awk '{print $1 "~" $2 "~" $3 "\t" $4}' \
     | sort -k 1b,1 \
     > $reducedFile
done

# Join the two sets of scores.
echo -e "scan\tcharge\tpeptide\ttide-search xcorr\tcomet xcorr" > xcorr.txt
join \
    tide-search/tide-search.target.reduced.txt \
    comet/comet.target.reduced.txt \
  | awk -F "~" '{print $1 " " $2 " " $3}' \
  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' \
  | sort -n \
  >> xcorr.txt

# Create a scatter plot of XCorr scores.
gnuplot=xcorr.gnuplot
echo set output \"/dev/null\" > $gnuplot
echo set terminal png >> $gnuplot
echo set xlabel \"Tide XCorr\" >> $gnuplot
echo set ylabel \"Comet XCorr\" >> $gnuplot
echo plot x notitle with lines >> $gnuplot
echo replot \"xcorr.txt\" using 4\:5 notitle >> $gnuplot
echo set output >> $gnuplot
echo replot >> $gnuplot
gnuplot $gnuplot > xcorr.png
