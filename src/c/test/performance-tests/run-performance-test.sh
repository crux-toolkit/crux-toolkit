#!/bin/bash -efx
# AUTHOR: William Stafford Noble

# This script runs a performance test on the crux toolkit.  Its usage
# and output are described in ./performance.html.  The MS2 file comes
# from the PAnDA paper by Hoopman et al. (JPR 2009).  Extensive
# documentation of this particular data set is available from
# crux-projects/panda-data.

if [[ -e /etc/profile.d/modules.sh ]]; then
  source /etc/profile.d/modules.sh
  module load modules modules-init modules-gs modules-noble
  module load mpc/0.8.2 mpfr/3.0.0 gmp/5.0.2 gcc/4.8.1
  module load protobuf/2.5.0
fi

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

# Create the parameter file.
parameters=crux.param

# Enzymatic digestion rules.
echo enzyme=trypsin >> $parameters
echo search_enzyme_number=1 >> $parameters 
echo digestion=full-digest > $parameters
echo num_enzyme_termini=2 >> $parameters
echo missed-cleavages=0 >> $parameters
echo allowed_missed_cleavage=0 >> $parameters

# Minimums
echo minimum_peaks=10 >> $parameters
echo min-peaks=10 >> $parameters
echo min_length=2 >> $parameters

# Precursor selection rules.
echo precursor-window=3 >> $parameters
echo precursor-window-type=mass >> $parameters
echo peptide_mass_tolerance=3 >> $parameters
echo peptide_mass_units=0 >> $parameters # 0=amu, 1=mmu, 2=ppm
echo precursor_tolerance_type=0 >> $parameters # 0=MH+ (default), 1=precursor m/z

# Precursor mass type.
echo isotopic-mass=mono >> $parameters
echo monoisotopic-precursor=T >> $parameters
echo mass_type_parent=1 >> $parameters # 1=monoisotopic

# Fragment mass type.  Tides uses only monoisotopic.
echo fragment-mass=mono >> $parameters
echo mass_type_fragment=1 >> $parameters # 1=monoisotopic

# Decoys.
echo decoy-format=peptide-reverse >> $parameters
echo num-decoys-per-target=1 >> $parameters
echo keep-terminal-aminos=C >> $parameters
echo decoy_search=2 >> $parameters  # 2 = separate decoy search

# Report the top 5 matches.
echo num_results=6 >> $parameters
echo num_output_lines=5 >> $parameters
echo top-match=5 >> $parameters

# Precursor removal.
echo remove_precursor_peak=1 >> $parameters
echo remove_precursor_tolerance=15 >> $parameters
echo remove-precursor-peak=T >> $parameters
echo remove-precursor-tolerance=15 >> $parameters

# Flanking peaks.
echo use-flanking-peaks=F >> $parameters
echo theoretical_fragment_ions=1 >> $parameters # 0 = flanks; 1 = no flanks
echo use-neutral-loss-peaks=F >> $parameters 
# Fragment m/z discretization.  This is fixed in Tide.
echo fragment_bin_offset=0.4 >> $parameters
echo fragment_bin_tol=1.0005079 >> $parameters
echo mz-bin-offset=0.4 >>$parameters
echo mz-bin-width=1.0005079 >>$parameters

# Other Crux parameters.
echo compute-sp=F >> $parameters
echo verbosity=40 >> $parameters
echo overwrite=T >> $parameters
echo peptide-list=T >> $parameters

# Comet parameters
echo add_C_cysteine=57.021464 >> $parameters
echo num_threads=1 >> $parameters # Multithreaded sometimes dumps core.
echo digest_mass_range=200 7200 >> $parameters
echo max_fragment_charge=2 >> $parameters
echo isotope_error=0 >> $parameters
echo use_A_ions=0 >> $parameters
echo use_B_ions=1 >> $parameters
echo use_C_ions=0 >> $parameters
echo use_X_ions=0 >> $parameters
echo use_Y_ions=1 >> $parameters
echo use_Z_ions=0 >> $parameters
echo use_NL_ions=0 >> $parameters
echo variable_mod1=0.0 X 0 3 >> $parameters
echo variable_mod2=0.0 X 0 3 >> $parameters
echo "[COMET_ENZYME_INFO]" >> $parameters
echo "0.  No_enzyme              0      -           -" >> $parameters
echo "1.  Trypsin                1      KR          P" >> $parameters

# Create the index.
db=worm+contaminants
fasta=$db.fa
if [[ -e $db ]]; then
  echo Skipping create-index.
else
  $CRUX tide-index --output-dir tide-index --parameter-file $parameters \
     $fasta $db
fi

ms2=051708-worm-ASMS-10.ms2

# Do the whole test twice, once for each search tool.
for searchtool in comet tide-search; do

  # Do we use an index or the fasta?
  if [[ $searchtool == "comet" ]]; then
    proteins=$fasta
  else
    proteins=$db
  fi

  # Run the search.
  if [[ -e $searchtool/$searchtool.target.txt ]]; then
    echo Skipping search-for-matches.
  else  
    $CRUX $searchtool \
      --parameter-file crux.param --output-dir $searchtool \
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

  $CRUX sort-by-column --column-type real --ascending true $searchtool/qvalues.target.txt "decoy q-value (xcorr)" | $CRUX extract-columns - "decoy q-value (xcorr)" > $searchtool/qvalues.xcorr.txt
  echo replot \"$searchtool/qvalues.xcorr.txt\" using 1:0 title \"$searchtool xcorr\" with lines >> $gnuplot

  if [[ $searchtool == "comet" ]]; then
    $CRUX extract-columns $searchtool/qvalues.target.txt "decoy q-value (e-value)" > $searchtool/qvalues.evalue.txt 
    echo replot \"$searchtool/qvalues.evalue.txt\" using 1:0 title \"$searchtool e-value\" with lines >> $gnuplot
  fi

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
