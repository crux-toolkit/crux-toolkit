#!/bin/env python
# CREATE DATE: 11 March 2015
# AUTHOR: William Stafford Noble
# Derived from run-performance-test.sh.

# This script runs a performance test on the crux toolkit.  Its usage
# and output are described in ./performance.html.  The MS2 file comes
# from the PAnDA paper by Hoopman et al. (JPR 2009).  Extensive
# documentation of this particular data set is available from
# crux-projects/panda-data.
import sys
import subprocess
import os

# The location of the crux binary.
CRUX = "../../src/crux"

# Input files.
database = "worm+contaminants"
ms2 = "051708-worm-ASMS-10.ms2"

#############################################################################
# Run a command with error checking.
def runCommand(command, outputFileName):

  # Skip the command if the output file already exists.
  if (outputFileName != "") and (os.path.exists(outputFileName)):
    sys.stderr.write("%s exists.\n" % outputFileName)
    return
  
  sys.stderr.write("RUN: %s\n" % command)
  try:
    returnCode = subprocess.call(command, shell=True)
    if (returnCode != 0):
      sys.stderr.write("Child was terminated by signal %d\n" % -returnCode)
      sys.exit(1)
  except OSError, e:
    sys.stderr.write("Execution failed: %s\n" % e)
    sys.exit(1)

#############################################################################
def createParameterFile(parameterFileName):
  parameterFile = open(parameterFileName, "w")

  # Enzymatic digestion rules.
  parameterFile.write("enzyme=trypsin\n")
  parameterFile.write("search_enzyme_number=1\n") 
  parameterFile.write("digestion=full-digest\n")
  parameterFile.write("num_enzyme_termini=2\n")
  parameterFile.write("missed-cleavages=0\n")
  parameterFile.write("allowed_missed_cleavage=0\n")
  
  # Minimums
  parameterFile.write("minimum_peaks=10\n")
  parameterFile.write("min-peaks=10\n")
  
  # Precursor selection rules.
  parameterFile.write("precursor-window=3\n")
  parameterFile.write("precursor-window-type=mass\n")
  parameterFile.write("peptide_mass_tolerance=3\n")
  parameterFile.write("peptide_mass_units=0\n") # 0=amu, 1=mmu, 2=ppm
  
  # Precursor mass type.
  parameterFile.write("isotopic-mass=mono\n")
  parameterFile.write("monoisotopic-precursor=T\n")
  parameterFile.write("mass_type_parent=1\n") # 1=monoisotopic
  
  # Fragment mass type.  Tides uses only monoisotopic.
  parameterFile.write("fragment-mass=mono\n")
  parameterFile.write("mass_type_fragment=1\n") # 1=monoisotopic
  
  # Decoys.
  parameterFile.write("decoy-format=peptide-reverse\n")
  parameterFile.write("num-decoys-per-target=1\n")
  parameterFile.write("keep-terminal-aminos=C\n")
  parameterFile.write("decoy_search=1\n")  # 1 = concatenated decoy search
  
  # Report the top 5 matches.
  parameterFile.write("num_results=6\n")
  parameterFile.write("num_output_lines=5\n")
  parameterFile.write("top-match=5\n")
  
  # Precursor removal.
  parameterFile.write("remove_precursor_peak=1\n")
  parameterFile.write("remove_precursor_tolerance=15\n")
  parameterFile.write("remove-precursor-peak=T\n")
  parameterFile.write("remove-precursor-tolerance=15\n")
  
  # Flanking peaks.
  parameterFile.write("use-flanking-peaks=F\n")
  parameterFile.write("theoretical_fragment_ions=1\n") # 0 = flanks; 1 = no flanks
  parameterFile.write("use-neutral-loss-peaks=F\n") 

  # Fragment m/z discretization.
  parameterFile.write("fragment_bin_offset=0.68\n")
  parameterFile.write("fragment_bin_tol=1.0005079\n")
  parameterFile.write("mz-bin-offset=0.68\n")
  parameterFile.write("mz-bin-width=1.0005079\n")
  
  # Other Crux parameters.
  parameterFile.write("concat=T\n")
  parameterFile.write("compute-sp=T\n")
  parameterFile.write("verbosity=40\n")
  parameterFile.write("overwrite=T\n")
  parameterFile.write("peptide-list=T\n")
  
  # Comet parameters
  parameterFile.write("output_pepxml=0\n")
  parameterFile.write("add_C_cysteine=57.021464\n")
  parameterFile.write("num_threads=1\n") # Multithreaded sometimes dumps core.
  parameterFile.write("digest_mass_range=200 7200\n")
  parameterFile.write("max_fragment_charge=2\n")
  parameterFile.write("isotope_error=0\n")
  parameterFile.write("use_A_ions=0\n")
  parameterFile.write("use_B_ions=1\n")
  parameterFile.write("use_C_ions=0\n")
  parameterFile.write("use_X_ions=0\n")
  parameterFile.write("use_Y_ions=1\n")
  parameterFile.write("use_Z_ions=0\n")
  parameterFile.write("use_NL_ions=0\n")
  parameterFile.write("variable_mod01=0.0 X 0 3\n")
  parameterFile.write("variable_mod02=0.0 X 0 3\n")
  parameterFile.write("[COMET_ENZYME_INFO]\n")
  parameterFile.write("0.  No_enzyme              0      -           -\n")
  parameterFile.write("1.  Trypsin                1      KR          P\n")
  parameterFile.close()

#############################################################################
def extractData(inputFileName, columnName, outputFileName):
  runCommand("%s extract-columns %s \"%s\" > %s" % 
             (CRUX, inputFileName, columnName, outputFileName), "")

#############################################################################
def runSearch(outputDirectory, searchName, searchParam, database, 
              psmFile, scoreColumn, confidenceParam):

  runCommand("%s %s --output-dir %s --parameter-file %s %s %s %s"
             % (CRUX, searchName, outputDirectory, parameterFileName, 
                searchParam, ms2, database),
             psmFile)

  confidenceFile = "%s/assign-confidence.target.txt" % outputDirectory
  runCommand("%s assign-confidence --output-dir %s %s %s" % 
             (CRUX, outputDirectory, confidenceParam, psmFile), confidenceFile)

  qFile = "%s/%s.q.txt" % (outputDirectory, searchName)
  extractData(confidenceFile, "tdc q-value", qFile)

  percolatorFile = "%s/percolator.target.psms.txt" % outputDirectory
  runCommand("%s percolator --output-dir %s %s"
             % (CRUX, outputDirectory, psmFile), percolatorFile)

  qFile = "%s/%s.percolator.q.txt" % (outputDirectory, searchName)
  extractData(percolatorFile, "percolator q-value", qFile)

  qrankerFile = "%s/q-ranker.target.psms.txt" % outputDirectory
  runCommand("%s q-ranker --decoy-prefix decoy_ --output-dir %s %s %s"
             % (CRUX, outputDirectory, ms2, psmFile), qrankerFile)

  qFile = "%s/%s.q-ranker.q.txt" % (outputDirectory, searchName)
  extractData(qrankerFile, "q-ranker q-value", qFile)

  reducedFile = "%s/%s.target.reduced.txt" % (outputDirectory, searchName)
  runCommand("%s extract-columns %s \"scan,charge,sequence,%s\" | awk 'NR > 1' | awk '{print $1 \"~\" $2 \"~\" $3 \"\t\" $4}' | sort -k 1b,1 > %s"
             % (CRUX, psmFile, scoreColumn, reducedFile), "")


# Create a scatter plot of XCorr scores.
def makeScatterPlot(xData, xLabel, yData, yLabel, outputRoot):

  runCommand("join %s %s | awk -F \"~\" '{print $1 \" \" $2 \" \" $3}' | awk '{print $1 \"\t\" $2 \"\t\" $3 \"\t\" $4 \"\t\" $5}' | sort -n > %s.txt"
             % (xData, yData, outputRoot), "")

  gnuplotFileName = "%s.gnuplot" % outputRoot
  gnuplotFile = open(gnuplotFileName, "w")
  gnuplotFile.write("set output \"/dev/null\"\n")
  gnuplotFile.write("set terminal png\n")
  gnuplotFile.write("set xlabel \"%s\"\n" % xLabel)
  gnuplotFile.write("set ylabel \"%s\"\n" % yLabel)
  gnuplotFile.write("plot x notitle with lines\n")
  gnuplotFile.write("replot \"%s.txt\" using 4:5 notitle\n" % outputRoot)
  gnuplotFile.write("set output\n")
  gnuplotFile.write("replot\n")
  gnuplotFile.close()

  runCommand("gnuplot %s > %s.png" % (gnuplotFileName, outputRoot), "")


#############################################################################
# MAIN
#############################################################################

# Create the parameter file.
parameterFileName = "crux.param"
createParameterFile(parameterFileName)

# Create the index.
runCommand("%s tide-index --output-dir %s --parameter-file %s %s.fa %s"
           % (CRUX, database, parameterFileName, database, database), 
           "%s/tide-index.peptides.target.txt" % database)

# Run three searches (Comet, Tide XCorr, and Tide p-value).
runSearch("tide-xcorr", "tide-search", "--exact-p-value T", database, 
          "tide-xcorr/tide-search.txt", "xcorr score", "")
runSearch("tide-p-value", "tide-search", "", database,
          "tide-p-value/tide-search.txt", "refactored xcorr",
          "--smaller-is-better T --score \"exact p-value\"")
runSearch("comet", "comet", "", "%s.fa" % database, 
          "comet/comet.target.txt", "xcorr score", "")
# FIXME: The last option should be "--score \"e-value\""

# Make the performance plot.
gnuplotFileName = "performance.gnuplot"
gnuplotFile = open(gnuplotFileName, "w")
gnuplotFile.write("set output \"/dev/null\"\n")
gnuplotFile.write("set terminal png\n")
gnuplotFile.write("set xlabel \"q-value threshold\"\n")
gnuplotFile.write("set ylabel \"Number of accepted PSMs\"\n")
gnuplotFile.write("set xrange [0:0.1]\n")
gnuplotFile.write("set key center right\n")
gnuplotFile.write("plot \"comet/comet.q.txt\" using 1:0 title \"Comet E-value\" with lines\n")
gnuplotFile.write("replot \"tide-p-value/tide-search.q.txt\" using 1:0 title \"Tide p-value\" with lines\n")
gnuplotFile.write("replot \"tide-xcorr/tide-search.q.txt\" using 1:0 title \"Tide XCorr\" with lines\n")
gnuplotFile.write("replot \"comet/comet.percolator.q.txt\" using 1:0 title \"Comet Percolator\" with lines\n")
gnuplotFile.write("replot \"tide-p-value/tide-search.percolator.q.txt\" using 1:0 title \"Tide p-value Percolator\" with lines\n")
gnuplotFile.write("replot \"tide-xcorr/tide-search.percolator.q.txt\" using 1:0 title \"Tide XCorr Percolator\" with lines\n")
gnuplotFile.write("replot \"comet/comet.q-ranker.q.txt\" using 1:0 title \"Comet q-ranker\" with lines\n")
gnuplotFile.write("replot \"tide-p-value/tide-search.q-ranker.q.txt\" using 1:0 title \"Tide p-value q-ranker\" with lines\n")
gnuplotFile.write("replot \"tide-xcorr/tide-search.q-ranker.q.txt\" using 1:0 title \"Tide XCorr q-ranker\" with lines\n")
gnuplotFile.write("set output\n")
gnuplotFile.write("replot\n")
gnuplotFile.close() 
runCommand("gnuplot %s > performance.png" % gnuplotFileName, "")


# Make the XCorr scatter plots.
makeScatterPlot("tide-xcorr/tide-search.target.reduced.txt", 
                "Tide XCorr",
                "tide-p-value/tide-search.target.reduced.txt", 
                "Refactored XCorr",
                "xcorr.refactored")
makeScatterPlot("tide-xcorr/tide-search.target.reduced.txt", 
                "Tide XCorr",
                "comet/comet.target.reduced.txt", 
                "Comet XCorr",
                "xcorr.comet")


