#!/usr/bin/env python
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
import pandas as pd

# The location of the crux binary.
CRUX = "../../src/crux"

# Input files.
database = "worm+contaminants"
ms2 = "051708-worm-ASMS-10.ms2"
dia_index = "cerevisiae_orf_trans_all"
dia_ms = "e01306.mzXML"
#############################################################################
# Run a command with error checking.
def runCommand(command, outputFileName, overwrite=False):

  # Skip the command if the output file already exists.
  if overwrite == False: 
    if (outputFileName != "") and (os.path.exists(outputFileName)):
      sys.stderr.write("%s exists.\n" % outputFileName)
      return
  
  sys.stderr.write("RUN: %s\n" % command)
  try:
    returnCode = subprocess.call(command, shell=True)
    if (returnCode != 0):
      sys.stderr.write("Child was terminated by signal %d\n" % -returnCode)
      sys.exit(1)
  except OSError as e:
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
  parameterFile.write("spectrum-parser=mstoolkit\n")

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
  parameterFile.write("mass_type_parent=1\n") # 1=monoisotopic
  
  # Fragment mass type.  Tides uses only monoisotopic.
  parameterFile.write("fragment-mass=mono\n")
  parameterFile.write("mass_type_fragment=1\n") # 1=monoisotopic
  
  # Decoys.
  parameterFile.write("decoy-format=peptide-reverse\n")
  parameterFile.write("decoy_search=1\n")  # 1 = concatenated decoy search
  parameterFile.write("concat=T\n")
  parameterFile.write("keep-terminal-aminos=C\n") # No corresponding Comet param
  
  # Report the top 5 matches.
  parameterFile.write("top-match=5\n")
  parameterFile.write("num_results=6\n")
  parameterFile.write("num_output_lines=5\n")
  
  # Precursor removal.
  parameterFile.write("remove-precursor-peak=T\n")
  parameterFile.write("remove-precursor-tolerance=15\n")
  parameterFile.write("remove_precursor_peak=1\n")
  parameterFile.write("remove_precursor_tolerance=15\n")
  
  # Flanking peaks.
  parameterFile.write("use-flanking-peaks=F\n")
  parameterFile.write("theoretical_fragment_ions=1\n") # 0 = flanks; 1 = no flanks
  parameterFile.write("use-neutral-loss-peaks=F\n") 

  # Fragment m/z discretization.
  parameterFile.write("mz-bin-offset=0.68\n")
  parameterFile.write("mz-bin-width=1.0005079\n")
  parameterFile.write("fragment_bin_offset=0.68\n")
  parameterFile.write("fragment_bin_tol=1.0005079\n")
  
  # Peptide mass range.
  parameterFile.write("min-mass=200\n")
  parameterFile.write("max-mass=7200\n")
  parameterFile.write("digest_mass_range=200 7200\n")
  
  # Other Crux parameters.
  parameterFile.write("num-threads=12\n")
  parameterFile.write("compute-sp=T\n")
  parameterFile.write("overwrite=T\n")
  parameterFile.write("peptide-list=T\n")
  
  # Comet parameters
  parameterFile.write("output_pepxmlfile=0\n")
  parameterFile.write("add_C_cysteine=57.021464\n")
  parameterFile.write("num_threads=12\n") 
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
#  dataset = pd.read_csv(inputFileName, sep='\t')
#  data_frame = pd.DataFrame(dataset)
  data_frame = pd.read_csv(inputFileName, sep='\t')
  data_frame = data_frame[[columnName]]
  data_frame.to_csv(outputFileName, sep='\t', index=False)

#############################################################################
def runSearch(outputDirectory, searchName, searchParam, database, 
              psmFile, scoreColumn, confidenceParam):

  runCommand("%s %s --output-dir %s --parameter-file %s %s %s %s"
             % (CRUX, searchName, outputDirectory, parameterFileName, 
                searchParam, ms2, database),
             psmFile)

  percolatorFile = "%s/percolator.target.psms.txt" % outputDirectory
  runCommand("%s percolator --output-dir %s %s"
             % (CRUX, outputDirectory, psmFile), percolatorFile)

def calculateQvalues(outputDirectory, searchName, searchParam, database, 
              psmFile, scoreColumn, confidenceParam):    
               
  confidenceFile = "%s/assign-confidence.target.txt" % outputDirectory
  runCommand("%s assign-confidence --overwrite T --output-dir %s %s %s" % 
             (CRUX, outputDirectory, confidenceParam, psmFile), confidenceFile, overwrite=True)

  qFile = "%s/%s.q.txt" % (outputDirectory, searchName)
  extractData(confidenceFile, "tdc q-value", qFile)

  percolatorFile = "%s/percolator.target.psms.txt" % outputDirectory
  qFile = "%s/%s.percolator.q.txt" % (outputDirectory, searchName)
  extractData(percolatorFile, "q-value", qFile)

  fourColFile = "%s/%s.target.four-col.txt" % (outputDirectory, searchName)
  dataset = pd.read_csv(psmFile, sep='\t')
  data_frame = pd.DataFrame(dataset)
  data_frame = data_frame[["scan", "charge", "sequence", scoreColumn]]
  data_frame.to_csv(fourColFile, sep='\t', index=False)

  reducedFile = "%s/%s.target.reduced.txt" % (outputDirectory, searchName)
  runCommand("awk 'NR > 1 {print $1 \"~\" $2 \"~\" $3 \"\t\" $4}' %s | LANG=en_EN sort -k 1b,1 > %s" %
             (fourColFile, reducedFile), "")

def getPercolatorQvalues(outputDirectory, searchName):    
               
  percolatorFile = "%s/percolator.target.psms.txt" % outputDirectory
  qFile = "%s/%s.percolator.q.txt" % (outputDirectory, searchName)
  extractData(percolatorFile, "q-value", qFile)

# Create a scatter plot of XCorr scores.
def makeScatterPlot(xData, xLabel, yData, yLabel, outputRoot):

  runCommand("LANG=en_EN join %s %s | awk -F \"~\" '{print $1 \" \" $2 \" \" $3}' | awk '{print $1 \"\t\" $2 \"\t\" $3 \"\t\" $4 \"\t\" $5}' | sort -n > %s.txt"
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

# Create a gnuplot of a list of methods.
def makePerformancePlot(title, listOfMethods):
  if not os.path.isdir("plots"):
    os.mkdir("plots")
  gnuplotFileName = "plots/%s.gnuplot" % title
  gnuplotFile = open(gnuplotFileName, "w")
  gnuplotFile.write("set output \"/dev/null\"\n")
  gnuplotFile.write("set terminal png\n")
  gnuplotFile.write("set title \"%s\"\n" % title)
  gnuplotFile.write("set xlabel \"q-value threshold\"\n")
  gnuplotFile.write("set ylabel \"Number of accepted PSMs\"\n")
  gnuplotFile.write("set xrange [0:0.1]\n")
  gnuplotFile.write("set key bottom right\n")
  firstOne = True
  for myTuple in listOfMethods:
    (dataFileName, seriesTitle) = myTuple

    if firstOne:
      firstOne = False
    else:
      gnuplotFile.write("re")

    gnuplotFile.write("plot \"%s\" using 1:0 title \"%s\" with lines lw 1\n" 
                      % (dataFileName, seriesTitle))

  gnuplotFile.write("set output\n")
  gnuplotFile.write("replot\n")
  gnuplotFile.close() 
  runCommand("gnuplot %s > plots/%s.png" % (gnuplotFileName, title), "")
  
  
#############################################################################
# MAIN
#############################################################################

# Create the parameter file.
parameterFileName = "crux.param"
# createParameterFile(parameterFileName)

# # Create the index for diameter.
# runCommand("%s tide-index --output-dir %s --parameter-file %s %s.fasta %s"
#            % (CRUX, dia_index, parameterFileName, dia_index, dia_index),
#            "%s/tide-index.peptides.txt" % dia_index)

# # Run diameter.
# runSearch("diameter", "diameter", "", dia_index,
#           "diameter/percolator-output/percolator.target.psms.txt",
#           "tailor score", "--score \"tailor score\"")
# getPercolatorQvalues("diameter", "diameter")

# Add performance tests for cascade search too.


# Create the index.
runCommand("%s tide-index --output-dir %s --parameter-file %s %s.fa %s"
           % (CRUX, database, parameterFileName, database, database),
           "%s/tide-index.peptides.txt" % database)

# Run searches (Comet, and Tide XCorr (includes tailor), XCorr combined (includes exact p-value and res-ev).

runSearch("tide-xcorr", "tide-search", "", database,
          "tide-xcorr/tide-search.txt",
          "xcorr score", "--score \"xcorr score\"")

calculateQvalues("tide-xcorr", "tide-xcorr", "", database,
          "tide-xcorr/tide-search.txt",
          "xcorr score", "--score \"xcorr score\"")

calculateQvalues("tide-xcorr", "tide-tailor", "", database,
          "tide-xcorr/tide-search.txt",
          "tailor score", "--score \"tailor score\"")

runSearch("tide-combined", "tide-search",
          "--score-function combined-p-values --fragment-tolerance 1.0005079",
          database, "tide-combined/tide-search.txt",
          "combined p-value", "--score \"combined p-value\"")

calculateQvalues("tide-combined", "tide-combined",
          "--score-function combined-p-values --fragment-tolerance 1.0005079",
          database, "tide-combined/tide-search.txt",
          "combined p-value", "--score \"combined p-value\"")

calculateQvalues("tide-combined", "tide-xpv",
          "--score-function combined-p-values --fragment-tolerance 1.0005079",
          database, "tide-combined/tide-search.txt",
          "exact p-value", "--score \"exact p-value\"")

calculateQvalues("tide-combined", "tide-refact",
          "--score-function combined-p-values --fragment-tolerance 1.0005079",
          database, "tide-combined/tide-search.txt",
          "refactored xcorr", "--score \"exact p-value\"")

calculateQvalues("tide-combined", "tide-res-ev",
          "--score-function combined-p-values --fragment-tolerance 1.0005079",
          database, "tide-combined/tide-search.txt",
          "res-ev score", "--score \"res-ev p-value\"")

runSearch("comet", "comet", "", "%s.fa" % database,
          "comet/comet.txt",
          "xcorr score", "--score e-value")

calculateQvalues("comet", "comet", "", "%s.fa" % database,
          "comet/comet.txt",
          "xcorr score", "--score e-value")


# Make the performance plots, segregated by search method..
makePerformancePlot("comet", [("comet/comet.q.txt", "Comet E-value"),
                     ("comet/comet.percolator.q.txt", "Comet Percolator")])
makePerformancePlot("tide.xcorr",
                    [("tide-xcorr/tide-xcorr.q.txt", "Tide XCorr"),
                     ("tide-xcorr/tide-xcorr.percolator.q.txt", "Tide XCorr Percolator")])
makePerformancePlot("tide.tailor",
                    [("tide-xcorr/tide-tailor.q.txt", "Tide Tailor "),
                     ("tide-xcorr/tide-tailor.percolator.q.txt", "Tide Tailor Percolator")])

makePerformancePlot("tide.p-value",
                    [("tide-combined/tide-xpv.q.txt", "Tide p-value"),
                     ("tide-combined/tide-xpv.percolator.q.txt", "Tide p-value Percolator")])
makePerformancePlot("tide.res-ev",
                    [("tide-combined/tide-res-ev.q.txt", "Tide res-ev"),
                     ("tide-combined/tide-res-ev.percolator.q.txt", "Tide res-ev Percolator")])
makePerformancePlot("tide.combined",
                    [("tide-combined/tide-combined.q.txt", "Tide combined p-value"),
                     ("tide-combined/tide-combined.percolator.q.txt", "Tide combined p-value Percolator")])

# Make the performance plots, segregated by post-processor.
makePerformancePlot("assign-confidence",
                    [("comet/comet.q.txt", "Comet E-value"),
                     ("tide-combined/tide-xpv.q.txt", "Tide p-value"),
                     ("tide-xcorr/tide-xcorr.q.txt", "Tide XCorr"),
                     ("tide-combined/tide-res-ev.q.txt", "Tide res-ev"),
                     ("tide-combined/tide-combined.q.txt", "Tide combined p-value"),
                     ("tide-xcorr/tide-tailor.q.txt", "Tide Tailor")])
makePerformancePlot("percolator",
                    [("comet/comet.percolator.q.txt", "Comet Percolator"),
                     ("tide-combined/tide-xpv.percolator.q.txt", "Tide p-value Percolator"),
                     ("tide-xcorr/tide-xcorr.percolator.q.txt", "Tide XCorr Percolator"),
                     ("tide-combined/tide-res-ev.percolator.q.txt", "Tide res-ev Percolator"),
                     ("tide-combined/tide-combined.percolator.q.txt", "Tide combined p-value Percolator"),
                     ("tide-xcorr/tide-tailor.percolator.q.txt", "Tide Tailor Percolator")])


# Make the XCorr scatter plots.
makeScatterPlot("tide-xcorr/tide-xcorr.target.reduced.txt", 
                "Tide XCorr",
                "comet/comet.target.reduced.txt", 
                "Comet XCorr",
                "plots/xcorr.comet")
makeScatterPlot("tide-xcorr/tide-xcorr.target.reduced.txt", 
                "Tide XCorr",
                "tide-combined/tide-refact.target.reduced.txt", 
                "Refactored XCorr",
                "plots/xcorr.refactored")
makeScatterPlot("tide-combined/tide-res-ev.target.reduced.txt", 
                "Current Tide Res-Ev score",
                "./tide-res-ev.target.reduced.txt", 
                "Stored Tide Res-Ev score",
                "plots/resev")


