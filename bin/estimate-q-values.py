#!/usr/bin/env python
# AUTHOR: William Stafford Noble
# CREATE DATE: 21 April 2014
import sys
import math

USAGE = """USAGE: estimate-q-values.py [options] <file>

  Read a Crux tab-delimited file of PSMs and estimate a q-value for
  each PSM.  This can be done using Benjamini-Hochberg if p-values are
  provided, or using target-decoy competition if a set of decoys are
  provided.

  Options:

    --method tdc|b-h       Specify whether to use Benjamini-Hochberg or
                           target-decoy competition to estimate q-values.
                           Default = b-h.

    --decoys <file>        Provide a parallel file containing the decoy
                           scores.  Default = no file.

    --score-column <name>  Specify the name of the column containing
                           the scores.  Default = "XCorr p-value"

    --best-score high|low  Specify whether high-scoring or low-scoring
                           PSMs are considered good matches.
                           Default = low.

    --rank-column <name>   Specify the name of the column containing
                           the ranks.  The program will ignore all PSMs
                           for which the value in this column is not 1.
                           Default = no rank column.

"""

#############################################################################
def estimateQvalues(estimationMethod, sortedPSMs):
  returnValue = [] # List of tuples.

  # Estimate FDRs using target-decoy competition.
  if (estimationMethod == "tdc"):
    numDecoys = 1
    numTargets = 0
    for index in range(0, len(sortedPSMs)):
      (score, isTarget, line) = sortedPSMs[index]
  
      if isTarget:
        numTargets += 1
        if (numTargets == 0):
          fdr = 1.0
        else:
          fdr = float(numDecoys) / float(numTargets)
#          sys.stderr.write("%d/%d=%g (%g)\n" % (numDecoys, numTargets, fdr, score));
        if (fdr > 1.0):
          fdr = 1.0
        returnValue.append((score, line, fdr))
      else:
#        sys.stderr.write("decoy (%g)\n" % score);
        numDecoys += 1

  # Estimate FDRs using target-decoy competition.
  else:
    sys.stderr.write("Sorry! Benjamini-Hochberg is not yet implemented.\n")
    sys.exit(1)

  # Compute q-values.
  minFDR = 1.0
  for index in range(len(returnValue)-1, -1, -1):
    (score, line, fdr) = returnValue[index]
    if (fdr < minFDR):
      minFDR = fdr
    returnValue[index] = (score, line, minFDR)

  return(returnValue)

#############################################################################
# MAIN
#############################################################################

# Set default parameters.
estimationMethod = "b-h"
decoyFileName = ""
scoreColumn = "XCorr p-value"
highIsGood = False

# Parse the optional arguments.
sys.argv = sys.argv[1:]
while (len(sys.argv) > 1):
  nextArg = sys.argv[0]
  sys.argv = sys.argv[1:]
  if (nextArg == "--method"):
    if (sys.argv[0] in ["tdc", "b-h"]):
      estimationMethod = sys.argv[0]
    else:
      sys.stderr.write("Unrecognized estimation method (%s).\n" %
                       sys.argv[0])
      sys.exit(1)
    sys.argv = sys.argv[1:]
  elif (nextArg == "--decoys"):
    decoyFileName = sys.argv[0]
    sys.argv = sys.argv[1:]
  elif (nextArg == "--score-column"):
    scoreColumn = sys.argv[0]
    sys.argv = sys.argv[1:]
  elif (nextArg == "--best-score"):
    if (sys.argv[0] == "high"):
      highIsGood = True
    elif (sys.argv[0] == "low"):
      highIsGood = False
    else:
      sys.stderr.write("Unrecognized best-score option (%s).\n" %
                       sys.argv[0])
      sys.exit(1)
    sys.argv = sys.argv[1:]
  elif (nextArg == "--rank-column"):
    sys.stderr.write("--rank-column is not yet implemented.\n")
    sys.exit(1)
  else:
    sys.stderr.write("Unrecognized option (%s).\n" % nextArg)
    sys.exit(1)

# Parse the required argument.
if (len(sys.argv) != 1):
  sys.stderr.write(USAGE)
  sys.exit(1)
psmFileName = sys.argv[0]

# Complain about parameter mismatches.
if (decoyFileName != "") and (estimationMethod == "b-h"):
  sys.stderr.write("Error: You asked for Benjamini-Hochberg estimates but provided decoys.\n")
  sys.stderr.exit(1)
if (decoyFileName == "") and (estimationMethod == "tdc"):
  sys.stderr.write("Error: You asked for target-decoy estimates but did not provide decoys.\n")
  sys.stderr.exit(1)
if (highIsGood and (estimationMethod == "b-h")):
  sys.stderr.write("Benjamini-Hochberg is incompatible with --best-score high.\n")
  sys.exit(1)

# Parse the header line and identify the target column.
psmFile = open(psmFileName, "r")
headerLine = psmFile.readline().rstrip()

# Locate the desired columns.
try:
  scanColumnIndex = headerLine.split("\t").index("scan")
except:
  sys.stderr.write("Can't find column \"scan.\"\n")
  sys.exit(1)
sys.stderr.write("Reading scan numbers from column %d.\n" % scanColumnIndex)
try:
  scoreColumnIndex = headerLine.split("\t").index(scoreColumn)
except:
  sys.stderr.write("Can't find column \"%s.\"\n" % scoreColumn)
  sys.exit(1)
sys.stderr.write("Reading scores from column %d (\"%s\").\n" % 
                 (scoreColumnIndex, scoreColumn))

# If a decoy file was provided, open it.
if (decoyFileName != ""):
  decoyFile = open(decoyFileName, "r")
  decoyFile.readline()

# Read line by line.
psms = [] # (score, isTargetScore, line)
lineNumber = 0
numDecoys = 0
for line in psmFile:
  lineNumber += 1
  targetWords = line.rstrip().split("\t")

  score = float(targetWords[scoreColumnIndex])
  isTargetScore = True

  if (decoyFileName != ""):
    decoyWords = decoyFile.readline().rstrip().split("\t")

    # Check that the target and decoy scan numbers match.
    if (targetWords[scanColumnIndex] != decoyWords[scanColumnIndex]):
      sys.stderr.write("Target-decoy scan number mismatch (%s != %s) at line %d.\n" %
                       (targetWords[scanColumnIndex],
                        decoyWords[scanColumnIndex],
                        lineNumber))
      continue

    # Do the target-decoy competition.
    decoyScore = float(decoyWords[scoreColumnIndex])
#    sys.stdout.write("%g versus %g\n" % (score, decoyScore))
    if ((not highIsGood and (decoyScore < score)) or
        (highIsGood and (decoyScore > score))):
      score = decoyScore
      isTargetScore = False
      numDecoys += 1

  psms.append((score, isTargetScore, line.strip()))
psmFile.close()
sys.stderr.write("Read %d PSMs from %s.\n" % (len(psms), psmFileName))
if (decoyFileName != ""):
  sys.stderr.write("%d decoys won target-decoy competition.\n" % numDecoys)

# Sort by score.
sortedPSMs = sorted(psms, key=lambda tup: tup[0], reverse=highIsGood)
sys.stderr.write("Scores range from %g to %g.\n" 
                 % (sortedPSMs[0][0], sortedPSMs[-1][0]))

# Output is list of tuples: (score, line, q-value)
reportedPSMs = estimateQvalues(estimationMethod, sortedPSMs)

# Print!
sys.stdout.write("%s\tq-value\n" % (headerLine))
for psm in reportedPSMs:
  sys.stdout.write("%s\t%g\n" % (psm[1], psm[2]))
