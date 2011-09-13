#!/bin/env python
# CREATE DATE: 13 Sep 2011
# AUTHOR: William Stafford Noble
import sys

PEPTIDE_COLUMN = "sequence"
PROTEIN_COLUMN = "protein id"

USAGE = """USAGE: ppnet.py [options] <file> <root>

  Takes as input a Crux text output file, and produces as output a
  series of plots showing connected components in the bipartite graph
  that connects peptides to proteins.

  The input file is tab-delimited, with peptides in a column named
  "%s" and comma-delimited lists of protein IDs in a column
  named "%s".

  The output is a series of pairs of files with names of the form
  <root>.<int>.gvz and <root>.<int>.png, where <root> is given on the
  command line and <int> is an ascending integer.  The gvz file
  contains a graphviz description of one component of the graph, and
  the png file contains a picture of the graph.

  Options:

    --min-component <int>   Skip components with fewer than <int> nodes.
""" % (PEPTIDE_COLUMN, PROTEIN_COLUMN)

###############################################################################
# MAIN
###############################################################################

# Parse the command line.
minComponent = 1
sys.argv = sys.argv[1:]
while (len(sys.argv) > 2):
  nextArg = sys.argv[1]
  sys.argv = sys.argv[1:]
  if (nextArg == "--min-component"):
    minComponent = int(sys.argv[0])
    sys.argv = sys.argv[1:]
  else:
    sys.stderr.write("Invalid option (%s).\n" % nextArg)
    sys.exit(1)
if (len(sys.argv) != 2):
  sys.stderr.write(USAGE)
  sys.exit(1)
inputFileName = sys.argv[0]
outputFileRoot = sys.argv[1]

# Read the header line and identify the target columns.
inputFile = open(inputFileName, "r")
headerLine = inputFile.readline().rstrip()
colIndex = 0
peptideColumn = -1
proteinColumn = -1
for column in headerLine.split("\t"):
  if (column == PEPTIDE_COLUMN):
    peptideColumn = colIndex
  elif (column == PROTEIN_COLUMN):
    proteinColumn = colIndex
  colIndex += 1
if (peptideColumn == -1):
  sys.stderr.write("Cannot find column with header %s.\n" % PEPTIDE_COLUMN)
  sys.exit(1)
sys.stderr.write("Reading peptides from column %d.\n" % peptideColumn)
if (proteinColumn == -1):
  sys.stderr.write("Cannot find column with header %s.\n" % PROTEIN_COLUMN)
  sys.exit(1)
sys.stderr.write("Reading protein IDs from column %d.\n" % proteinColumn)



# Print the graphviz files.

