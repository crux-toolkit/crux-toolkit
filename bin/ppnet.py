#!/usr/bin/env python
# CREATE DATE: 13 Sep 2011
# AUTHOR: William Stafford Noble
import sys
import os
import subprocess

PEPTIDE_COLUMN = "sequence"
PROTEIN_COLUMN = "protein id"

USAGE = """USAGE: ppnet.py [options] <file> <root>

  Takes as input a Crux text output file, and produces as output a
  series of plots showing connected components in the bipartite graph
  that connects peptides to proteins.

  The input file is tab-delimited, with peptides in a column named
  "%s" and comma-delimited lists of protein IDs in a column
  named "%s".  Specifying "-" will read from stdin.

  The output is a series of pairs of files with names of the form
  <root>.<int>.gvz and <root>.<int>.png, where <root> is given on the
  command line and <int> is an ascending integer.  The gvz file
  contains a graphviz description of one component of the graph, and
  the png file contains a picture of the graph.  An HTML file named
  <root>.html is also created, showing all of the PNG images.

  Options:

    --min-nodes <int>   Skip components with fewer than <int> nodes.
    --min-peptides <int>  Skip components with fewer than <int> peptides.
    --min-proteins <int>  Skip components with fewer than <int> proteins.
    --protein-mapping <file>  Convert protein names.

""" % (PEPTIDE_COLUMN, PROTEIN_COLUMN)

# Define the graph data structure.
# N.B. These are global variables.
edges = {}    # Key = peptide sequence or protein ID, value = list of 
              #       peptide sequences or protein IDs
isPrinted = {}# Key = peptide sequence, protein ID or (peptide, protein) pair.
              # Value = True
  
###############################################################################
# Add one edge to the graph.  Each edge is represented twice (once for each
# direction) with a node as the key and a list of nodes as the value.
def addEdge(node1, node2):
  global edges
  
  if (edges.has_key(node1)):
    edges[node1].append(node2)
  else:
    edges[node1] = [node2]

  if (edges.has_key(node2)):
    edges[node2].append(node1)
  else:
    edges[node2] = [node1]

###############################################################################
def printConnectedPeptides(startProtein, graphvizString, stats):
  global edges, isPrinted
  
  if (startProtein not in isPrinted):
    graphvizString = "%s\"%s\";\n" % (graphvizString, startProtein)
    isPrinted[startProtein] = True # Mark this protein as printed.
    stats[0] += 1

    for peptide in edges[startProtein]:
      if ((peptide, startProtein) not in isPrinted):
        graphvizString = "%s\"%s\" -- \"%s\";\n" % (graphvizString, 
                                                    startProtein, peptide)
        isPrinted[(peptide, startProtein)] = True # Mark this edge as printed.
        stats[2] += 1
        graphvizString = printConnectedProteins(peptide, graphvizString, stats)

  return(graphvizString)

###############################################################################
def printConnectedProteins(startPeptide, graphvizString, stats):
  global edges, isPrinted
  
  if (startPeptide not in isPrinted):
    graphvizString = "%s\"%s\";\n" % (graphvizString, startPeptide)
    isPrinted[startPeptide] = True # Mark this peptide as printed.
    stats[1] += 1

    for proteinID in edges[startPeptide]:
      if ((startPeptide, proteinID) not in isPrinted):
        graphvizString = "%s\"%s\" -- \"%s\";\n" % (graphvizString,
                                                    proteinID, startPeptide)
        isPrinted[(startPeptide, proteinID)] = True # Mark this edge as printed.
        stats[2] += 1
        graphvizString = \
            printConnectedPeptides(proteinID, graphvizString, stats)

  return(graphvizString)

#############################################################################
# Run a command with error checking.
def runCommand(command):
  sys.stderr.write("RUN: %s\n" % command)
  try:
    returnCode = subprocess.call(command, shell=True)
    if (returnCode != 0):
      sys.stderr.write("Child was terminated by signal %d\n" % -returnCode)
      sys.exit(1)
  except OSError, e:
    sys.stderr.write("Execution failed: %s\n" % e)
    sys.exit(1)

###############################################################################
# MAIN
###############################################################################

def main():
  global USAGE, PROTEIN_COLUMN, PEPTIDE_COLUMN
  global edges, isPrinted
  
  # Parse the command line.
  minNodes = 1
  minPeptides = 1
  minProteins = 1
  proteinMappingFileName = ""
  sys.argv = sys.argv[1:]
  while (len(sys.argv) > 2):
    nextArg = sys.argv[0]
    sys.argv = sys.argv[1:]
    if (nextArg == "--min-nodes"):
      minNodes = int(sys.argv[0])
      sys.argv = sys.argv[1:]
    elif (nextArg == "--min-peptides"):
      minPeptides = int(sys.argv[0])
      sys.argv = sys.argv[1:]
    elif (nextArg == "--min-proteins"):
      minProteins = int(sys.argv[0])
      sys.argv = sys.argv[1:]
    elif (nextArg == "--protein-mapping"):
      proteinMappingFileName = sys.argv[0]
      sys.argv = sys.argv[1:]
    else:
      sys.stderr.write("Invalid option (%s).\n" % nextArg)
      sys.exit(1)
  if (len(sys.argv) != 2):
    sys.stderr.write(USAGE)
    sys.exit(1)
  inputFileName = sys.argv[0]
  outputFileRoot = sys.argv[1]
  
  # If provided, read the protein mapping into a dictionary.
  proteinMapping = {} # Key = old name, value = new name
  if (proteinMappingFileName != ""):
    proteinMappingFile = open(proteinMappingFileName, "r")
    for line in proteinMappingFile:
      (oldName, newName) = line.rstrip().split()
      proteinMapping[oldName] = newName
    proteinMappingFile.close()
  
  # Read the header line and identify the target columns.
  if (inputFileName == "-"):
    inputFile = sys.stdin
  else:
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

  # Nodes of the graph.
  peptides = {} # Key = peptide sequence, value = True
  proteins = {} # Key = protein ID, value = True
  # N.B. Edges are stored in a global variable.

  # Read the graph from the input file.
  lineNum = 0
  numEdges = 0
  for line in inputFile:
    line = line.rstrip()
    words = line.split("\t")
    peptideSequence = words[peptideColumn]
    proteinIDs = words[proteinColumn].split(",")
    peptides[peptideSequence] = False
    for proteinID in proteinIDs:
      proteinID = proteinID.split("(")[0] # Get rid of (<int>) on each ID.
      if (proteinID in proteinMapping):
        proteinID = proteinMapping[proteinID]
      proteins[proteinID] = False
      addEdge(peptideSequence, proteinID)
      numEdges += 1
    lineNum += 1
  inputFile.close()
  sys.stderr.write("Read %d peptides, %d proteins and %d edges from %s.\n"
                   % (len(peptides), len(proteins), numEdges, inputFileName))
  
  # Initialize the HTML
  htmlFileName = "%s.html" % outputFileRoot
  htmlFile = open(htmlFileName, "w")
  htmlFile.write("<html><body>\n")
  
  graphNumber = 0
  for startProtein in proteins:
    if (startProtein not in isPrinted):
      sys.stderr.write("Starting with %s.\n" % startProtein)
  
      # Initialize a counter of number of proteins, peptides, edges.
      stats = [0, 0, 0]
  
      # Create the graphviz-formated graph.
      graphvizString = "graph G {\n"
      graphvizString = printConnectedPeptides(startProtein, 
                                              graphvizString, stats)
      graphvizString = "%s}\n" % graphvizString
  
      # Is this component big enough?
      if ((stats[0] + stats[1] >= minNodes) and
          (stats[0] >= minProteins) and
          (stats[1] >= minPeptides)):
  
        # Print the graphviz file.
        graphFileName = "%s.%d.gvz" % (outputFileRoot, graphNumber)
        graphFile = open(graphFileName, "w")
        graphFile.write(graphvizString)
        graphFile.close()
  
        # Create the PNG file.
        pngFileName = "%s.%d.png" % (outputFileRoot, graphNumber)
        runCommand("dot -Tpng %s > %s" % (graphFileName, pngFileName))
  
        # Add it to the HTML page.
        htmlFile.write("<img src=\"%s\"><br>\n" % pngFileName)
        message = "%d: %d proteins, %d peptides, %d edges.\n" % \
           (graphNumber, stats[0], stats[1], stats[2])
        htmlFile.write("%s<hr></hr>\n" % message)
  
        sys.stderr.write(message)
        graphNumber += 1
  
      else:
        sys.stderr.write("Skipping component with %d proteins and %d peptides.\n"
                         % (stats[0], stats[1]))
  
  sys.stderr.write("Printed %d graphs.\n" % graphNumber)
  htmlFile.write("</body></html>\n")
  htmlFile.close()

if __name__ == "__main__":
  main()
  
