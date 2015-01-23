#!/usr/bin/env python
import sys
import random

# Minimum tryptic peptide length.
MIN_PEPTIDE = 7

# Maximum number of times to attempt to shuffle each peptide.
NUM_SHUFFLES = 10

USAGE = """USAGE: crux-generate-decoys.py <file> <output>

  This program takes as input a protein FASTA file and produces as
  output four files:

    - <output>.peptide.target.txt: A list of tryptic peptides (using
      the simple KR cleavage rule) derived from the proteins.

    - <output>.peptide.decoy.txt: A matched list of shuffled tryptic
      peptides, where the N-term and C-term amino acids remain in
      place.

    - <output>.protein.decoy.fa: A decoy protein FASTA file that
      matches the input file, but with shuffled peptides in place of
      the original peptides.

    - <output>.log.txt: A log file with various information about the
      run.

  The program attempts to ensure that there are no duplicate peptides
  in the union of the target and decoy peptide lists.  Peptides for
  which a decoy is not successfully created (e.g., homopolymers) are
  indicated in the log file. Peptides shorter than %d amino acids are
  not shuffled.

""" % MIN_PEPTIDE

##############################################################################
def read_fasta_sequence (fasta_file):

  # Read 1 byte.  
  first_char = fasta_file.read(1)

  # If it's empty, we're done.
  if (first_char == ""):
    return(["", ""])
  # If it's ">", then this is the first sequence in the file.
  elif (first_char == ">"):
    line = ""
  else:
    line = first_char

  # Read the rest of the header line.
  line = line + fasta_file.readline()

  # Get the rest of the ID.
  words = line.split()
  if (len(words) == 0):
    sys.stderr.write("No words in header line (%s)\n" % line)
    sys.exit(1)
  id = words[0]
      
  # Read the sequence, through the next ">".
  first_char = fasta_file.read(1)
  sequence = ""
  while ((first_char != ">") and (first_char != "")):
    if (first_char != "\n"): # Handle blank lines.
      line = fasta_file.readline()
      sequence = sequence + first_char + line
    first_char = fasta_file.read(1)

  # Remove EOLs.
  clean_sequence = ""
  for letter in sequence:
    if (letter != "\n"):
      clean_sequence = clean_sequence + letter
  sequence = clean_sequence

  # Remove spaces.
  clean_sequence = ""
  for letter in sequence:
    if (letter != " "):
      clean_sequence = clean_sequence + letter
      
  return([id, sequence.upper()])

##############################################################################
# Write a message both to stderr and a globally defined log file.
def log(message):

  sys.stdout.write(message)
  logFile.write(message)

##############################################################################
# Convert an amino acid sequence into a list of tryptic peptides,
# cleaving at every K or R (irrespective of P).
def cleaveTryptically(sequence):

  returnValue = []   

  peptideStart = 0
  for position in range(0, len(sequence)):
    if (sequence[position] == "K") or (sequence[position] == "R"):
      returnValue.append(sequence[peptideStart:position+1])
      peptideStart = position+1
  returnValue.append(sequence[peptideStart:position+1])

  #log("%s -> %s\n" % (sequence, "|".join(returnValue)))
  return(returnValue)

##############################################################################
### MAIN
##############################################################################

# Parse the command line.
if (len(sys.argv) != 3):
  sys.stderr.write(USAGE)
  sys.exit(1)
inputFileName = sys.argv[1]
root = sys.argv[2]

# Open the log file for output.
logFileName = "%s.log.txt" % root
logFile = open(logFileName, "w")

# Ordered list of protein IDs.
proteinIDs = []

proteinSeqs = {} # Key = ID, value = list of tryptic peptides
targetPeptides = {} # Key = peptide, value = True

# Read the file sequence by sequence.
inputFile = open(inputFileName, "r")
[id, sequence] = read_fasta_sequence(inputFile)
while (id != ""):

  proteinIDs.append(id)
  proteinSeqs[id] = cleaveTryptically(sequence)

  # Read the target peptides into a hash.
  for peptide in proteinSeqs[id]:
    targetPeptides[peptide] = True

  # Read the next sequence.
  [id, sequence] = read_fasta_sequence(inputFile)

inputFile.close()
log("Read %d proteins and %d peptides from %s.\n" % 
    (len(proteinIDs), len(targetPeptides), inputFileName))




# Open the peptide output lists.
targetPeptideFile = open("%s.peptide.target.txt" % root, "w")
decoyPeptideFile = open("%s.peptide.decoy.txt" % root, "w")

# Create the decoys.
decoyPeptides = {} # Key = peptide, value = True
targetToDecoy = {} # Key = target peptide, value = decoy peptide
for targetPeptide in targetPeptides.keys():

  # Don't bother with short peptides.
  if (len(targetPeptide) <= MIN_PEPTIDE):
    continue

  success = False
  for shuffle in range(0, NUM_SHUFFLES):
    decoyList = list(targetPeptide)[1:-1] # Don't shuffle terminal AAs.
    random.shuffle(decoyList)
    decoyList.insert(0,targetPeptide[0])
    decoyList.append(targetPeptide[-1])
    decoyPeptide = ''.join(decoyList)
    if ((not targetPeptides.has_key(decoyPeptide)) and
        (not decoyPeptides.has_key(decoyPeptide))):
      decoyPeptides[decoyPeptide] = True
      success = True
      break

  if not success:
    if targetPeptides.has_key(decoyPeptide):
      log("overlap decoy %s\n" % decoyPeptide)
    if decoyPeptides.has_key(decoyPeptide):
      log("duplicate decoy %s\n" % decoyPeptide)
      continue

  targetPeptideFile.write("%s\n" % targetPeptide)
  decoyPeptideFile.write("%s\n" % decoyPeptide)
  targetToDecoy[targetPeptide] = decoyPeptide

log("Printed %d peptides.\n" % len(targetToDecoy))




# Print the decoy proteins.
decoyProteinFileName = "%s.protein.decoy.fa" % root
decoyProteinFile = open(decoyProteinFileName, "w")
for proteinID in proteinIDs:
  decoyProteinFile.write(">%s\n" % proteinID)
  for targetPeptide in proteinSeqs[proteinID]:
    if targetPeptide in targetToDecoy:
      decoyProteinFile.write("%s" % targetToDecoy[targetPeptide])
    else:
      decoyProteinFile.write("%s" % targetPeptide)
  decoyProteinFile.write("\n")
decoyProteinFile.close()
log("Printed %d decoy proteins to %s.\n" % 
    (len(proteinIDs), decoyProteinFileName))


logFile.close()
