#!/usr/bin/env python
# Benjamin Diament

# Takes raw_proteins.proto format from standard input, performs a tryptic digest,
# and outputs resulting peptides in peptides.proto format to standard output.

# Example command-line:
# cat <intput file> | make_peptides.py  [options] > <output file>
# Options are:
# --min_mass=<float>
# --max_mass=<float>
# --min_length=<int>
# --max_length=<int>
#    Peptides that fall outside the range of the preceding parameters are
#    discarded.
# --mono 
#    Use monoisotopic masses for aminio acids.
import sys, params, records
from protoobj import raw_proteins_pb2, peptides_pb2

PARAMS, ARGS = params.ParseParameters()

MIN_MASS = float(PARAMS.get('min_mass', 200.0))
MAX_MASS = float(PARAMS.get('max_mass', 7200.0))
MIN_LENGTH = int(PARAMS.get('min_length', 6))
MAX_LENGTH = int(PARAMS.get('max_length', 50))

import aamass

AA_MASS = PARAMS.has_key('mono') and aamass.MONOISOTOPIC or aamass.AVG
AA_MASS['C'] += 57

def GetPeptides(seq):
  # Iterator over tryptic peptides of amino acid sequence seq.
  # Yields tuples (mass, peptide, start pos., length).
  seq_masses = [0] + [AA_MASS[s] for s in seq]
  seqlen = len(seq_masses) - 1

  # seq_masses to get partial sums.
  for i in range(2,seqlen+1):
    seq_masses[i] += seq_masses[i-1]

  # breakpoints is list of tryptic cut points, beginning and end inclusive.
  breakpoints = [0] + [i + 1 for i in range(seqlen-1) if seq[i] in "KR" and seq[i+1] != "P"] + [seqlen]
  bp_len = len(breakpoints)
  for i in range(bp_len - 1):
    #for j in range(i + 1, bp_len):
    start = breakpoints[i]
    #end = breakpoints[j]
    end = breakpoints[i+1]
    if end-start < MIN_LENGTH:
      continue
    elif end-start > MAX_LENGTH:
      continue
      #break
    mass = seq_masses[end] - seq_masses[start] + AA_MASS['H2O']
    if mass < MIN_MASS:
      continue
    if mass > MAX_MASS:
      continue
      #break
    yield mass, seq[start:end], start, end-start


def seq_main():
  # Alternate functionality to get protein sequence from single command line
  # parameter and display peptides. Run thus: 
  #  make_peptides.py --seq YRKLMNN...
  for seq in GetPeptides(ARGS[0]):
    print seq


def main():
  peptides = []
  for protein in records.Records(sys.stdin, raw_proteins_pb2.Protein):
    for mass, pep_seq, start, length in GetPeptides(protein.residues):
      peptides += [(pep_seq, mass, protein.id, start, length)]
  # Sort by sequence so we notice dups.
  peptides.sort()

  writer = records.RecordWriter(sys.stdout)
  peptide = None
  last_pep_seq = ""
  # Build a peptide protocol buffer for each peptide sequence, including all
  # positions. Call write.write(peptide) whenever we've processed all 
  # occurences of a given peptide sequence.
  for pep_seq, mass, id, position, length in peptides:
    if pep_seq != last_pep_seq:
      if peptide:
        writer.write(peptide)
      peptide = peptides_pb2.Peptide()
      # no id assigned until mass-sorting, which happens in another index
      # processing stage
      peptide.mass = mass
      peptide.length = length
      locations = peptide.location
    location = locations.add()
    location.protein_id = id
    location.pos = position
    last_pep_seq = pep_seq
  if peptide:
    writer.write(peptide)


if PARAMS.has_key('seq'):
  seq_main()
else:
  main()
