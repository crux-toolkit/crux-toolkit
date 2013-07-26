#!/usr/bin/env python
# Benjamin Diament

# script to convert a fasta file of proteins from stdin to a
# set of records of type raw_proteins.proto format on stdout

# TODO: We are ignoring non-residues -- simply eliminating them from
# the input. It might be more appropriate to divide sequences that
# have odd characters into two separate proteins.

import sys, records, re
from protoobj import raw_proteins_pb2
g_writer = records.RecordWriter(sys.stdout)
g_id = 0


def Add(name, residues):
  if not residues:
    return
  global g_id
  protein = raw_proteins_pb2.Protein()
  protein.name = name
  protein.id = g_id
  g_id += 1
  protein.residues = residues
  g_writer.write(protein)


def ParseFasta():
  # Proteins are separated by comment lines begining with a '>'. 
  # We collect lines containing residues until a comment or EOF, at which
  # point we call Add().
  # Blank lines are ignored.
  non_residues = re.compile('[^ACDEFGHIKLMNPQRSTVWY]')
  residues = []
  name = ""
  for line in sys.stdin:
    line = line.strip()
    if len(line) > 0:
      if line[0] == '>':
        Add(name, ''.join(residues))
        # The name of the protein is the first word after '>'
        name = line[1:].split()[0] 
        residues = []
      else:
        residues += non_residues.split(line)
  Add(name, ''.join(residues))


ParseFasta()
