#! /usr/bin/python
# Benjamin Diament
# Testing aid

import sys, records
from protoobj import raw_proteins_pb2

header, proteins = records.HeadedRecords(sys.stdin, raw_proteins_pb2.Protein)
print header

for protein in proteins:
  print "id: %s name: %s" % (protein.id, protein.name)
  print protein.residues
  print
