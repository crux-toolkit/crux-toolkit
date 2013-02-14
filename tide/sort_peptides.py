#!/usr/bin/env python

import sys, records
from protoobj import peptides_pb2

peptides = list(records.Records(sys.stdin, peptides_pb2.Peptide))
mass_list = [(peptide.mass, pos) for peptide, pos in zip(peptides, range(len(peptides)))]
mass_list.sort()

writer = records.RecordWriter(sys.stdout)
peptide = peptides_pb2.Peptide()
id = 0
for mass, pos in mass_list:
  peptide.CopyFrom(peptides[pos])
  peptide.id = id
  writer.write(peptide)
  id += 1
