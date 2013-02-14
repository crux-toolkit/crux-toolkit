#!/usr/bin/env python

# Read a file of records of a type given on the command line.
# Usage: readpb.py <filename> <protobuf.Type>
#        where <filename> is the file of records,
#        protobuf is the protocol buffer definition file,
#        and Type is the individual message type of the records.
#        e.g.
#        readpb.py tmp.pepix peptides.Peptide

import sys, records

module, pb_type = sys.argv[2].split('.')
module += "_pb2"

__import__("protoobj", globals(), locals(), [module], -1)
module = sys.modules["protoobj." + module]
pb_type = module.__dict__[pb_type]

header, recs = records.HeadedRecords(open(sys.argv[1], "rb"), pb_type)

print header
for r in recs:
  print r
