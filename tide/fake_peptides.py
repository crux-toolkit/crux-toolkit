#!/usr/bin/env python
# Benjamin Diament

import sys, params, records
from protoobj import header_pb2, raw_proteins_pb2, peptides_pb2

PARAMS, ARGS = params.ParseParameters()

INFILE = PARAMS.get('file', None)

def Write(outfile):
  header = header_pb2.Header()
  header.file_type = header_pb2.Header.PEPTIDES
  peptides_header = header.peptides_header
  peptides_header.has_peaks = False
  
  writer = records.HeadedRecordWriter(outfile, header)
  for m in ARGS:
    peptide = peptides_pb2.Peptide()
    peptide.mass = float(m)
    writer.write(peptide)    

def Read(infile):
  header, reader = records.HeadedRecords(infile, peptides_pb2.Peptide)
  print "====== HEADER: ======"
  print header
  print "====================="
  print
  for record in reader:
    print record


def main():
  if (INFILE):
    Read(open(INFILE, "rb"))
  else:
    Write(sys.stdout)


main()
