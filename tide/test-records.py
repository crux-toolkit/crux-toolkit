#!/usr/bin/env python

# Test reading and writing of HeadedRecords

import sys
import records
from protoobj import header_pb2
from protoobj import raw_proteins_pb2

FILENAME = '/tmp/test_records'

def TestHeader():
  header = header_pb2.Header()
  header.file_type = header_pb2.Header.RAW_PROTEINS
  source = header.source.add()
  source.filename = 'this would be the name of a fasta file'
  source.filetype = 'fasta'
  return header

def TestProteins():
  for id in range(100):
    protein = raw_proteins_pb2.Protein()
    protein.id = id
    protein.name = 'some name'
    protein.residues = 'RESID' + 'LETTERAFTERT' + 'ES'
    yield protein

def TestWrite():
  outfile = open(FILENAME, 'w')
  writer = records.HeadedRecordWriter(outfile, TestHeader())
  for protein in TestProteins():
    writer.write(protein)

def TestRead():
  infile = open(FILENAME)
  header, iterator = records.HeadedRecords(infile, raw_proteins_pb2.Protein)
  assert header == TestHeader()
  for a, b in zip(TestProteins(), iterator):
    assert a == b

TestWrite()
TestRead()
