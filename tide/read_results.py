#! /usr/bin/python

import sys, params, records
from protoobj import raw_proteins_pb2, peptides_pb2, results_pb2, spectrum_pb2

PARAMS, ARGS = params.ParseParameters()

RESULTS_FILE = PARAMS['results']
PROTEIN_FILE = PARAMS['proteins']

protheader, proteins = records.HeadedRecords(open(PROTEIN_FILE, 'rb'), raw_proteins_pb2.Protein)
proteins = list(proteins)

resheader, results = records.HeadedRecords(open(RESULTS_FILE, 'rb'), results_pb2.Results)

# check that protein ids are exactly in serialized order
protein_ids = [protein.id for protein in proteins]
assert(protein_ids == range(len(protein_ids)))

def IntLog(x):
  if x <= 1: return x
  result = 0
  x -= 1
  while x > 0:
    x >>= 1
    result += 1
  return result

UNIQUE_DELTAS = [i for i in resheader.results_header.peptides_header.mods.unique_deltas]
LOG_UNIQUE_DELTAS = IntLog(len(UNIQUE_DELTAS))
UNIQUE_DELTAS_MASK = (1 << LOG_UNIQUE_DELTAS) - 1

def Decode(code):
  return (code >> LOG_UNIQUE_DELTAS, UNIQUE_DELTAS[code & UNIQUE_DELTAS_MASK])

def ModFormat(pep_seq, mod_vec):
  #print pep_seq, mod_vec
  seq = [i for i in pep_seq]
  for pos, delta in mod_vec:
    seq[pos] += "[+%.1f]" % delta
  return ''.join(seq)

for result in results:
  spectrum = result.spectrum
  #outvec = (spectrum.spectrum_number, spectrum.precursor_m_z)
  #print "%d\t%.2f" % outvec
  matches = result.matches
  for match in matches:
    peptide = match.peptide
    location = peptide.first_location
    protein = proteins[location.protein_id]
    pos = location.pos
    length = peptide.length
    assert(length > 0)
    pep_seq = protein.residues[pos:pos+length]
    assert(pep_seq)
    mod_vec = [Decode(m) for m in peptide.modifications]
    pep_seq = ModFormat(pep_seq, mod_vec)
    #outvec = (match.charge, match.xcorr, pep_seq, protein.name)
    #print "%d\t%.6f\t%s\t%s" % outvec
    outvec = (spectrum.spectrum_number, spectrum.precursor_m_z, match.charge, match.xcorr, pep_seq)
    print "%d\t%.2f\t%d\t%.6f\t%s" % outvec
