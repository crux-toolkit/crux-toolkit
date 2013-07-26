#! /usr/bin/python

import sys, params, records
from protoobj import raw_proteins_pb2, peptides_pb2, header_pb2

PARAMS, ARGS = params.ParseParameters()

MODS_FILE = PARAMS['mods']
PROTEIN_FILE = PARAMS['proteins']
PEPTIDE_FILE = PARAMS['peptides']
OUT_FILE = PARAMS['out_peptides']
ID = PARAMS.get('id', None)

modheader, mods = records.HeadedRecords(open(MODS_FILE, 'rb'), header_pb2.ModTable)
mods = [i for i in mods][0]
protheader, proteins = records.HeadedRecords(open(PROTEIN_FILE, 'rb'), raw_proteins_pb2.Protein)
proteins = list(proteins)
pepheader, peptides = records.HeadedRecords(open(PEPTIDE_FILE, 'rb'), peptides_pb2.Peptide)

def UseOnce(d, val):
  assert(not d.has_key(val))
  d[val] = None

MAX_COUNTS = []
UNIQUE_DELTAS = []
MOD_DICT = {}
LOG_UNIQUE_DELTAS = 0

def IntLog(x):
  if x <= 1: return x
  result = 0
  x -= 1
  while x > 0:
    x >>= 1
    result += 1
  return result


def MakeModDict(v_mods):
  # Check no repeats
  d = {}
  for i in v_mods:
    for j in i.amino_acids:
      UseOnce(d, (j, i.delta))
  global MAX_COUNTS
  MAX_COUNTS = [i.max_count for i in v_mods]
  global UNIQUE_DELTAS
  deltas = {}
  for i in v_mods:
    deltas[i.delta] = None
  UNIQUE_DELTAS = sorted(deltas.keys())
  for i, j in zip(UNIQUE_DELTAS, range(len(UNIQUE_DELTAS))):
    deltas[i] = j
  # Dictionary maps each amino_acid to a list of (delta_index, counter_index)
  global MOD_DICT
  counter_index = 0
  for i in v_mods:
    for aa in i.amino_acids:
      MOD_DICT[aa] = MOD_DICT.get(aa, []) + [(deltas[i.delta], counter_index)]
    counter_index += 1
  global LOG_UNIQUE_DELTAS
  LOG_UNIQUE_DELTAS = IntLog(len(UNIQUE_DELTAS))


MakeModDict(mods.variable_mod)
mods.ClearField("unique_deltas")
for i in UNIQUE_DELTAS:
  mods.unique_deltas.append(i)
pepheader.peptides_header.mods.CopyFrom(mods)
out_peptides = records.HeadedRecordWriter(open(OUT_FILE, 'wb'), pepheader)


def WriteAllModVersions(writer, peptide, pep_seq, pos, counters):
  writer.write(peptide)
  while pos < len(pep_seq):
    possibles = MOD_DICT.get(pep_seq[pos], [])
    for delta_index, counter_index in possibles:
      if counters[counter_index] < MAX_COUNTS[counter_index]:
        counters[counter_index] += 1
        old_mass = peptide.mass
        peptide.modifications.append((pos << LOG_UNIQUE_DELTAS) + delta_index)
        peptide.mass += UNIQUE_DELTAS[delta_index]
        WriteAllModVersions(writer, peptide, pep_seq, pos+1, counters)
        del peptide.modifications[-1]
        peptide.mass = old_mass
        counters[counter_index] -= 1
    pos += 1

peptides = [p for p in peptides]
if ID:
  peptides = [peptides[int(ID)]]
for peptide in peptides:
  location = peptide.location[0]
  protein = proteins[location.protein_id]
  pos = location.pos
  length = peptide.length
  assert(length > 0)
  pep_seq = protein.residues[pos:pos+length]
  assert(pep_seq)

  WriteAllModVersions(out_peptides, peptide, pep_seq, 0, [0]*len(MAX_COUNTS))
