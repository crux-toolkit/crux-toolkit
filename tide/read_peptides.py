#! /usr/bin/python

import sys, params, records
from protoobj import raw_proteins_pb2, peptides_pb2

PARAMS, ARGS = params.ParseParameters()

EXTRA_CHECK = PARAMS.has_key('check')
SHOW_LOCATIONS = PARAMS.has_key('locations')
SHOW_FRAGMENTS = PARAMS.has_key('fragments')

PROTEIN_FILE = PARAMS['proteins']
PEPTIDE_FILE = PARAMS['peptides']
AUX_LOCS_FILE = PARAMS.get('aux_locations', None)
SEEK_POS = long(PARAMS.get('seek_pos', 0))

peptideFile = open(PEPTIDE_FILE, 'rb')
pepheader, peptides = records.HeadedRecords(peptideFile, peptides_pb2.Peptide)

if (SEEK_POS > peptideFile.tell()):
  peptideFile.seek(SEEK_POS, 0)

protheader, proteins = records.HeadedRecords(open(PROTEIN_FILE, 'rb'), raw_proteins_pb2.Protein)
proteins = list(proteins)

if AUX_LOCS_FILE:
  auxlocheader, aux_locations = records.HeadedRecords(open(AUX_LOCS_FILE, 'rb'), peptides_pb2.AuxLocation)

# check that protein ids are exactly in serialized order
protein_ids = [protein.id for protein in proteins]
assert(protein_ids == range(len(protein_ids)))

series_lookup = {
  peptides_pb2.Peptide.B: 'B',
  peptides_pb2.Peptide.Y: 'Y',
}

def IntLog(x):
  if x <= 1: return x
  result = 0
  x -= 1
  while x > 0:
    x >>= 1
    result += 1
  return result


UNIQUE_DELTAS = [i for i in pepheader.peptides_header.mods.unique_deltas]
LOG_UNIQUE_DELTAS = IntLog(len(UNIQUE_DELTAS))
UNIQUE_DELTAS_MASK = (1 << LOG_UNIQUE_DELTAS) - 1

def Decode(code):
  return (code >> LOG_UNIQUE_DELTAS, UNIQUE_DELTAS[code & UNIQUE_DELTAS_MASK])

def ModFormat(pep_seq, mod_vec):
  #print pep_seq, mod_vec
  seq = [i for i in pep_seq]
  for pos, delta in mod_vec:
    seq[pos] += "[+%g]" % delta
  return ''.join(seq)


aux_loc_index = -1
for peptide in peptides:
  last_pep_seq = None
  protein = proteins[peptide.first_location.protein_id]
  pos = peptide.first_location.pos
  length = peptide.length
  assert(length > 0)  
  pep_seq = protein.residues[pos:pos+length]
  assert(pep_seq)
  if EXTRA_CHECK and last_pep_seq:
    assert(last_pep_seq == pep_seq)
    last_pep_seq = pep_seq
  mod_vec = [Decode(m) for m in peptide.modifications]
  pep_seq = ModFormat(pep_seq, mod_vec)
  outvec = (peptide.id, peptide.mass, protein.name, pos+1, length, pep_seq)
  print "%s\t%.5f\t%s\t%s\t%s\t%s" % outvec
  if AUX_LOCS_FILE and peptide.HasField("aux_locations_index"):
    aux_loc_index += 1
    assert aux_loc_index == peptide.aux_locations_index, "%s!=%s" % (aux_loc_index, peptide.aux_locations_index)
    aux_loc = aux_locations.next()
    if EXTRA_CHECK:
      # To get the total number of locations, we need to add one because the
      # first_location is included as part of the peptide
      locs = len(aux_loc.location) + 1
      assert(locs > 0)
      if locs > 1:
        print "Peptide %s in %s locations" % (peptide.id, locs)
    if SHOW_LOCATIONS:
      for loc in aux_loc.location:
        protein = proteins[loc.protein_id]
        pos = loc.pos
        mod_vec = [Decode(m) for m in peptide.modifications]
        pep_seq = protein.residues[pos:pos+length]
        outvec = (peptide.id, peptide.mass, protein.name, pos+1, length, pep_seq)
        print "%s\t%.5f\t%s\t%s\t%s\t%s" % outvec
  if SHOW_FRAGMENTS:
    for fragment in peptide.fragment:
      print "\t%s%d:%.5f" % (series_lookup[fragment.series], fragment.length, fragment.mass)
