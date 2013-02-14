#! /usr/bin/python

import sys, params, records
from protoobj import raw_proteins_pb2, peptides_pb2

PARAMS, ARGS = params.ParseParameters()

AUX_LOCS_FILE = PARAMS['aux_locs']
PROTEIN_FILE = PARAMS['proteins']

protheader, proteins = records.HeadedRecords(open(PROTEIN_FILE, 'rb'), raw_proteins_pb2.Protein)
proteins = list(proteins)

aux_loc_header, aux_locs = records.HeadedRecords(open(AUX_LOCS_FILE, 'rb'), peptides_pb2.AuxLocation)

# check that protein ids are exactly in serialized order
protein_ids = [protein.id for protein in proteins]
assert(protein_ids == range(len(protein_ids)))

index = 0
for aux_loc in aux_locs:  
  locations = aux_loc.location
  index = index + 1
  print "aux_loc_index=%d" % (index)
  for location in locations:
    protein = proteins[location.protein_id]
    pos = location.pos    
    print "%s\t%d" % (protein.name, pos+1)
  print "\n"
