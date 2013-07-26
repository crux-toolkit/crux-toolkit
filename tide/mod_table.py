#!/usr/bin/env python
# Benjamin Diament

import sys, os, re, records
from protoobj import header_pb2

def ShowModfile(f):
  header, recs = records.HeadedRecords(f, header_pb2.ModTable)
  mod_table = [mt for mt in recs]
  assert(len(mod_table) == 1)
  print mod_table[0]


def ParseSpec(spec):
  var_mod_spec_re = re.compile(r'\s*([1-9][0-9]*)\s*([ACDEFGHIKLMNPQRSTVWY]*)\s*[+]\s*([0-9].*)\s*')
  static_mod_spec_re = re.compile(r'\s*([ACDEFGHIKLMNPQRSTVWY])\s*[+]\s*([0-9].*)\s*')
  static_mods = []
  var_mods = []
  for i in spec.split(','):
    m = var_mod_spec_re.match(i)
    if m:
      max_ct, aa, delta = m.groups()
      var_mods += [(int(max_ct), aa, float(delta))]
    else:
      m = static_mod_spec_re.match(i)
      if not m:
        raise Exception("Can't parse expression '%s' as either variable or static modification specification." % i)
      aa, delta = m.groups()
      static_mods += [(aa, float(delta))]
  return (static_mods, var_mods)


def MakeModfile(f, spec):
  header = header_pb2.Header()
  header.file_type = header_pb2.Header.MOD_TABLE
  writer = records.HeadedRecordWriter(f, header)
  static_mods, var_mods = spec
  mod_table = header_pb2.ModTable()
  
  for aa, delta in static_mods:
    assert(delta > 0.0)
    assert(len(aa) == 1)
    mod = mod_table.static_mod.add()
    mod.amino_acids = aa
    mod.delta = delta
  for max_ct, aa, delta in var_mods:
    assert(delta > 0.0)
    mod = mod_table.variable_mod.add()
    mod.amino_acids = aa
    mod.delta = delta
    mod.max_count = max_ct

  writer.write(mod_table)


def Main():
  filename = sys.argv[1]
  if os.path.exists(filename):
    f = open(filename, "rb")
    ShowModfile(f)
  else:
    if len(sys.argv) > 2:
        spec = ' '.join(sys.argv[2:])
    else:
        spec = sys.stdin.read()
    spec = ParseSpec(spec)
    f = open(filename, "wb")
    MakeModfile(f, spec)
    f.close()

Main()
