#!/usr/bin/env python

# script to convert an ms2 file of spectra from stdin to a
# set of records of type spectrum.proto format on stdout

import sys, records
from protoobj import spectrum_pb2
spectrum = None

writer = records.RecordWriter(sys.stdout)

for line in sys.stdin:
  line = line.split()
  if line[0] == 'S':
    if spectrum:
      writer.write(spectrum)
    spectrum = spectrum_pb2.Spectrum()
    assert(line[1] == line[2])
    spectrum.spectrum_number = int(line[1])
    spectrum.precursor_m_z = float(line[3])
  elif line[0] == 'Z':
    spectrum.charge_state.append(int(line[1]))
  elif line[0][0] >= '0' and line[0][0] <= '9':
    peak = spectrum.peak.add()
    peak.m_z = float(line[0])
    peak.intensity = float(line[1])

if spectrum:
  writer.write(spectrum)
