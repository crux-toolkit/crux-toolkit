#! /usr/bin/python

import sys, records
from protoobj import spectrum_pb2

header, spectra = records.HeadedRecords(sys.stdin, spectrum_pb2.Spectrum)

print "=== HEADER ===="
print header
print "=== SPECTRA ==="
for spectrum in spectra:
  print "Spectrum Number: %s  Precursor m/z: %g" % (spectrum.spectrum_number, spectrum.precursor_m_z)
  print "Charge states: %s" % ', '.join([str(i) for i in spectrum.charge_state])
  for peak in spectrum.peak:
    print "%.1f %.1f" % (peak.m_z, peak.intensity)
