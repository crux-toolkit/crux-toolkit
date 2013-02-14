#!/usr/bin/env python
# Benjamin Diament
#
# Program to recompute the XCorr computation for output from search.
#
# Run search with --debug_peaks to show peaks from preprocessed
# spectra; save results to an output file to be supplied on the
# command line here.

import aamass, re, sys

aamass.MONOISOTOPIC['C'] += 57.0

aa_pat = re.compile("([A-Z])(?:[[][+]([0-9.]+)[]])?")

def GetAAMasses(peptide):
  return [aamass.MONOISOTOPIC[aa] + float("0"+additional) for aa, additional in aa_pat.findall(peptide)]


bin_width = 1.0005079
co = 27.9949

def AddToTheor(theor, peak_pos, val):
  float_bin = peak_pos / bin_width + 0.5
  bin = int(float_bin)
  if bin + 1 - float_bin < 0.0001:
    print "WARNING"
  if val == 2:
    bins = [(bin, val)]
  elif val == 10:
    bins = [(bin-1, 5), (bin, 10), (bin+1, 5)]
  else:
    raise Exception("val = %s should have been 2 or 10" % val)
  for b, v in bins:
    theor[b] = max(theor.get(b, 0), v)


def AddToTheorWithCharge(theor, peak_pos, val, charge):
  h = aamass.MONOISOTOPIC['H_ATOM']
  AddToTheor(theor, peak_pos + h, val)
  if charge == 3:
    AddToTheor(theor, peak_pos/2.0 + h, val)


def TheorSpec(aa_masses, charge):
  theor = {}
  total = 0.0
  for aa in aa_masses[:-1]:
    total += aa
    AddToTheorWithCharge(theor, total, 10, charge)
    AddToTheorWithCharge(theor, total - aamass.MONOISOTOPIC['H2O'], 2, charge)
    AddToTheorWithCharge(theor, total - aamass.MONOISOTOPIC['NH3'], 2, charge)
    AddToTheorWithCharge(theor, total - co, 2, charge)
  total = aamass.MONOISOTOPIC['H2O']
  for aa in reversed(aa_masses[1:]):
    total += aa
    AddToTheorWithCharge(theor, total, 10, charge)
    AddToTheorWithCharge(theor, total - aamass.MONOISOTOPIC['NH3'], 2, charge)
  return theor


def ConfirmXcorr(charge, xcorr, peptide):
  theor = TheorSpec(GetAAMasses(peptide), charge)
  #print "theor", sorted(theor.keys())
  dotprod = 0.0
  for bin, val in theor.items():
    dotprod += val * peaks_[bin]
  dotprod /= 2000.0

  #print dotprod
  if abs(xcorr-dotprod) > 0.0001:
    raise Exception("Recomputed XCorr = %g" % dotprod)


results_file = open(sys.argv[1])

peaks_dirty = True
count = 0
for line in results_file:
  if line[:5] == "peaks":
    if peaks_dirty:
      peaks_ = [0.0] * 10000
      peaks_dirty = False
    exec line
  else:
    count += 1
    peaks_dirty = True
    fields = line.split()
    a, b, charge, xcorr, peptide = fields
    charge = int(charge)
    xcorr = float(xcorr)
    print count, line[:-1]
    ConfirmXcorr(charge, xcorr, peptide)
