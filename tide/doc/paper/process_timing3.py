#!/usr/bin/env python

import re

def Fields(s):
  prog = re.search(r'crux|tide|sequest_orig|sequest_recent', s).group()
  org = re.search(r'worm|yeast', s).group()
  window = 0.25 if re.search(r'smallwindow', s) else 3.0
  digest = 'partial' if re.search(r'partial', s) else 'full'
  #spectra = int(re.search(r'.*?/.*?(10*)', s).group(1))
  spectra = 100
  return (prog, org, window, digest, spectra)

def Secs(t):
  m = re.match(r'([0-9]+)m([0-9.]+)s', t)
  return 60*int(m.group(1)) + float(m.group(2))

def Summary(list):
  return (sum(list)/len(list), len(list), list)

def AddToStats(line, stats):
  broken = line.split()
  fields = Fields(broken[0])
  secs = Secs(broken[1])
  stats[fields] = stats.get(fields, []) + [secs]

def GetStats(lines):
  stats = {}
  for line in lines:
    AddToStats(line, stats)
  results = []
  for fields, times in stats.items():
    (prog, org, window, digest, num_spectra) = fields
    num_runs = len(times)
    avg_secs = sum(times)/num_runs
    results += [(prog, org, window, digest, num_spectra/avg_secs, num_spectra, avg_secs, num_runs, times)]
  results.sort()
  return results

stats = GetStats(open('summary').readlines())

# This replaces all floats in stats with a two-digit (%.2f) corresponding respresentation, and gives a string
two_digits = re.sub(r'[0-9]+[.][0-9]+', lambda s: "%.2f" % float(s.group()), "%s" % stats)

print re.sub(r'[)], ', '),\n', two_digits)

#open('timing_results4','w').write("%s" % y)
