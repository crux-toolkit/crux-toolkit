import MatchCollection
import plotFDR
import pylab
from Dataset import *
import copy
# import pdb

import os
 
_proc_status = '/proc/%d/status' % os.getpid()

_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}

def _VmB(VmKey):
  '''Private.
  '''
  global _proc_status, _scale
  # get pseudo file  /proc/<pid>/status
  try:
    t = open(_proc_status)
    v = t.read()
    t.close()
  except:
    return 0.0  # non-Linux?
  # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
  i = v.index(VmKey)
  v = v[i:].split(None, 3)  # whitespace
  if len(v) < 3:
    return 0.0  # invalid format?
  # convert Vm value to bytes
  return float(v[1]) * _scale[v[2]]


def memory(since=0.0):
  '''Return memory usage in bytes.
  '''
  return _VmB('VmSize:') - since


def resident(since=0.0):
  '''Return resident memory usage in bytes.
  '''
  return _VmB('VmRSS:') - since


def stacksize(since=0.0):
  '''Return stack size in bytes.
  '''
  return _VmB('VmStk:') - since

print >>sys.stderr, "Current memory usage %s" % memory()

sqts = filter(lambda x: x.endswith(".sqt"), os.listdir("."))
decoys = filter(lambda x: "decoy" in x, sqts)
targets = filter(lambda x: "decoy" not in x, sqts)

negatives = []
positives = []

for decoy in decoys:
  random = MatchCollection.SqtMatchCollection()
  random.saveSpace = True
  random.parseFromFilename(decoy)
  random.allMatches()
  sqtPsmColl = SqtPsmInfoCollection()
  sqtPsmColl.addSqtCollection(random, real=False)
  negatives += sqtPsmColl.getScores(real=False)
  print >>sys.stderr, "Current memory usage %s" % memory()


for target in targets:
  random = MatchCollection.SqtMatchCollection()
  random.saveSpace = True
  random.parseFromFilename(target)
  random.allMatches()
  sqtPsmColl = SqtPsmInfoCollection()
  sqtPsmColl.addSqtCollection(random, real=True)
  positives += sqtPsmColl.getScores(real=True)
  print >>sys.stderr, "Current memory usage %s" % memory()

xs, ys = plotFDR.calcQ(positives, negatives)

fh = open("pq.xy", "w")
for idx in range(len(xs)):
  fh.write("%.6f\t%.6f\n" % (xs[idx], ys[idx]))
fh.close()
