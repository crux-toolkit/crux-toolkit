import MatchCollection
import plotFDR
import pylab
from Dataset import *
import copy
# import pdb

forwardSequest = "forward/spectra_random_1000.sqt"
randomSequest  = "random/spectra_random_1000.sqt"

#------- SEQUEST analysis -----------------------------------------
forward = MatchCollection.SqtMatchCollection()
forward.parseFromFilename(forwardSequest)
forward.saveSpace = True
forward.allMatches()

random = MatchCollection.SqtMatchCollection()
random.saveSpace = True
random.parseFromFilename(randomSequest)
random.allMatches()

sqtPsmColl = SqtPsmInfoCollection()
sqtPsmColl.addSqtCollection(forward, real=True)
sqtPsmColl.addSqtCollection(random, real=False)

positives = sqtPsmColl.getScores(real=True)
negatives = sqtPsmColl.getScores(real=False)
xs, ys = plotFDR.calcQ(positives, negatives)

fh = open("pq.xy", "w")
for idx in range(len(xs)):
  fh.write("%.6f\t%.6f\n" % (xs[idx], ys[idx]))
fh.close()
