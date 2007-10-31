import copy
import math
import os
 
filename = "match_analysis.txt"

xs = []
fh = open(filename)
for line in fh:
  fields = line.rstrip("\n").split()
  value = fields[4]
  if (value == "inf"):
    xs.append(0.0)
  else:
    xs.append(pow(math.e, -float(value)))
  
xs.sort()
ys = map(float, range(len(xs)+1))

fh = open("pq.xy", "w")
for idx in range(len(xs)):
  fh.write("%.6f\t%.6f\n" % (xs[idx], ys[idx]))
fh.close()
