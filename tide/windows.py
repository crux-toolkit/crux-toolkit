#!/usr/bin/env python

x=open('bar').readlines()
w=[line.split() for line in x]
x=[(int(a), int(float(b)*100)) for a,b in w]

window = 600

for low in xrange(len(x)):
  low_count, low_mass = x[low]
  for high_count, high_mass in x[low:]:
    if high_mass - low_mass > window:
      print high_count - low_count, low_mass
      break
