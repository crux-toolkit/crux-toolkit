min_mass = 400.00
max_mass = 1400.0
step = 50.0 
for i in range(int(min_mass),int(max_mass),int(step)):
  print "%.2f\t%.2f" % (i, i+step)
