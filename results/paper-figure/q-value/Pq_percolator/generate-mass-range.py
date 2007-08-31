min_mass = 0.00
max_mass = 7200.0
step = 200.0 
for i in range(int(min_mass),int(max_mass),int(step)):
  print "%.2f\t%.2f" % (i, i+step)
