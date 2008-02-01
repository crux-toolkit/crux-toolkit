from pylab import *
windows = ["0.1", "1", "3"]
names = ["crux", "crux.no", "sequest", "sequest.no"]
# names = ["crux"]
styles = ["k-o", "k--o", "m-^", "m--^"]
labels = ["Crux (w/ index)", "Crux (w/o index)", "Sequest (w/ index)", "Sequest (w/o) index"]

idx = 0
for name in names:
  ys = []
  xs = []
  for window in windows:
    filename = "%s.%s.time" % (window, name)
    fh = open(filename)
    for line in fh:
      if line.startswith("real"):
        time = float(line.split("m")[1].split("s")[0])
        ys.append(time)
        xs.append(float(window))
    fh.close()

  loglog(xs, ys, styles[idx], label=labels[idx])
  idx += 1

legend(loc="lower right")
xlabel("Mass window (Da)")
xlim(0.05, 5)
ylim(0.1, 50000)
ylabel("Runtime for 425 spectra (s)")
savefig("indexing-windows.eps")
savefig("indexing-windows.png")
      
