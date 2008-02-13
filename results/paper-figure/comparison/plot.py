import pylab


# '../index/sequest-full-xcorr-yeast.fasta.xy', '../index/crux-no-index-yeast.fasta.xy', '../index/crux-fast-yeast.fasta.xy', 



def parseXY(file):
  fh = open(file)
  times = []
  for line in fh:
    times.append(float(line.split()[1]))
  fh.close()
  return times

xs = [0.1, 1.0, 3.0]
def plotRatio(indexed, noindex, style, label):
  ratios = []
  for idx in range(len(indexed)):
    ratios.append(noindex[idx] / indexed[idx])
  pylab.loglog(xs, ratios, style, label=label)

# '../index/sequest-full-xcorr-human.fasta.xy', 


def parseTime(filename):
  masses = ['0.1', '1', '3']
  
  times = []
  for mass in masses:
    file = filename % mass
    fh = open(file)
    for line in fh:
      if line.startswith("real"):
        secondsTime = float(line.split("m")[1].split("s")[0])
        minutesTime = float(line.split()[1].split("m")[0])
        times.append(secondsTime + 60.0 * minutesTime)
        break
    fh.close()
  return times 


cruxWindowsNone = parseTime('../turbo-no-missed-human/%s.crux.no.time')
cruxWindowsIndexed = parseTime('../turbo-no-missed-human/%s.crux.time')

plotRatio(cruxWindowsIndexed, cruxWindowsNone, style='k-o', label = "Crux on Windows")

cruxLinuxNone    = parseXY('../index/crux-no-index-human.fasta.xy')
cruxLinuxIndexed = parseXY('../index/crux-fast-human.fasta.xy')

plotRatio(cruxLinuxIndexed, cruxLinuxNone, style='k--o', label="Crux on Linux")

sequestWindowsNone = parseTime('../turbo-no-missed-human/%s.sequest.no.time')
sequestWindowsIndexed = parseTime('../turbo-no-missed-human/%s.sequest.time')
plotRatio(sequestWindowsIndexed, sequestWindowsNone, style='m-o', label = "Sequest on Windows")



pylab.xlim(0.05, 5)
pylab.ylim(0.1, 5000.0)
pylab.xticks(size=20)
pylab.yticks(size=20)
pylab.ylabel("Runtime (No index / index)", size=20)
pylab.xlabel("Mass window (Da)", size=20)
pylab.legend()
pylab.savefig("ratio.png")
pylab.savefig("ratio.eps")
pylab.show()



