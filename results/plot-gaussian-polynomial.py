from plotFuncs import *
import IdMap
from pylab import *
LABEL_FONT = 40
TITLE_FONT = 50
TIC_FONT = 20
FDR_THRESHOLDS = [ i / 200.0 for i in range(20) ]
UNFILTERED = 240
RT_THRESHOLDS = range(10,120,5) + [UNFILTERED]
styles = ['k-','k--','k-.'] 
 

def process_file(file):
  id = IdMap.file2id[file]

  figure(figsize=(11,8.5), dpi=300, facecolor='w', edgecolor='k')
  title(id, size=TITLE_FONT)
  xlabel("FDR (%)", size=LABEL_FONT)
  ylabel("TP (1000s)", size=LABEL_FONT)

  hold(True)

  plotXYsFromFileRelabel(
      "gaussian-floor/%s.out/%s-filtered.xy" % (file, id), \
      label="Gaussian", style=styles[0])

  plotXYsFromFileRelabel(
      "polynomial-floor/%s.out/%s-filtered.xy" % (file, id), \
      label="Polynomial", style=styles[1])

  plotXYsFromFileRelabel(
      "gaussian-floor/%s.out/%s-unfiltered.xy" % (file, id), \
      label="Unfiltered", style=styles[2])
 
  legend(loc='lower right')
  # axis([0, 0.10, 0, 1000 * 1.2])

  locs, labels = xticks()
  xtics = [ float(i) * 100.0 for i in locs ] 
  xticks(locs, map(str, xtics), size=TIC_FONT)

  locs, labels = yticks()
  ytics = [ float(i) / 1000.0 for i in locs ] 
  yticks(locs, map(str, ytics), size=TIC_FONT)

  savefig(id + ".eps")
  # savefig(id + ".png", dpi=300)
  # show()

files = ['090306-20cm-yeast-2h-04','090306-40cm-yeast-2h-04', ]

for file in files:
  process_file(file)
  
