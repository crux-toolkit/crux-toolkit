import plotFuncs
import pylab

plotFuncs.plotXYsFromFile("Pq_percolator/pq", label="Percolator", style="b-")
plotFuncs.plotXYsFromFile("Pq_pvalue/pq", label="p-value", style="m-")
plotFuncs.plotXYsFromFile("Pq_onthefly/pq", label="On the fly null", style="g-")
plotFuncs.plotXYsFromFile("Pq_static/pq", label="Static null", style="r--")
plotFuncs.plotXYsFromFile("Pq_sequest/pq", label="Sequest", style="y--")
pylab.legend(loc="lower right")
pylab.xlabel("q-value")
pylab.xlim(0.0,0.10)
pylab.ylabel("positives")
pylab.ylim(0.0,8000.0)
pylab.savefig("q-value.eps")
pylab.savefig("q-value.png")
