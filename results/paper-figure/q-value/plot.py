import plotFuncs
import pylab

plotFuncs.plotXYsFromFile("Pq_percolator/pq", label="Percolator", style="b-")
plotFuncs.plotXYsFromFile("Pq_percolator-control/pq", label="Percolator control", style="b--")
plotFuncs.plotXYsFromFile("Pq_onthefly/pq", label="On the fly null", style="k-")
plotFuncs.plotXYsFromFile("Pq_static/pq", label="Static null", style="g-")
plotFuncs.plotXYsFromFile("Pq_sequest/pq", label="Sequest", style="m-")
# plotFuncs.plotXYsFromFile("Pq_percolator-advanced/pq", label="Percolator extended", style="o-")
pylab.legend()
pylab.xlabel("q-value")
pylab.xlim(0.0,0.10)
pylab.ylabel("positives")
pylab.ylim(0.0,1000.0)
pylab.savefig("q-value.eps")
pylab.savefig("q-value.png")
