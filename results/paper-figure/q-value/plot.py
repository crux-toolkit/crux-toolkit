import plotFuncs
import pylab

plotFuncs.plotXYsFromFile("Pq_onthefly/pq")
plotFuncs.plotXYsFromFile("Pq_static/pq")
plotFuncs.plotXYsFromFile("Pq_sequest/pq")
plotFuncs.plotXYsFromFile("Pq_percolator/pq")
pylab.savefig("p-value.eps")
pylab.savefig("p-value.png")
pylab.show()
