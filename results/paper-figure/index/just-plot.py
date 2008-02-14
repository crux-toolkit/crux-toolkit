import plotFuncs, pylab
plotFuncs.plotLogLogFromFile("sequest-human.fasta", style="m-x",
    label="Sequest (w/o index)", yCorrectionFactor=0.01)
plotFuncs.plotLogLogFromFile("crux-no-index-human.fasta", style="k-o",
    label="Crux (w/o index)", yCorrectionFactor=0.01)
plotFuncs.plotLogLogFromFile("crux-fast-human.fasta", style="k--o", 
    label="Crux (w/ index)", yCorrectionFactor=0.01)
pylab.legend(loc="lower right")
pylab.xlabel("Mass window (Da)", size=20)
pylab.xlim(0.05, 5)
pylab.ylim(0.001, 50.0)
pylab.xticks(size=20)
pylab.yticks(size=20)
pylab.ylabel("Runtime per spectrum (s)", size=20)
pylab.savefig("indexing-human.eps")
pylab.savefig("indexing-human.png")
pylab.clf()

plotFuncs.plotLogLogFromFile("sequest-yeast.fasta", style="m-x",
    label="Sequest (w/o index)", yCorrectionFactor=0.01)
plotFuncs.plotLogLogFromFile("crux-no-index-yeast.fasta", style="k-o",
    label="Crux (w/o index)", yCorrectionFactor=0.01)
plotFuncs.plotLogLogFromFile("crux-fast-yeast.fasta", style="k--o", 
    label="Crux (w/ index)", yCorrectionFactor=0.01)
pylab.legend(loc="lower right")
pylab.xlabel("Mass window (Da)", size=20)
pylab.xlim(0.05, 5)
pylab.ylim(0.001, 50.0)
pylab.xticks(size=20)
pylab.yticks(size=20)
pylab.ylabel("Runtime per spectrum (s)", size=20)
pylab.savefig("indexing-yeast.eps")
pylab.savefig("indexing-yeast.png")
