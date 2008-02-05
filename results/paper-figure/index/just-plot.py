import plotFuncs, pylab
plotFuncs.plotLogLogFromFile("crux-no-index-human.fasta", style="k-o",
    label="Crux (w/o index)")
plotFuncs.plotLogLogFromFile("sequest-human.fasta", style="m-s",
    label="Sequest (w/o index)")
plotFuncs.plotLogLogFromFile("crux-fast-human.fasta", style="k--o", 
    label="Crux (w/ index)")
pylab.legend(loc="lower right")
pylab.xlabel("Mass window (Da)", size=20)
pylab.xlim(0.05, 5)
pylab.ylim(0.1, 5000.0)
pylab.xticks(size=20)
pylab.yticks(size=20)
pylab.ylabel("Runtime for 100 spectra (s)", size=20)
pylab.savefig("indexing-human.eps")
pylab.savefig("indexing-human.png")
pylab.clf()

plotFuncs.plotLogLogFromFile("crux-no-index-yeast.fasta", style="k-o",
    label="Crux (w/o index)")
plotFuncs.plotLogLogFromFile("sequest-yeast.fasta", style="m-s",
    label="Sequest (w/o index)")
plotFuncs.plotLogLogFromFile("crux-fast-yeast.fasta", style="k--o", 
    label="Crux (w/ index)")
pylab.legend(loc="lower right")
pylab.xlabel("Mass window (Da)", size=20)
pylab.xlim(0.05, 5)
pylab.ylim(0.1, 5000.0)
pylab.xticks(size=20)
pylab.yticks(size=20)
pylab.ylabel("Runtime for 100 spectra (s)", size=20)
pylab.savefig("indexing-yeast.eps")
pylab.savefig("indexing-yeast.png")
