import plotFuncs, pylab
plotFuncs.plotLogLogFromFile("crux-no-index-human.fasta", style="k-o",
    label="Crux (w/o index)")
plotFuncs.plotLogLogFromFile("sequest-human.fasta", style="m-s",
    label="Sequest")
plotFuncs.plotLogLogFromFile("crux-fast-human.fasta", style="g-^", 
    label="Crux (w/ index)")
pylab.legend(loc="lower right")
pylab.xlabel("Mass window (Da)")
pylab.xlim(0.05, 5)
pylab.ylim(0.1, 50000)
pylab.ylabel("Runtime for 100 spectra (s)")
pylab.savefig("indexing-human.eps")
pylab.savefig("indexing-human.png")
pylab.clf()

plotFuncs.plotLogLogFromFile("crux-no-index-yeast.fasta", style="k-o",
    label="Crux (w/o index)")
plotFuncs.plotLogLogFromFile("sequest-yeast.fasta", style="m-s",
    label="Sequest")
plotFuncs.plotLogLogFromFile("crux-fast-yeast.fasta", style="g-^", 
    label="Crux (w/ index)")
pylab.legend(loc="lower right")
pylab.xlabel("Mass window (Da)")
pylab.xlim(0.05, 5)
pylab.ylim(0.1, 50000)
pylab.ylabel("Runtime for 100 spectra (s)")
pylab.savefig("indexing-yeast.eps")
pylab.savefig("indexing-yeast.png")
