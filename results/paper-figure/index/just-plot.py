import plotFuncs, pylab
plotFuncs.plotLogLogFromFile("crux-no-index-human.fasta", style="k-o",
    label="Crux w/o index (human)")
plotFuncs.plotLogLogFromFile("crux-no-index-human-dummy.fasta",
    style="k--o", label="Crux w/o index (yeast)")
plotFuncs.plotLogLogFromFile("sequest-human.fasta", style="m-s",
    label="Sequest (human)")
plotFuncs.plotLogLogFromFile("sequest-human-dummy.fasta", style="m--s",
    label="Sequest (yeast)")
plotFuncs.plotLogLogFromFile("crux-fast-human.fasta", style="g-^", 
    label="Crux (w/ index) human")
plotFuncs.plotLogLogFromFile("crux-fast-human-dummy.fasta", style="g--^", 
    label="Crux (w/ index) yeast")
pylab.legend(loc="lower right")
pylab.xlabel("Mass window (Da)")
pylab.xlim(0.05, 5)
pylab.ylim(0.1, 50000)
pylab.ylabel("Runtime for 1000 spectra (s)")
pylab.savefig("indexing.eps")
pylab.savefig("indexing.png")
