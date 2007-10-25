all: crux.pdf

# Check out the 'refs' CVS module and add this directory to your BIBTEX path.
REFS = $(HOME)/refs/refs.bib

crux.bbl: crux.tex \
	plos.bst \
	$(REFS) \
	schematic.eps \
	./Images/indexing.eps \
	./Images/random-sp.eps \
	./Images/random-xcorr.eps

%.eps: %.fig
	fig2dev -L eps schematic.fig > schematic.eps

%.bbl: %.tex \
	$(REFS)
	latex $*
	latex $*
	bibtex $*

%.dvi: %.bbl
	latex $*
	latex $*

%.pdf: %.dvi
	dvipdf $*.dvi

# Clean up.
clean:
	rm -f *.aux *.log *.blg *.dvi *.bbl *.ps *.eps
