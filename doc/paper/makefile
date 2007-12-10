all: crux.pdf \
	plos-cover.pdf \
	jpr-cover.pdf \
	response.pdf

# Check out the 'refs' CVS module and add this directory to your BIBTEX path.
REFDIR = $(subst ::,,$(BIBINPUTS))
REFS = $(REFDIR)/refs.bib

crux.bbl: crux.tex \
	plos.bst \
	$(REFS) \
	schematic.eps \
	./Images/indexing.eps \
	./Images/random-sp.eps \
	./Images/random-xcorr.eps

plos-cover.dvi: plos-cover.tex
	latex plos-cover

jpr-cover.dvi: jpr-cover.tex
	latex jpr-cover

response.dvi: response.tex
	latex response

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
