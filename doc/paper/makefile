all: crux.pdf \
	plos-cover.pdf \
	jpr-cover.pdf \
	response.pdf \
	supp.pdf \
	crux.rtf

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

corrected-xcorr.dvi: corrected-xcorr.tex
	latex corrected-xcorr

plos-cover.dvi: plos-cover.tex
	latex plos-cover

jpr-cover.dvi: jpr-cover.tex
	latex jpr-cover

response.dvi: response.tex
	latex response

supp.dvi: supp.tex
	latex supp

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

%.rtf: %.tex
	latex2rtf $*
	latex2rtf $*

# Clean up.
clean:
	rm -f *.aux *.log *.blg *.dvi *.bbl *.ps *.eps
