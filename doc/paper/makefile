all: crux.pdf

# Check out the 'refs' CVS module and add this directory to your BIBTEX path.
REFS = $(HOME)/refs/refs.bib

crux.bbl: crux.tex \
	$(REFS)

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
