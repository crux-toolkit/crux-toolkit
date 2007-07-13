all: hmmseg-appnote.pdf \
	cover.pdf

# Check out the 'refs' CVS module and add this directory to your BIBTEX path.
REFS = $(HOME)/refs/refs.bib

hmmseg-appnote.bbl: hmmseg-appnote.tex \
	$(REFS)

cover.dvi: cover.tex
	latex cover
	latex cover

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
