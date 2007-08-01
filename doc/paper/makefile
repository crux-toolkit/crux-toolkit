all: crux.pdf

# Check out the 'refs' CVS module and add this directory to your BIBTEX path.
REFS = $(HOME)/refs/refs.bib

crux.bbl: crux.tex \
	$(REFS) \
	indexing.eps

indexing.eps: ../../results/paper-figure/index/runall \
	../../results/paper-figure/index/test.ms2 \
	../../results/paper-figure/index/test.fasta \
	../../results/paper-figure/index/ms22dta.pl \
	../../results/paper-figure/index/compare_runtime.py
	cd ../../results/paper-figure/index; ./runall
	cp ../../results/paper-figure/index/indexing.eps .

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
