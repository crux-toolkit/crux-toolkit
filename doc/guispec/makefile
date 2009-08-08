all: guispec.pdf

%.eps: %.fig
	fig2dev -L eps $*.fig $*.eps

%.dvi: %.tex
	latex $*
	if ( grep -q citation $*.aux ) ; then \
		bibtex $* ; \
		latex $*.tex ; \
	fi
	latex $*

%.ps: %.dvi
	dvips -t letter $* -o $*.ps

%.pdf: %.ps
	ps2pdf -sPAPERSIZE=letter $*.ps


clean:
	rm -f *.dvi *.pdf *.aux *.blg *.log *.ps

