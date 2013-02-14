all: realtime.pdf

FILES = realtime.tex tide.tex tide-supplement.tex

realtime.pdf: $(FILES)

%.pdf: %.tex
	pdflatex $*.tex
	if ( grep -q citation $*.aux ) ; then \
		bibtex $* ; \
		pdflatex $*.tex ; \
	fi
	pdflatex $*.tex

%.pdf: %.fig
	fig2dev -L pdftex -p dummyarg $*.fig $*.pdf

clean:
	rm -f *.dvi *.aux *.blg *.log *.ps *.bbl *.pyc *~