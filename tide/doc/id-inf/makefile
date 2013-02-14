all: id-inf.pdf \
	jpr-editor.pdf \
	rcms-editor.pdf

FILES=id-inf.tex \
	plots/margins.pdf \
	plots/e10-vs-e12.pdf \
	plots/1_excl2.pdf \
	plots/2b_variable2.pdf \
	plots/2_variable2.pdf \
	plots/3b_weib_triage2.pdf \
	plots/3_weib_triage2.pdf \
	plots/4b_chargestates2.pdf \
	plots/4_chargestates2.pdf

id-inf.pdf: $(FILES)

id-inf.zip: $(FILES) \
	id-inf.bbl
	zip id-inf.zip $(FILES)


%.pdf: %.tex
	pdflatex $*.tex
	if ( grep -q citation $*.aux ) ; then \
		bibtex $* ; \
		pdflatex $*.tex ; \
	fi
	pdflatex $*.tex

%.pdf: %.fig
	fig2dev -L pdftex -p dummyarg $*.fig $*.pdf

%.pdf: %.gnuplot
	gnuplot $*.gnuplot > $*.pdf

clean:
	rm -f *.dvi *.aux *.blg *.log *.ps
