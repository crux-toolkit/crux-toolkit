all: tide.pdf \
	tide-supplement.pdf \
	search-engines.pdf \
	coverletter.pdf \
	response.pdf \
	response2.pdf \
	response3.pdf \
	tidediff.pdf

FILES = tide.tex \
	tidediff.tex \
	tide-header.tex \
	table-incremental-optimizations.tex \
	timing_chart_main.pdf \
	timing_chart_mods.pdf \
	timing_chart_tocs.pdf \
	profiles.png \
	Diagrams_p1-p1.pdf \
	Diagrams_p2-p2.pdf \
	scatter/xcorr_xcorr_worm_seqorig_seqrecent.pdf \
	scatter/xcorr_xcorr_worm_seqorig_tide.pdf \
	scatter/xcorr_xcorr_worm_seqrecent_tide.pdf \
	scatter/xcorr_xcorr_yeast_seqorig_seqrecent.pdf \
	scatter/xcorr_xcorr_yeast_seqorig_tide.pdf \
	scatter/xcorr_xcorr_yeast_seqrecent_tide.pdf

tide.pdf: $(FILES)

tide-supplement.pdf: tide.pdf \
	tide-supplement.tex

tide.zip: $(FILES)
	zip tide.zip $(FILES)

table-incremental-optimizations.tex: table1.py
	python ./table1.py

profiles.png: profiles.py
	python ./profiles.py

Diagrams_p1-p1.pdf: Diagrams.pdf pdfextr
	. ./pdfextr 1 1 Diagrams.pdf

Diagrams_p2-p2.pdf: Diagrams.pdf pdfextr
	. ./pdfextr 2 2 Diagrams.pdf

%.pdf: %.gnuplot %.txt
	gnuplot $*.gnuplot > $*.pdf

%.gnuplot: %.txt
	cd scatter && ./make_gnuplots.py

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
	rm -f *.dvi *.aux *.blg *.log *.ps *.bbl scatter/*.gnuplot scatter/*.pdf *.pyc *~