
schematic.eps: schematic.fig
	fig2dev -L eps schematic.fig schematic.eps

schematic.png: schematic.eps
	convert -density 144 -resize 40%x40% schematic.eps schematic.png

htmldir = $(prefix)/doc
cruximgdir = $(prefix)/doc

cruximg_DATA = \
  schematic.png

dist_html_DATA  =  \
  crux-analyze-matches.html \
  crux-create-index.html \
  crux-generate-peptides.html \
  crux-get-ms2-spectrum.html \
  crux.html \
  crux-predict-peptide-ions.html \
  crux-score-peptides-spectrum.html \
  crux-search-for-matches.html \
  customize-tutorial.html \
  home-tutorial.html \
  index-tutorial.html \
  install-tutorial.html \
  parameter-file.html \
  search-tutorial.html \
  start-tutorial.html \
  stylesheet.css \
  data/demo.ms2 \
  data/small-yeast.fasta

EXTRA_DIST = schematic.fig

WEBDIR=/noble/www/htdocs/proj/crux
webpage:
	cp crux.html $(WEBDIR)/index.html
	cp -f *.html $(WEBDIR)/
	cp *.css $(WEBDIR)/
	cp default.params $(WEBDIR)/
	rm -fr $(WEBDIR)/data/
	cp -R data/ $(WEBDIR)/
	cp -R download/ $(WEBDIR)/
	rm -fr $(WEBDIR)/*/CVS

