
#schematic.eps: schematic.fig
#	fig2dev -L eps schematic.fig schematic.eps

schematic.png: schematic.eps
	convert -density 144 -resize 40%x40% schematic.eps schematic.png

htmldir = $(prefix)/doc
cruximgdir = $(prefix)/doc

if SEARCH
  SEARCH_ENABLED = true
else
  SEARCH_ENABLED =
endif

cruximg_DATA = \
  schematic.png

dist_html_DATA  =  \
	broken-link.html \
	compute-q-values.html \
	crux-analyze-matches.html \
	crux-create-index.html \
	crux-generate-peptides.html \
	crux-get-ms2-spectrum.html \
	crux.css \
	crux.html \
	crux-predict-peptide-ions.html \
	crux-search-for-matches.html \
	default.params \
	customize-tutorial.html \
	fasta-format.html \
	index-tutorial.html \
	install-tutorial.html \
	ms2-format.html \
	nibib.gif \
	parameter-file.html \
	percolator.html \
	q-ranker.html \
	release-notes.html \
	sample.fasta \
	schematic.png \
	search-tutorial.html \
	sqt-format.html \
	start-tutorial.html \
	txt-format.html \
	UWlogo150p.gif \
	data/demo.ms2 \
	data/small-yeast.fasta \
	data/search.target.txt \
	data/target.sqt

#EXTRA_DIST = schematic.fig
EXTRA_DIST = schematic.eps schematic.png

WEBDIR=/noble/www/htdocs/proj/crux
webpage:
	cp crux.html $(WEBDIR)/index.html
	cp -f *.html $(WEBDIR)/
	cp -f *.png $(WEBDIR)/
	cp -f *.gif $(WEBDIR)/
	cp *.css $(WEBDIR)/
	cp default.params $(WEBDIR)/
	rm -fr $(WEBDIR)/data/
	cp -R data/ $(WEBDIR)/
	cp -R open_source_download/ $(WEBDIR)/
	cp -R download/ $(WEBDIR)/
	rm -fr $(WEBDIR)/*/.svn

# If SEARCH_ENABLED is not defined this is an open-source installation
# and we replace the usual crux-search-for-matches.html file with
# the version for the open source distribution.
# Only needed when installing from SVN with search disabled.
install-data-hook:
	if test -z "$(SEARCH_ENABLED)"; \
    then  sed -e '/<!-- open source warning/d' \
              -e '/open source warning -->/d' \
              < crux-search-for-matches.html \
              > $(htmldir)/crux-search-for-matches.html; fi

# Similarly, subsitute files for open-source distribution
dist-hook:
	if test -z "$(SEARCH_ENABLED)"; \
    then  sed -e '/<!-- open source warning/d' \
              -e '/open source warning -->/d' \
              < crux-search-for-matches.html \
              > $(distdir)/crux-search-for-matches.html; fi
