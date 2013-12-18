cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

# Only create webpage directory if in svn workspace
if (EXISTS /noble/www/htdocs/proj/crux)
  set(WEBDIR /noble/www/htdocs/proj/crux)
  add_custom_target(
    webpage
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/user
    COMMAND cp crux.html ${WEBDIR}/index.html
    COMMAND cp -f *.html ${WEBDIR}/
    COMMAND cp -f *.png ${WEBDIR}/
    COMMAND cp -f *.gif ${WEBDIR}/
    COMMAND cp *.css ${WEBDIR}/
    COMMAND cp default.params ${WEBDIR}/
    COMMAND rm -fr ${WEBDIR}/data/
    COMMAND mkdir ${WEBDIR}/data
    COMMAND cp data/demo.ms2 ${WEBDIR}/data
    COMMAND cp data/search.target.txt ${WEBDIR}/data
    COMMAND cp data/small-yeast.fasta ${WEBDIR}/data
    COMMAND cp data/target.sqt ${WEBDIR}/data
    COMMAND mkdir -p ${WEBDIR}/download
    COMMAND cp download/index.html ${WEBDIR}/download
  )
endif (EXISTS /noble/www/htdocs/proj/crux)

set (
  crux_documentation
	UWlogo150p.gif 
	broken-link.html
	compute-q-values.html 
	crux-analyze-matches.html 
	crux-create-index.html 
	crux-get-ms2-spectrum.html 
	crux-predict-peptide-ions.html 
	crux-search-for-matches.html 
	crux-sequest-search.html 
	crux.css 
	crux.html 
	customize-tutorial.html 
	data/demo.ms2 
	data/search.target.txt 
	data/small-yeast.fasta 
	data/target.sqt
	default.params 
	fasta-format.html 
	generate-peptides.html 
	install-tutorial.html 
	ms2-format.html 
	nibib.gif 
	parameter-file.html 
	percolator.html 
	q-ranker.html 
	release-notes.html 
	sample.fasta 
	schematic.png 
	search-tutorial.html 
	sqt-format.html 
	start-tutorial.html 
	txt-format.html 
	xlink-txt-format.html 
        crux-search-for-xlinks.html 
        spectral-counts-tutorial.html
        spectral-counts.html
)
install (FILES ${crux_documentation} DESTINATION doc)
