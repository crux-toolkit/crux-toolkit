cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

set (
  crux_documentation
	broken-link.html
	compute-q-values.html 
	crux-analyze-matches.html 
	crux-create-index.html 
	crux-generate-peptides.html 
	crux-get-ms2-spectrum.html 
	crux.css 
	crux.html 
	crux-predict-peptide-ions.html 
	crux-search-for-matches.html 
	crux-sequest-search.html 
  crux-search-for-xlinks.html 
	default.params 
	customize-tutorial.html 
	fasta-format.html 
	index-tutorial.html 
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
	UWlogo150p.gif 
	data/demo.ms2 
	data/small-yeast.fasta 
	data/search.target.txt 
	data/target.sqt
)
install (FILES ${crux_documentation} DESTINATION doc)
