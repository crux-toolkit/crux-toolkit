Feature: search-for-xlinks
  search-for-xlinks should search a collection of spectra against a sequence
    database, returning a collection of matches corresponding to linear and
    cross-linked peptides scored by XCorr

Scenario Outline: User runs search-for-xlinks
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments <args> <spectra> <fasta> <sites> <mass>
  When I run search-for-xlinks
  Then the return value should be 0
  And crux-output/<actual_output> should match good_results/<expected_output>

Examples:
  |test_name        |args                            |spectra  |fasta         |sites  |mass  |actual_output     |expected_output        |
  |xlink-db         |--parameter-file params/xlink.db|xlink.ms2|xlink.db.fasta|K:K    |222   |xlink_peptides.txt|xlink_peptides.txt     |
  |search-for-xlinks|--parameter-file params/xlink   |xlink.ms2|xlink.fasta   |E,D:K|-18.01|search-for-xlinks.target.txt|search-xlink.target.txt|
  #|search-for-xlinks-new|--parameter-file params/xlink-new|xlink.ms2|xlink.fasta   |E,D:K|-18.01|search-for-xlinks.txt|search-for-xlinks.new.txt|
  |search-for-xlinks-cz-ions|--parameter-file params/xlink-cz|xlink.ms2|xlink.fasta   |E,D:K|-18.01|search-for-xlinks.target.txt|search-for-xlinks.cz.txt|
  #|search-for-xlinks-ribo|--parameter-file params/xlink-ribo|good3.mgf|good3.fasta|K,nterm:K,nterm|136.100049|search-for-xlinks.txt|search-for-xlinks.ribo.txt|

# The search-for-xlinks-ribo test consists of three cross-linked spectra 
# with validated peptides from a ribosomal data set, provided by Jeff Howbert.
# For details, see the 29 June 2016 and 7 July 2016 entries here:
# http://noble.gs.washington.edu/~wnoble/proj/crux-projects/2016gaussian/results/bill/results.html
