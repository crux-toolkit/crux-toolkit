Feature: search-for-xlinks
  search-for-xlinks should search a collection of spectra against a sequence
    database, returning a collection of matches corresponding to linear and
    cross-linked peptides scored by XCorr

Scenario Outline: User runs search-for-xlinks
  Given the path to Crux is ../../crux
  And I want to run a test named <test_name>
  And I pass the arguments <args> <spectra> <fasta> <sites> <mass>
  When I run search-for-xlinks
  Then the return value should be 0
  And crux-output/<actual_output> should match good_results/<expected_output>

Examples:
  |test_name        |args                            |spectra  |fasta         |sites  |mass  |actual_output     |expected_output        |
  |xlink-db         |--parameter-file params/xlink.db|xlink.ms2|xlink.db.fasta|K:K    |222   |xlink_peptides.txt|xlink_peptides.txt     |
  |search-for-xlinks|--parameter-file params/xlink   |xlink.ms2|xlink.fasta   |E:K,D:K|-18.01|search.target.txt |search-xlink.target.txt|

