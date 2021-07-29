Feature: diameter
  diameter should search a collection of DIA spectra against a sequence database or its index, returning a collection of peptide-spectrum matches (PSMs)

Scenario Outline: User runs diameter
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T <args> <spectra> <fasta>
  When I run diameter
  Then the return value should be 0
  And All lines in crux-output/<actual_output> should be in good_results/<expected_output> with 5 digits precision

Examples:
  |test_name           |args                                                      |spectra            |fasta            |actual_output               |expected_output           |
  |diameter_test1      |--top-match 5 --prec-ppm 10 --frag-ppm 10 --psm-filter F  |diameter_test.mzXML|small-yeast.fasta|diameter-search.tmp.txt     |diameter-search.tmp.txt   |
  |diameter_test2      |--top-match 5 --prec-ppm 10 --frag-ppm 10 --psm-filter F  |diameter_test.mzXML|small-yeast.fasta|diameter-search.scaled.txt  |diameter-search.scaled.txt|
  |diameter_test3      |--top-match 5 --prec-ppm 10 --frag-ppm 10 --psm-filter F  |diameter_test.mzXML|small-yeast.fasta|diameter-search.filtered.txt|diameter-search.filtered.txt|
