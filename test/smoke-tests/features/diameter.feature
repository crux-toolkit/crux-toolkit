Feature: diameter
  diameter should search a collection of DIA spectra against a sequence database or its index, returning a collection of peptide-spectrum matches (PSMs)

Scenario Outline: User runs diameter
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T --seed 7 <index_args> <fasta> <index>
  When I run tide-index as an intermediate step
  Then the return value should be 0
  And I pass the arguments --overwrite T <search_args> <spectra> <index>
  When I run diameter
  Then the return value should be 0
  And All lines in crux-output/<actual_output> should be in good_results/<expected_output> with 2 digits precision

Examples:
  |test_name       |index_args                    |search_args                              |spectra            |fasta            |index         |actual_output                     |expected_output             |
  |diameter_scale  |--decoy-format peptide-reverse|--top-match 5 --prec-ppm 10 --frag-ppm 10|diameter_test.mzXML|small-yeast.fasta|diameter_index|diameter.psm-features.txt         |diameter-search.scaled.txt  |
  |diameter_filter |--decoy-format peptide-reverse|--top-match 5 --prec-ppm 10 --frag-ppm 10|diameter_test.mzXML|small-yeast.fasta|diameter_index|diameter.psm-features.filtered.txt|diameter-search.filtered.txt|
