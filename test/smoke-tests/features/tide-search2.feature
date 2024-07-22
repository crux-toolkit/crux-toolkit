# Please add any tests that have been added to this file to tide-search.feature
# Note that if tide-search2.feature is run before tide-search2.feature, then
# cucumber will fail due to Tide not dealing with ties well.

Feature: tide-search
  tide-index will not be called and instead tide-search will be directly given a fasta file
    subsequent calls to tide-search
  tide-search should search a collection of spectra against a sequence database,
    returning a collection of peptide-spectrum matches (PSMs)

Scenario Outline: User runs tide-search
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T --seed 7 <index_args> <search_args> <spectra> <fasta>
  When I run tide-search
  Then the return value should be 0
  And All lines in crux-output/<actual_output> should be in good_results/<expected_output> with 2 digits precision

Examples:
  |test_name      |index_args                    |search_args                                             |fasta            |spectra |actual_output         |expected_output    |
  # Tests that vary tide-search options
  |tide-ppmwin-fa    |                           |--precursor-window 5 --precursor-window-type ppm  --mz-bin-width 1.0005079     |small-yeast.fasta|demo.ms2|tide-search.target.txt|tide-ppmwin.txt    |
