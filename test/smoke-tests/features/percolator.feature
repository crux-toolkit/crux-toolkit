Feature: percolator
  percolator should re-rank a collection of PSMs using the percolator algorithm

Scenario Outline: User runs percolator
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T <args> <pin_file>
  When I run percolator
  Then the return value should be 0
  And All lines in crux-output/<actual_output> should be in good_results/<expected_output> with 5 digits precision

Examples:
  |test_name        |args                            |pin_file                     |actual_output                 |expected_output                       |
  |percolator-simple|--train-fdr 0.05 --test-fdr 0.05|sample2.search.target.txt.pin|percolator.target.peptides.txt|percolator.txt.pin.target.peptides.txt|
  |percolator-search-output|--search-input separate --train-fdr 0.05 --test-fdr 0.05|sample2.search.target.txt.pin|percolator.target.peptides.txt|percolator.txt.pin.target.peptides.txt|
  |percolator-static|--test-fdr 0.05 --init-weights percolator.weights.txt --static T|sample2.search.target.txt.pin|percolator.target.peptides.txt|static.percolator.target.peptides.txt|
