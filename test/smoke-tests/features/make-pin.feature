Feature: make-pin
  make-pin should, given a set of search results files, generate a pin file for
    input to crux percolator

Scenario Outline: User runs make-pin
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T --output-file <output_file> <target_input>
  When I run make-pin
  Then the return value should be 0
  And crux-output/<output_file> should match good_results/<expected_output>

Examples:
  |test_name   |output_file     |target_input                 |expected_output |
  |make-pin-txt|make-pin_txt.pin|sample2.search.target.txt    |make-pin_txt.pin|
  |make-pin-pep|make-pin_pep.pin|sample2.search.target.pep.xml|make-pin_pep.pin|

