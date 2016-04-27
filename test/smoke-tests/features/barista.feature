Feature: barista
  barista should re-rank a collection of PSMs using the barista algorithm

Scenario Outline: User runs barista
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T <args> <spectra> <search_results>
  When I run barista
  Then the return value should be 0
  And crux-output/<actual_output> should match good_results/<expected_output>

Examples:
  |test_name	   |args 	|spectra |search_results|actual_output |expected_output|
  |barista-simple|--separate-searches barista-args/barista-tide-search.decoy.txt barista-args/barista-concat.fasta|barista-args/barista-demo.ms2|barista-args/barista-tide-search.target.txt|barista.target.psms.txt|barista-barista.target.psms.txt|

