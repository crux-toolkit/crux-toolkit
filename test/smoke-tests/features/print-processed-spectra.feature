Feature: print-processed-spectra
  print-processed-spectra should process spectra as for scoring xcorr and print
    the results to a file

Scenario Outline: User runs print-processed-spectra
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments <args> <input_spectra> <output_spectra>
  When I run print-processed-spectra
  Then the return value should be 0
  And crux-output/<output_spectra> should match good_results/<expected_output>

Examples:
  |test_name               |args                                         |input_spectra|output_spectra    |expected_output   |
  |print_processed_spectrum|--overwrite T --remove-precursor-tolerance 15|test.ms2     |processed-test.ms2|processed-test.ms2|

