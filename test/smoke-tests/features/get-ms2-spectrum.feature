Feature: get-ms2-spectrum
  get-ms2-spectrum should extract one or more fragmentation spectra, specified
    by scan number, from an MS2 file

Scenario Outline: User runs get-ms2-spectrum
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments <args> <spectra>
  When I run get-ms2-spectrum
  Then the return value should be 0
  And <actual_output> should match good_results/<expected_output>

Examples:
  |test_name                      |args                     |spectra       |actual_output|expected_output                                          |
  |get_ms2_spectrum               |--scan-number 2          |test.ms2      |stdout       |standard_get_ms2_spectrum.out               |
  |get_ms2_spectrum_stats         |--stats T --scan-number 2|test.ms2      |stdout       |standard_get_ms2_spectrum_stats.out         |
  |get_ms2_spectrum_without_zlines|--scan-number 2-16       |test-no-z1.ms2|stdout       |standard_get_ms2_spectrum_without_zlines.out|

