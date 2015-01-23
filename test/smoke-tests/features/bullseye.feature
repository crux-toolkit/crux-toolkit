Feature: bullseye
  bullseye should assign high resolution precursor m/z values to MS/MS data
    using the Hardklor algorithm

Scenario Outline: User runs bullseye
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T --parameter-file <param_file> <ms1_spectra> <ms2_spectra>
  When I run bullseye
  Then the return value should be 0
  And crux-output/bullseye.pid.ms2 should match good_results/<expected_output>

Examples:
  |test_name       |param_file             |ms1_spectra      |ms2_spectra      |expected_output |
  |bullseye-default|params/default-bullseye|hardklor.test.ms1|bullseye.test.ms2|bullseye.pid.ms2|

