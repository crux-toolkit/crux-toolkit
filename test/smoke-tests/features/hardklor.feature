Feature: hardklor
  hardklor should identify isotopic distributions from high-resolution mass
    spectra

Scenario Outline: User runs hardklor
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T <spectra>
  When I run hardklor
  Then the return value should be 0
  And crux-output/hardklor.mono.txt should match good_results/<expected_output>

Examples:
  |test_name       |spectra          |expected_output  |
  |hardklor-default|hardklor.test.ms1|hardklor.mono.txt|

