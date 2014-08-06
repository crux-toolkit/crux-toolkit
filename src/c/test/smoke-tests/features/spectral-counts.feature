Feature: spectral-counts
  spectral-counts should quantify peptides or proteins using one of three
    spectral counting methods

Scenario Outline: User runs spectral-counts
  Given the path to Crux is ../../crux
  And I want to run a test named <test_name>
  And I pass the arguments --protein-database test.fasta --parameter-file <param_file> test.target.txt
  When I run spectral-counts
  Then the return value should be 0
  And crux-output/spectral-counts.target.txt should match good_results/<expected_output>

Examples:
  |test_name            |param_file  |expected_output          |
  |spectral-counts-raw  |params/raw  |spectral-counts.raw.txt  |
  |spectral-counts-sin  |params/sin  |spectral-counts.sin.txt  |
  |spectral-counts-nsaf |params/nsaf |spectral-counts.nsaf.txt |
  |spectral-counts-empai|params/empai|spectral-counts.empai.txt|
  |spectral-counts-dnsaf|params/dnsaf|spectral-counts.dnsaf.txt|

