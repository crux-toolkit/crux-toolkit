Feature: q-ranker
  q-ranker should re-rank a collection of PSMs using the Q-ranker algorithm

Scenario Outline: User runs q-ranker
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments <args> <spectra> <search_results>
  When I run q-ranker
  Then the return value should be 0
  And crux-output/<actual_output> should match good_results/<expected_output>

Examples:
  |test_name       |args                                                                                   |spectra |search_results           |actual_output           |expected_output       |
  |qranker-separate|--parameter-file params/set_rand_seed_only --separate-searches sample3.search.decoy.txt|demo.ms2|sample3.search.target.txt|q-ranker.target.psms.txt|qranker-sep.target.txt|

