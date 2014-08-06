Feature: compute-q-values
  compute-q-values should assign two types of statistical confidence measures
    (q-values and posterior error probabilities) to each PSM in a given set

Scenario Outline: User runs compute-q-values
  Given the path to Crux is ../../crux
  And I want to run a test named <test_name>
  And I pass the arguments <args> <target_input>
  When I run compute-q-values
  Then the return value should be 0
  And <actual_output> should match good_results/<expected_output>

Examples:
  |test_name       |args                                   |target_input             |actual_output                     |expected_output|
  |compute-q-values|--parameter-file params/pval           |sample.search.target.txt |crux-output/qvalues.target.txt    |qvalues.txt    |
  |decoy_qval      |--parameter-file params/decoy-qval-pval|sample4.search.target.txt|decoy-qval-pval/qvalues.target.txt|decoy-qval.txt |

