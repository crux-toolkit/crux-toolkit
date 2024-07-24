Feature: assign-confidence
  assign-confidence evaluates the statistical confidence for each PSMs obtained 
  from tide-search.

Scenario Outline: User runs tide-index / tide-search
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T --seed 7 <args> <input>
  When I run assign-confidence
  Then the return value should be 0
  And crux-output/<actual_output> should contain the same lines as good_results/<expected_output>

Examples:
  |test_name        |args                                                     |input                                  |actual_output               |expected_output                          |
  |xcorr            |  --score "xcorr score"                                  |assign-default-tide-search.target.txt  |assign-confidence.target.txt|assign-confidence-default.target.txt     |
  |concat           |  --score "xcorr score"                                  |assign-concat-tide-search.txt          |assign-confidence.target.txt|assign-confidence-concat.target.txt      |
  |auto-score       |                                                         |assign-default-tide-search.target.txt  |assign-confidence.target.txt|assign-confidence-tailor-tdc.target.txt  |
  |combined_pval    |--score "combined p-value"                               |assign-pval-tide-search.target.txt     |assign-confidence.target.txt|assign-confidence-combined-tdc.target.txt|
  |exact_pval       |--score "exact p-value"                                  |assign-pval-tide-search.target.txt     |assign-confidence.target.txt|assign-confidence-exact-tdc.target.txt   |
  |atdc             |                                                         |assign-5dec-tide-search.target.txt     |assign-confidence.target.txt|assign-confidence-atdc.target.txt        |

