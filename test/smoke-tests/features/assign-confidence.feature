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
  |test_name        |args                                                     |input                      |actual_output               |expected_output                          |
  |xcorr            |                                                         |assign-default.target.txt  |assign-confidence.target.txt|assign-confidence-default.target.txt     |
  |concat           |                                                         |assign-concat.txt          |assign-confidence.target.txt|assign-confidence-concat.target.txt      |
  |auto-score       |                                                         |assign-exactpval.target.txt|assign-confidence.target.txt|assign-confidence-exact-tdc.target.txt   |
  |exact_pval       |--score "exact p-value"                                  |assign-exactpval.target.txt|assign-confidence.target.txt|assign-confidence-exact-tdc.target.txt   |
  |exact_pval_mixmax|--score "exact p-value" --estimation-method mix-max      |assign-exactpval.target.txt|assign-confidence.target.txt|assign-confidence-mixmax.target.txt      |
  |sidak            |--score "exact p-value" --sidak T                        |assign-exactpval.target.txt|assign-confidence.target.txt|assign-confidence-sidak.target.txt       |
  |peptide-level    |--score "exact p-value" --estimation-method peptide-level|assign-exactpval.target.txt|assign-confidence.target.txt|assign-confidence-peptidelevel.target.txt|
  |atdc             |                                                         |tide-5d.target.txt         |assign-confidence.target.txt|assign-confidence-atdc.target.txt        |

