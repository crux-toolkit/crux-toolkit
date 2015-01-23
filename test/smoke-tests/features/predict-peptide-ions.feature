Feature: predict-peptide-ions
  predict-peptide-ions should, given a peptide and a charge state, predict the
    m/z values of the resulting fragment ions

Scenario Outline: User runs predict-peptide-ions
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments <args> <seq> <charge>
  When I run predict-peptide-ions
  Then the return value should be 0
  And stdout should match good_results/<expected_output>

Examples:
  |test_name                   |args                                             |seq    |charge|expected_output                          |
  |predict_ions_no_ops         |                                                 |IAMASEQ|2     |standard_predict_ions_no_ops.out         |
  |predict_ions_b_h2o          |--primary-ions b --precursor-ions T --h2o 1      |IAMASEQ|2     |standard_predict_ions_b_h2o.out          |
  |predict_ions_y_nh3          |--primary-ions y --nh3 1                         |IAMASEQ|2     |standard_predict_ions_y_nh3.out          |
  |predict_ions_by_flank_max_z1|--primary-ions by --max-ion-charge 1 --flanking T|IAMASEQ|3     |standard_predict_ions_by_flank_max_z1.out|

