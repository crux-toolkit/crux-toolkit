Feature: generate-peptides
  generate-peptides should extract from a given set of protein sequences a list
    of target and decoy peptides fitting the specified criteria

Scenario Outline: User runs generate-peptides
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments <args>
  When I run generate-peptides
  Then the return value should be 0
  And crux-output/generate-peptides.target.txt should match good_results/<expected_output>

Examples:
  |test_name                                  |args                                                 |expected_output                                           |
  |generate_peptides_default                  |--overwrite T test.fasta                             |standard_generate_peptides_default.out                    |
  |generate_peptides_specifying_miss_cleavages|--parameter-file params/params_uniq test.fasta       |standard_generate_peptides_specifying_missed_cleavages.out|
  |generate_peptides_change_minmax            |--parameter-file params/many_changes test.fasta      |standard_generate_peptides_change_minmax                  |
  |generate_peptides_chymo                    |--parameter-file params/chymo test.fasta             |standard_generate_peptides_chymo                          |
  #|generate-peptides-ambiguous-residues       |--min-mass 1464 --max-mass 1496 ambiguous.fasta 2>&1|generate-peptides-ambiguous.target.txt                    |

