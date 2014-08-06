Feature: generate-peptides
  generate-peptides should extract from a given set of protein sequences a list
    of target and decoy peptides fitting the specified criteria

Scenario Outline: User runs generate-peptides
  Given the path to Crux is ../../crux
  And I want to run a test named <test_name>
  And I pass the arguments <args>
  When I run generate-peptides
  Then the return value should be 0
  And stdout should match good_results/<expected_output>

Examples:
  |test_name                                  |args                                                                                                                               |expected_output                                                        |
  |generate_peptides_no_ops                   |test.fasta                                                                                                                         |standard_generate_peptides_no_ops.out                     |
  |generate_peptides_default                  |--output-sequence T test.fasta                                                                                                     |standard_generate_peptides_default.out                    |
  |generate_peptides_specifying_miss_cleavages|--parameter-file params/params_uniq --output-sequence T --enzyme trypsin --digestion partial-digest --missed-cleavages 3 test.fasta|standard_generate_peptides_specifying_missed_cleavages.out|
  |generate_peptides_change_minmax            |--parameter-file params/many_changes --output-sequence T --max-length 20 test.fasta                                                |standard_generate_peptides_change_minmax                  |
  |generate_peptides_mods                     |--parameter-file params/mods --output-sequence T --max-length 20 test.fasta                                                        |standard_generate_peptides_mods                           |
  |generate_peptides_fixed_mods               |--parameter-file params/fixed-mods --output-sequence T --max-length 20 test.fasta                                                  |standard_generate_peptides_fixed_mods                     |
  |generate_peptides_chymo                    |--parameter-file params/chymo --output-sequence T test.fasta                                                                       |standard_generate_peptides_chymo                          |
  #|generate-peptides-ambiguous-residues       |--min-mass 1464 --max-mass 1496 --output-sequence T ambiguous.fasta 2>&1                                                           |generate-peptides-ambiguous.target.txt                    |
  |generate_peptides_use_existing_index       |--output-sequence T --enzyme trypsin existing_crux_index                                                                           |standard_generate_peptides_use_exisiting_index.out        |

