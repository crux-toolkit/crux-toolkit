Feature: tide-index
  tide-index should create an index for all peptides in a fasta file, for use in
    subsequent calls to tide-search

Scenario Outline: User runs tide-index / tide-search
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T --peptide-list T <index_args> <fasta> <index>
  When I run tide-index
  Then the return value should be 0
  And crux-output/<actual_targets> should contain the same lines as good_results/<expected_targets>
  And crux-output/<actual_decoys> should contain the same lines as good_results/<expected_decoys>

Examples:
  |test_name      |index_args                                                                  |fasta            |index          |actual_targets                |expected_targets           |actual_decoys                |expected_decoys           |
  |tide-default   |                                                                            |small-yeast.fasta|tide_test_index|tide-index.peptides.target.txt|tide-default.target.txt    |tide-index.peptides.decoy.txt|tide-default.decoy.txt    |
  |tide-dups      |--allow-dups T                                                              |small-yeast.fasta|tide_test_index|tide-index.peptides.target.txt|tide-dups.target.txt       |tide-index.peptides.decoy.txt|tide-dups.decoy.txt       |
  |tide-proteinReverse|--decoy-format PROTEIN-REVERSE                                          |small-yeast.fasta|tide_test_index|tide-index.peptides.target.txt|tide-reverse.target.txt    |tide-index.peptides.decoy.txt|tide-reverse.decoy.txt    |
  |tide-temp-dir  |--temp-dir .                                                                |small-yeast.fasta|tide_test_index|tide-index.peptides.target.txt|tide-default.target.txt    |tide-index.peptides.decoy.txt|tide-default.decoy.txt    |
  |tide-no-enzyme |--enzyme no-enzyme                                                          |test.fasta       |tide_test_index|tide-index.peptides.target.txt|tide-no-enzyme.target.txt  |tide-index.peptides.decoy.txt|tide-no-enzyme.decoy.txt  |
  |tide-mods      |--mods-spec 2M+15.9949,2STY+79.9663 --max-mods 2                            |small-yeast.fasta|tide_test_index|tide-index.peptides.target.txt|tide-index-mods1.target.txt|tide-index.peptides.decoy.txt|tide-index-mods1.decoy.txt|
  |tide-mods-alt  |--mods-spec 2M+15.9949,2STY+79.9663 --max-mods 2 --modsoutputter-threshold 1|small-yeast.fasta|tide_test_index|tide-index.peptides.target.txt|tide-index-mods1.target.txt|tide-index.peptides.decoy.txt|tide-index-mods1.decoy.txt|
  |tide-multidecoy|--num-decoys-per-target 5                                                   |small-yeast.fasta|tide_test_index|tide-index.peptides.target.txt|tide-default.target.txt    |tide-index.peptides.decoy.txt|tide-index-multi.decoy.txt|

