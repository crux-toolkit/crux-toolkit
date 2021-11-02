Feature: tide-index
  tide-index should create an index for all peptides in a fasta file, for use in
    subsequent calls to tide-search

Scenario Outline: User runs tide-index / tide-search
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T --peptide-list T <index_args> <fasta> <index>
  When I run tide-index
  Then the return value should be 0
  And crux-output/<actual_tds> should contain the same lines as good_results/<expected_tds>

Examples:
  |test_name      |index_args                                                                  |fasta            |index          |actual_tds             |expected_tds                 |
  |tide-default   |                                                                            |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-default.peptides.txt    |
  |tide-dups      |--allow-dups T                                                              |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-dups.peptides.txt       |
  |tide-temp-dir  |--temp-dir .                                                                |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-default.peptides.txt    |
  |tide-no-enzyme |--enzyme no-enzyme                                                          |test.fasta       |tide_test_index|tide-index.peptides.txt|tide-no-enzyme.peptides.txt  |
  |tide-mods      |--mods-spec 2M+15.9949,2STY+79.9663 --max-mods 2                            |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-index-mods1.peptides.txt|
  |tide-mods-alt  |--mods-spec 2M+15.9949,2STY+79.9663 --max-mods 2 --modsoutputter-threshold 1|small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-index-mods2.peptides.txt|
  |tide-multidecoy|--num-decoys-per-target 5                                                   |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-multidecoy.peptides.txt |
