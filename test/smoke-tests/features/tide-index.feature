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
  |test_name       |index_args                                                                 |fasta            |index          |actual_tds             |expected_tds                 |
 # |tide-multidecoy|--num-decoys-per-target 5                                |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-multidecoy.peptides.txt |
  |tide-default    |                                                         |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-default.peptides.txt    |
  |tide-dups       |--enzyme no-enzyme --allow-dups T  --max-mass 3500       |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-dups.peptides.txt       |
  |tide-temp-dir   |--temp-dir .                                             |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-default.peptides.txt    |
  |tide-minmax-len |--min-length 10 --max-length 15                          |test.fasta       |tide_test_index|tide-index.peptides.txt|tide-minmaxlen.peptides.txt  |
  |tide-minmax-mass|--min-mass 1000 --max-mass 1500                          |test.fasta       |tide_test_index|tide-index.peptides.txt|tide-minmaxmass.peptides.txt |
  |tide-no-enzyme  |--enzyme no-enzyme                                       |test.fasta       |tide_test_index|tide-index.peptides.txt|tide-no-enzyme.peptides.txt  |
  |tide-mods       |--mods-spec 2M+15.9949,2STY+79.9663 --max-mods 2         |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-index-mods1.peptides.txt|
  |tide-nodecoy    |--decoy-format none                                      |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-index-nodecoys.peptides.txt|
  |tide-peptrev    |--decoy-format peptide-reverse                           |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-index-peptide-reverse.peptides.txt|
  |tide-shuffle    |--decoy-format shuffle                                   |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-index-shuffle.peptides.txt|
  |tide-full       |--digestion full-digest                                  |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-index-full.peptides.txt |
  |tide-partial    |--digestion partial-digest                               |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-index-partial.peptides.txt|
  |tide-non-spec   |--digestion non-specific-digest  --max-mass 3500         |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-index-nonspec.peptides.txt|
  |tide-mc4        |--missed-cleavages 4                                     |small-yeast.fasta|tide_test_index|tide-index.peptides.txt|tide-index-mc4.peptides.txt|
