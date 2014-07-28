Feature: search-for-matches
  search-for-matches should search a collection of spectra against a sequence
    database, returning a collection of peptide-spectrum matches (PSMs) scored
    by XCorr

Scenario Outline: User runs search-for-matches
  Given the path to Crux is ../../crux
  And I want to run a test named <test_name>
  And I pass the arguments --parameter-file <param_file> <args> <spectra> <fasta>
  When I run search-for-matches
  Then the return value should be 0
  And <actual_output> should match good_results/<expected_output>

Examples:
  |test_name                                |param_file                |args                                                |spectra            |fasta                     |actual_output                           |expected_output                                     |
  |search_for_matches_default               |params/set_rand_seed      |--decoys peptide-shuffle                            |test.ms2           |test.fasta                |crux-output/search.target.txt           |search_for_matches_default.target.txt               |
  |search_duplicates                        |params/set_rand_seed      |--decoys peptide-shuffle --fileroot duplicates      |test.ms2           |duplicates.fasta          |crux-output/duplicates.search.target.txt|search_for_matches_duplicates.target.txt            |
  |search_sp                                |params/sp                 |--decoys peptide-shuffle --output-dir crux-output-sp|test.ms2           |test.fasta                |crux-output-sp/search.target.txt        |sp.search.target.txt                                |
  |search_pvalues                           |params/pval               |                                                    |demo.ms2           |small-yeast.fasta         |crux-output/search.target.txt           |search_pvalues.target.txt                           |
  |search_mods_fasta                        |params/mods               |--decoys peptide-shuffle                            |test.ms2           |test.fasta                |crux-output/search.target.txt           |search_mods_fasta.target.txt                        |
  |reverse_sequence                         |params/reverse            |                                                    |test.ms2           |test-plus-palindrome.fasta|crux-output/reverse.search.decoy.txt    |reverse.decoy.txt                                   |
  |search_elastase                          |params/elastase           |                                                    |test.ms2           |test.fasta                |crux-output/search.target.txt           |search_elastase.target.txt                          |
  |search_custom_like_elastase              |params/custom-enzyme      |                                                    |test.ms2           |test.fasta                |crux-output/search.target.txt           |search_custom_like_elastase.target.txt              |
  |search_tdc                               |params/tdc                |                                                    |test.ms2           |test.fasta                |tdc/search.target.txt                   |search_tdc.target.txt                               |
  |decoys_one_file                          |params/one-decoy-file     |                                                    |test.ms2           |test.fasta                |one-decoy-file/search.decoy.txt         |decoys_one_file.decoy.txt                           |
  |mz_window                                |params/mz_window          |--fileroot mzwin                                    |demo.ms2           |small-yeast.fasta         |crux-output/mzwin.search.target.txt     |mz-window.target.txt                                |
  |ppm_window                               |params/ppm_window         |--fileroot ppmwin                                   |demo.ms2           |small-yeast.fasta         |crux-output/ppmwin.search.target.txt    |ppmwin.search.target.txt                            |
  |search_for_matches_default_existing_index|params/set_rand_seed      |--decoys peptide-shuffle                            |test.ms2           |existing_crux_index       |crux-output/search.target.txt           |search_for_matches_default_existing_index.target.txt|
  |search_mods_index                        |params/mods-high-precision|--decoys peptide-shuffle                            |test.ms2           |existing_crux_index       |crux-output/search.target.txt           |search_mods_index.target.txt                        |
  |search_select_scans                      |params/minimal            |--scan-number 150-153                               |demo.ms2           |small-yeast.fasta         |crux-output/search.target.txt           |search_select_scans.target.txt                      |
  |mzxml                                    |params/mstk               |--fileroot mzxml                                    |small.raw2xml.mzXML|small-yeast.fasta         |crux-output/mzxml.search.target.txt     |mzxml.search.target.txt                             |
  |mgf                                      |params/minimal            |--fileroot mgf                                      |test.mgf           |test.fasta                |crux-output/mgf.search.target.txt       |mgf.search.target.txt                               |

