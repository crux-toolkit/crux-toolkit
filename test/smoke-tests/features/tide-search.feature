Feature: tide-index / tide-search
  tide-index should create an index for all peptides in a fasta file, for use in
    subsequent calls to tide-search
  tide-search should search a collection of spectra against a sequence database,
    returning a collection of peptide-spectrum matches (PSMs)

Scenario Outline: User runs tide-index / tide-search
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T --seed 7 <index_args> <fasta> <index>
  When I run tide-index as an intermediate step
  Then the return value should be 0
  And I pass the arguments --overwrite T --file-column F <search_args> <spectra> <index>
  When I run tide-search
  Then the return value should be 0
  And crux-output/<actual_output> should contain the same lines as good_results/<expected_output>

Examples:
  |test_name      |index_args                                                   |search_args                                             |fasta            |index          |spectra |actual_output         |expected_output    |
  |tide-default   |                                                             |                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-default.txt   |
  # Tests that vary tide-index options
  |tide-peplen    |--min-length 5 --max-length 10                               |                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-peplen.txt    |
  |tide-pepmass   |--min-mass 1000 --max-mass 2000                              |                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-pepmass.txt   |
  |tide-avgmass   |--isotopic-mass average                                      |                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-avgmass.txt   |
  |tide-clipn     |--clip-nterm-methionine T                                    |                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-clipn.txt     |
  |tide-mods1     |--mods-spec C+57.02146,2M+15.9949,1STY+79.966331             |                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-mods1.txt     |
  |tide-mods1limit|--mods-spec C+57.02146,2M+15.9949,1STY+79.966331 --max-mods 1|                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-mods1limit.txt|
  |tide-mods1min  |--mods-spec C+57.02146,2M+15.9949,1STY+79.966331 --min-mods 1|                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-mods1min.txt  |
  |tide-modsn     |--nterm-peptide-mods-spec 1E-18.0106,C-17.0265               |                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-modsn.txt     |
  |tide-modsc     |--cterm-peptide-mods-spec X+21.9819                          |                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-modsc.txt     |
  |tide-chymo     |--enzyme chymotrypsin                                        |                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-chymo.txt     |
  |tide-partial   |--digestion partial-digest                                   |                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-partial.txt   |
  |tide-misscleave|--missed-cleavages 2                                         |                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-misscleave.txt|
  |tide-reverse   |--decoy-format peptide-reverse                               |                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-reverse.txt   |
  |tide-multidecoy|--num-decoys-per-target 5                                    |                                                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.decoy.txt |tide-5decoys.txt   |

  # Tests that vary tide-search options
  |tide-masswin   |                                                             |--precursor-window 5 --precursor-window-type mass       |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-masswin.txt   |
  |tide-mzwin     |                                                             |--precursor-window 5 --precursor-window-type mz         |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-mzwin.txt     |
  |tide-ppmwin    |                                                             |--precursor-window 5 --precursor-window-type ppm        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-ppmwin.txt    |
  |tide-computesp |                                                             |--compute-sp T                                          |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-computesp.txt |
  |tide-specmz    |                                                             |--spectrum-min-mz 800 --spectrum-max-mz 900             |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-specmz.txt    |
  |tide-minpeaks  |                                                             |--min-peaks 100                                         |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-minpeaks.txt  |
  |tide-speccharge|                                                             |--spectrum-charge 3                                     |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-speccharge.txt|
  |tide-scannums  |                                                             |--scan-number 30-36                                     |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-scannums.txt  |
  |tide-rempeaks  |                                                             |--remove-precursor-peak T --remove-precursor-tolerance 3|small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-rempeaks.txt  |
  |tide-useflank  |                                                             |--use-flanking-peaks T                                  |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-useflank.txt  |
  |tide-usenl     |                                                             |--use-neutral-loss-peaks T                              |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-usenl.txt     |
  |tide-mzbins    |                                                             |--mz-bin-width 0.02 --mz-bin-offset 0.34                |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-mzbins.txt    |
  |tide-exact-pval|                                                             |--exact-p-value T                                       |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-exact-pval.txt|
  |tide-1thread   |                                                             |--num-threads 1                                         |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-default.txt   |
  |tide-7thread   |                                                             |--num-threads 7                                         |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-default.txt   |
  |tide-exact-pval-1thread|                                                     |--exact-p-value T --num-threads 1                       |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-exact-pval.txt|
  |tide-exact-pval-7thread|                                                     |--exact-p-value T --num-threads 7                       |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-exact-pval.txt|
  |tide-concat    |                                                             |--concat T                                              |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.txt       |tide-concat.txt    |
  |tide-isoerr    |                                                             |--isotope-error 1,2,3                                   |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-isoerr.txt    |
  |tide-isoerrpval|                                                             |--isotope-error 1,2,3 --exact-p-value T                 |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-isoerrpval.txt|
  |tide-combine-pval-1thread|                                                   |--score-function both --exact-p-value T --num-threads 1 |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-combine-1.txt   |
  |tide-combine-pval-7thread|                                                   |--score-function both --exact-p-value T --num-threads 7 |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-combine-7.txt   |
  |tide-resEv-pval-1thread|                                                     |--score-function residue-evidence --exact-p-value T --num-threads 1 --use-neutral-loss-peaks F |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-resEv.txt   |
  |tide-resEv-pval-7thread|                                                     |--score-function residue-evidence --exact-p-value T --num-threads 7 --use-neutral-loss-peaks F |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-resEv.txt   |
  |tide-deiso     |                                                             |--deisotope 10                                          |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-deiso.txt     |
  |tide-deiso-pval|                                                             |--deisotope 10 --exact-p-value t                        |small-yeast.fasta|tide_test_index|demo.ms2|tide-search.target.txt|tide-deiso-pval.txt|
