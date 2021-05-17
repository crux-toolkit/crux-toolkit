Feature: tide-search
  tide-index will not be called and instead tide-search will be directly given a fasta file
    subsequent calls to tide-search
  tide-search should search a collection of spectra against a sequence database,
    returning a collection of peptide-spectrum matches (PSMs)

Scenario Outline: User runs tide-search
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T --seed 7 --file-column F <index_args> <search_args> <spectra> <fasta>
  When I run tide-search
  Then the return value should be 0
  And All lines in crux-output/<actual_output> should be in good_results/<expected_output> with 5 digits precision

Examples:
  |test_name      |index_args                                                   |search_args                                             |fasta            |spectra |actual_output         |expected_output    |
  # Tests that vary tide-search options
  |tide-masswin-fa   |                                                             |--precursor-window 5 --precursor-window-type mass --mz-bin-width 1.0005079     |tide_test_index|demo.ms2|tide-search.target.txt|tide-masswin.txt   |
  |tide-mzwin-fa     |                                                             |--precursor-window 5 --precursor-window-type mz   --mz-bin-width 1.0005079     |tide_test_index|demo.ms2|tide-search.target.txt|tide-mzwin.txt     |
  |tide-ppmwin-fa    |                                                             |--precursor-window 5 --precursor-window-type ppm  --mz-bin-width 1.0005079     |tide_test_index|demo.ms2|tide-search.target.txt|tide-ppmwin.txt    |
  |tide-computesp-fa |                                                             |--precursor-window 3 --precursor-window-type mass --compute-sp T --mz-bin-width 1.0005079                                         |tide_test_index|demo.ms2|tide-search.target.txt|tide-computesp.txt |
  |tide-specmz-fa    |                                                             |--precursor-window 3 --precursor-window-type mass --spectrum-min-mz 800 --spectrum-max-mz 900 --mz-bin-width 1.0005079            |tide_test_index|demo.ms2|tide-search.target.txt|tide-specmz.txt    |
  |tide-minpeaks-fa  |                                                             |--precursor-window 3 --precursor-window-type mass --min-peaks 100 --mz-bin-width 1.0005079                                        |tide_test_index|demo.ms2|tide-search.target.txt|tide-minpeaks.txt  |
  |tide-speccharge-fa|                                                             |--precursor-window 3 --precursor-window-type mass --spectrum-charge 3 --mz-bin-width 1.0005079                                    |tide_test_index|demo.ms2|tide-search.target.txt|tide-speccharge.txt|
  |tide-scannums-fa  |                                                             |--precursor-window 3 --precursor-window-type mass --scan-number 30-36 --mz-bin-width 1.0005079                                    |tide_test_index|demo.ms2|tide-search.target.txt|tide-scannums.txt  |
  |tide-rempeaks-fa  |                                                             |--precursor-window 3 --precursor-window-type mass --remove-precursor-peak T --remove-precursor-tolerance 3 --mz-bin-width 1.0005079|tide_test_index|demo.ms2|tide-search.target.txt|tide-rempeaks.txt  |
  |tide-useflank-fa  |                                                             |--precursor-window 3 --precursor-window-type mass --use-flanking-peaks T --mz-bin-width 1.0005079                                 |tide_test_index|demo.ms2|tide-search.target.txt|tide-useflank.txt  |
  |tide-usenl-fa     |                                                             |--precursor-window 3 --precursor-window-type mass --use-neutral-loss-peaks T --mz-bin-width 1.0005079                             |tide_test_index|demo.ms2|tide-search.target.txt|tide-usenl.txt     |
  |tide-mzbins-fa    |                                                             |--precursor-window 3 --precursor-window-type mass --mz-bin-width 0.02 --mz-bin-offset 0.34                                        |tide_test_index|demo.ms2|tide-search.target.txt|tide-mzbins.txt    |
  |tide-1thread-fa   |                                                             |--precursor-window 3 --precursor-window-type mass --num-threads 1 --mz-bin-width 1.0005079                                        |tide_test_index|demo.ms2|tide-search.target.txt|tide-default-1.txt   |
  |tide-7thread-fa   |                                                             |--precursor-window 3 --precursor-window-type mass --num-threads 7 --mz-bin-width 1.0005079                                        |tide_test_index|demo.ms2|tide-search.target.txt|tide-default-7.txt   |
  |tide-exact-pval-1thread-fa|                                                     |--precursor-window 3 --precursor-window-type mass --exact-p-value T --num-threads 1 --mz-bin-width 1.0005079                      |tide_test_index|demo.ms2|tide-search.target.txt|tide-exact-pval-1.txt|
  |tide-exact-pval-7thread-fa|                                                     |--precursor-window 3 --precursor-window-type mass --exact-p-value T --num-threads 7 --mz-bin-width 1.0005079                      |tide_test_index|demo.ms2|tide-search.target.txt|tide-exact-pval-7.txt|
  |tide-concat-fa    |                                                             |--precursor-window 3 --precursor-window-type mass --concat T --mz-bin-width 1.0005079                                             |tide_test_index|demo.ms2|tide-search.txt       |tide-concat.txt    |
  |tide-isoerr-fa    |                                                             |--precursor-window 3 --precursor-window-type mass --isotope-error 1,2,3 --mz-bin-width 1.0005079                                  |tide_test_index|demo.ms2|tide-search.target.txt|tide-isoerr.txt    |
  |tide-isoerrpval-fa|                                                             |--precursor-window 3 --precursor-window-type mass --isotope-error 1,2,3 --exact-p-value T --num-threads 1 --mz-bin-width 1.0005079 |tide_test_index|demo.ms2|tide-search.target.txt|tide-isoerrpval.txt|
  |tide-combine-pval-1thread-fa|                                                   |--precursor-window 3 --precursor-window-type mass --score-function both --exact-p-value T --num-threads 1 --mz-bin-width 1.0005079 |tide_test_index|demo.ms2|tide-search.target.txt|tide-combine-1.txt   |
  |tide-combine-pval-7thread-fa|                                                   |--precursor-window 3 --precursor-window-type mass --score-function both --exact-p-value T --num-threads 7 --mz-bin-width 1.0005079 |tide_test_index|demo.ms2|tide-search.target.txt|tide-combine-7.txt   |
  |tide-resEv-pval-1thread-fa|                                                     |--precursor-window 3 --precursor-window-type mass --score-function residue-evidence --exact-p-value T --num-threads 1 --use-neutral-loss-peaks F --mz-bin-width 1.0005079|tide_test_index|demo.ms2|tide-search.target.txt|tide-resEv-1.txt   |
  |tide-resEv-pval-7thread-fa|                                                     |--precursor-window 3 --precursor-window-type mass --score-function residue-evidence --exact-p-value T --num-threads 7 --use-neutral-loss-peaks F --mz-bin-width 1.0005079|tide_test_index|demo.ms2|tide-search.target.txt|tide-resEv-7.txt   |
  |tide-deiso-fa     |                                                             |--precursor-window 3 --precursor-window-type mass --deisotope 10 --mz-bin-width 1.0005079                                         |tide_test_index|demo.ms2|tide-search.target.txt|tide-deiso.txt     |
  |tide-deiso-pval-fa|                                                             |--precursor-window 3 --precursor-window-type mass --deisotope 10 --exact-p-value t --num-threads 1 --mz-bin-width 1.0005079       |tide_test_index|demo.ms2|tide-search.target.txt|tide-deiso-pval.txt|
  |tide-tailor-fa|                                                                 |--precursor-window 3 --precursor-window-type mass --num-threads 1 --mz-bin-width 1.0005079 --use-tailor-calibration T             |tide_test_index|demo.ms2|tide-search.target.txt|tide-tailor.txt|
  |tide-brief-fa|                                                                  |--precursor-window 3 --precursor-window-type mass --num-threads 1 --mz-bin-width 1.0005079 --brief-output T                       |tide_test_index|demo.ms2|tide-search.target.txt|tide-brief-output.txt|
  |tide-brief-centric-fa|                                                          |--precursor-window 3 --precursor-window-type mass --num-threads 1 --mz-bin-width 1.0005079 --brief-output T --peptide-centric-search T |tide_test_index|demo.ms2|tide-search.target.txt|tide-brief-peptide-centric.txt|
