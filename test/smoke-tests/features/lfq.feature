Feature: lfq
  crux lfq performs label-free quantification (LFQ) of peptides and proteins
  from MS1 spectra, following the FlashLFQ algorithm
  (Millikin et al., Journal of Proteome Research, 2018).

  Input: a PSM file (assign-confidence or percolator format) plus one or more
  mzML spectrum files.
  Output: crux-lfq-mod-pep.txt (per-peptide intensities) and
          crux-lfq-peaks.txt  (per-peak details).

  Note: Match-Between-Runs (MBR) is not yet implemented.

  Note: golden-file comparison steps are omitted here because the exact
  numerical output depends on peak-detection results and will be added in a
  follow-up commit once representative golden files have been generated.
  The scenarios below verify that the pipeline runs to completion without
  crashing for the two primary configurations (with and without normalization).

Scenario: User runs crux lfq without normalization
  Given the path to Crux is ../../src/crux
  And I want to run a test named lfq-no-norm
  And I pass the arguments --overwrite T --output-dir lfq-no-norm --psm-file-format assign-confidence --is-psm-filtered T lfq-psms-no-norm.tsv sliced-mzml.mzML
  When I run lfq
  Then the return value should be 0
  And lfq-no-norm/crux-lfq-peaks.txt should contain the same lines as good_results/lfq-no-norm-peaks.txt
  And lfq-no-norm/crux-lfq-mod-pep.txt should contain the same lines as good_results/lfq-no-norm-mod-pep.txt

Scenario: User runs crux lfq with normalization
  Given the path to Crux is ../../src/crux
  And I want to run a test named lfq-with-norm
  And I pass the arguments --overwrite T --output-dir lfq-with-norm --normalize T --psm-file-format assign-confidence --is-psm-filtered T --specfile-replicates lfq-specfile-replicates.tsv lfq-psms-with-norm.tsv 20100614_Velos1_TaGe_SA_K562_3.mzML 20100614_Velos1_TaGe_SA_K562_4.mzML
  When I run lfq
  Then the return value should be 0
  And lfq-with-norm/crux-lfq-peaks.txt should contain the same lines as good_results/lfq-with-norm-peaks.txt
  And lfq-with-norm/crux-lfq-mod-pep.txt should contain the same lines as good_results/lfq-with-norm-mod-pep.txt
