Feature: pipeline
  pipeline runs bullseye --> tide-search/comet --> percolator --> spectral-counts
  making up a complete analysis pipeline.

Scenario Outline: User runs crux pipeline
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T --seed 1 <args> <input>
  When I run pipeline
  Then the return value should be 0
  And crux-output/<actual_output> should contain the same lines as good_results/<expected_output>

Examples:
  |test_name                    |args                                               |input                        |actual_output              |expected_output                                                |
  |pipeline-percoator           |                                                   |pipeline-demo.ms2 small-yeast.fasta   |spectral-counts.target.txt |pipeline-percolator-spectral-counts.target.txt        |
  |pipeline-percoator-fileroot  | --fileroot foo                                    |pipeline-demo.ms2 small-yeast.fasta   |spectral-counts.target.txt |pipeline-percolator-spectral-counts.target.txt        |
  |pipeline-comet               | --search-engine comet --decoy_search 1            |pipeline-demo.ms2 small-yeast.fasta   |spectral-counts.target.txt |pipeline-comet-percolator-spectral-counts.target.txt        | 

# LFQ runs as the last pipeline step and needs MS1 scans, so it is only
# exercised with mzML input (the ms2 runs above skip it with a warning).
# lfq-pipeline-human.fasta: the 416 human proteins identified at q <= 0.01 when
# these two runs are searched against UniProtKB/Swiss-Prot human (release
# 2026_03, organism_id:9606 AND reviewed:true), plus 416 other proteins from
# that release drawn with random.seed(20261030), sorted by accession.
# UniProt data are licensed CC BY 4.0.
Scenario: User runs crux pipeline with label-free quantification
  Given the path to Crux is ../../src/crux
  And I want to run a test named pipeline-lfq
  And I pass the arguments --overwrite T --seed 1 --output-dir pipeline-lfq 20100614_Velos1_TaGe_SA_K562_3.mzML 20100614_Velos1_TaGe_SA_K562_4.mzML lfq-pipeline-human.fasta
  When I run pipeline
  Then the return value should be 0
  And pipeline-lfq/crux-lfq-peaks.txt should contain the same lines as good_results/pipeline-lfq-peaks.txt
  And pipeline-lfq/crux-lfq-mod-pep.txt should contain the same lines as good_results/pipeline-lfq-mod-pep.txt
