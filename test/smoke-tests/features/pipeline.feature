Feature: pipeline
  pipeline runs bullseye --> tide-search/comet --> percolator/assign-confidence --> spectral-counts
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
  |pipeline-assign-confidence   | --post-processor assign-confidence                |pipeline-demo.ms2 small-yeast.fasta   |spectral-counts.target.txt |pipeline-assign-confidence-spectral-counts.target.txt |
  |pipeline-assign-confidence   | --fileroot foo --post-processor assign-confidence |pipeline-demo.ms2 small-yeast.fasta   |spectral-counts.target.txt |pipeline-assign-confidence-spectral-counts.target.txt |

