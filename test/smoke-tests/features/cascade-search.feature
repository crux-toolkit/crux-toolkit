Feature: tide-index / cascade-search
  tide-index should create an index for all peptides in a fasta file, for use in
    subsequent calls to tide-search
  cascade-search should search a collection of spectra against a sequence database,
    returning a collection of peptide-spectrum matches (PSMs)

Scenario Outline: User runs tide-index / cascade-search
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T small-yeast.fasta small_yeast_index
  When I run tide-index as an intermediate step
  Then the return value should be 0
  And I pass the arguments --overwrite T --missed-cleavages 1 small-yeast.fasta small_yeast_index_mc1
  When I run tide-index as an intermediate step
  Then the return value should be 0
  And I pass the arguments --overwrite T <search_args> <spectra> <index>
  When I run cascade-search
  Then the return value should be 0
  And crux-output/<actual_output> should contain the same lines as good_results/<expected_output>


Examples:
  |test_name      |search_args                                                     |index                                  |spectra |actual_output            |expected_output     |
  |cascade-default|                                                                |small_yeast_index                      |demo.ms2|cascade-search.target.txt|cascade-default.txt |
  |cascade-2-index|                                                                |small_yeast_index,small_yeast_index_mc1|demo.ms2|cascade-search.target.txt|cascade-2-index.txt |
  |cascade-q-value|--q-value-threshold 0.0005                                      |small_yeast_index,small_yeast_index_mc1|demo.ms2|cascade-search.target.txt|cascade-q-value.txt |
  |cascade-min-max|--estimation-method min-max                                     |small_yeast_index,small_yeast_index_mc1|demo.ms2|cascade-search.target.txt|cascade-min-max.txt |
  |cascade-pep-lvl|--estimation-method peptide-level                               |small_yeast_index,small_yeast_index_mc1|demo.ms2|cascade-search.target.txt|cascade-pep-lvl.txt |
  |cascade-extpval|--score "exact p-value" --exact-p-value T                       |small_yeast_index,small_yeast_index_mc1|demo.ms2|cascade-search.target.txt|cascade-extpval.txt |
  |cascade-sidak  |--score "exact p-value" --exact-p-value T --sidak T             |small_yeast_index,small_yeast_index_mc1|demo.ms2|cascade-search.target.txt|cascade-sidak.txt   |
  |cascade-comb-cs|--estimation-method peptide-level --combine-charge-states T     |small_yeast_index,small_yeast_index_mc1|demo.ms2|cascade-search.target.txt|cascade-comb-cs.txt |
  |cascade-comb-mp|--estimation-method peptide-level --combine-modified-peptides T |small_yeast_index,small_yeast_index_mc1|demo.ms2|cascade-search.target.txt|cascade-comb-mp.txt |
  |cascade-file-column|--file-column F                                             |small_yeast_index                      |demo.ms2|cascade-search.target.txt|cascade-file-column.txt|
