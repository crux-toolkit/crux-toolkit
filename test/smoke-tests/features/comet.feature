Feature: comet
  comet should search a collection of spectra against a sequence database
    returning a collection of PSMs

Scenario Outline: User runs comet
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T --parameter-file <param_file> test_spectra/<spectra> test_databases/<fasta>
  When I run comet
  Then the return value should be 0
  And crux-output/comet.target.txt should match good_results/<expected_output>

Examples:
  |test_name         |param_file               |spectra                            |fasta                 |expected_output       |
  |comet_autodecoy   |params/comet_autodecoy   |JE102306_102306_18Mix4_Tube1_01.ms2|comet_autodecoy.fasta |comet_autodecoy.txt   |
  |comet_autodecoy1  |params/comet_autodecoy1  |JE102306_102306_18Mix4_Tube1_01.ms2|comet_autodecoy.fasta |comet_autodecoy1.txt  |
  |comet_binarymod2  |params/comet_binarymod2  |sh_1617_JX_070209p_KO410_run2.ms2  |comet_binarymod2.fasta|comet_binarymod2.txt  |
  |comet_binarymods  |params/comet_binarymods  |sh_1617_JX_070209p_KO410_run2.ms2  |comet_binarymods.fasta|comet_binarymods.txt  |
  |comet_commandline |params/comet_commandline |sh_1617_JX_070209p_KO410_run2.ms2  |comet_plain.fasta     |comet_commandline.txt |
  |comet_ctermmod    |params/comet_ctermmod    |sh_1617_JX_070209p_KO410_run1.ms2  |comet_term-mod.fasta  |comet_ctermmod.txt    |
  |comet_noenzyme    |params/comet_noenzyme    |sh_1617_JX_070209p_KO410_run1.ms2  |comet_noenzyme.fasta  |comet_noenzyme.txt    |
  |comet_plain       |params/comet_plain       |sh_1617_JX_070209p_KO410_run2.ms2  |comet_plain.fasta     |comet_plain.txt       |
  |comet_semi-tryptic|params/comet_semi-tryptic|sh_1617_JX_070209p_KO410_run1.ms2  |comet_tryptic.fasta   |comet_semi-tryptic.txt|
  |comet_stop-codon  |params/comet_stop-codon  |sh_1617_JX_070209p_KO410_run2.ms2  |comet_stop-codon.fasta|comet_stop-codon.txt  |
  |comet_term-mod    |params/comet_term-mod    |sh_1617_JX_070209p_KO410_run1.ms2  |comet_term-mod.fasta  |comet_term-mod.txt    |
  |comet_term-mod2   |params/comet_noenzyme    |sh_1617_JX_070209p_KO410_run1.ms2  |comet_term-mod.fasta  |comet_term-mod2.txt   |
  |comet_tryptic     |params/comet_tryptic     |sh_1617_JX_070209p_KO410_run1.ms2  |comet_tryptic.fasta   |comet_tryptic.txt     |

