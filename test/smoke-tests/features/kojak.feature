Feature: kojak
  kojak should search a collection of spectra against a sequence database returning a collection of cross-linked PSMs.

Scenario Outline: User runs kojak
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T --parameter-file params/kojak_params <args> test_spectra/kojak_arp.mzML.gz test_databases/arp.fasta
  When I run kojak
  Then the return value should be 0
  And crux-output/kojak_arp.kojak.txt should match good_results/<expected_output>

Examples:
  |test_name           |args                                                        |expected_output        |
  |kojak_autodecoy     |--decoy_filter "DECOY_ 1"                                   |kojak_autodecoy.txt    |
  |kojak_no-decoy      |--decoy_filter "DECOY_ 0"                                   |kojak_no-decoy.txt     |
  |kojak_multi-enzyme  |--kojak_enzyme "[KR] Tryp, [DE]\|{P} GluC"                  |kojak_multi-enzyme.txt |
  |kojak_multi-xlink   |--cross_link "nK nK 138.0681 BS3, nK DE -18.0106 EDC"       |kojak_multi-xlink.txt  |
  |kojak_param-medic   |--auto_ppm_tolerance_pre warn --auto_fragment_bin_size warn |kojak_param-medic.txt  |