Feature: psm-convert
  psm-convert should read in a file containing peptide-spectrum matches in one
    of the variety of supported formats and output the same PSMs in a different
    format

Scenario Outline: User runs psm-convert
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T test_results/<results> <out_format>
  When I run psm-convert
  When I ignore lines matching the pattern: /^H[ \t]+StartTime.*$/
  When I ignore lines matching the pattern: /^<msms_pipeline_analysis .*$/
  When I ignore lines matching the pattern: /^<parameter .*$/
  When I ignore lines matching the pattern: /^<MzIdentML .*$/
  Then the return value should be 0
  And crux-output/<actual_output> should match good_results/<expected_output>

# Generate results with static-anywhere, variable-anywhere, static-nterm, static-cterm (results1.txt)
# tide-search --mods-spec C+57.02,2M+15.99 --nterm-peptide-mods-spec X+1.7 --cterm-peptide-mods-spec X+1.8 demo.ms2 small-yeast.fasta

# Generate results with static-anywhere, variable-nterm, variable-cterm (results2.txt)
# tide-search --mods-spec C+57.02 --nterm-peptide-mods-spec 1A+1.7 --cterm-peptide-mods-spec 1K+1.8 demo.ms2 small-yeast.fasta

Examples:
  |test_name           |results         |out_format|actual_output      |expected_output          |
  |psmconv-txt-to-txt1 |results1.txt    |tsv       |psm-convert.txt    |psmconv-from-txt1.txt    |
  |psmconv-txt-to-txt2 |results2.txt    |tsv       |psm-convert.txt    |psmconv-from-txt2.txt    |
  # From txt
  |psmconv-txt-to-html1|results1.txt    |html      |psm-convert.html   |psmconv-from-txt1.html   |
  |psmconv-txt-to-html2|results2.txt    |html      |psm-convert.html   |psmconv-from-txt2.html   |
  #|psmconv-txt-to-sqt1 |results1.txt    |sqt       |psm-convert.sqt    |psmconv-from-txt1.sqt    | TODO: sqt term mods
  #|psmconv-txt-to-sqt2 |results2.txt    |sqt       |psm-convert.sqt    |psmconv-from-txt2.sqt    | TODO: sqt term mods
  |psmconv-txt-to-pin1 |results1.txt    |pin       |psm-convert.pin    |psmconv-from-txt1.pin    |
  |psmconv-txt-to-pin2 |results2.txt    |pin       |psm-convert.pin    |psmconv-from-txt2.pin    |
  |psmconv-txt-to-pep1 |results1.txt    |pepxml    |psm-convert.pep.xml|psmconv-from-txt1.pep.xml|
  |psmconv-txt-to-pep2 |results2.txt    |pepxml    |psm-convert.pep.xml|psmconv-from-txt2.pep.xml|
  |psmconv-txt-to-mzid1|results1.txt    |mzidentml |psm-convert.mzid   |psmconv-from-txt1.mzid   |
  |psmconv-txt-to-mzid2|results2.txt    |mzidentml |psm-convert.mzid   |psmconv-from-txt2.mzid   |
  # To txt
  #|psmconv-sqt-to-txt1 |results1.sqt    |tsv       |psm-convert.txt    |psmconv-from-sqt1.txt    | TODO: sqt term mods
  #|psmconv-sqt-to-txt2 |results2.sqt    |tsv       |psm-convert.txt    |psmconv-from-sqt2.txt    | TODO: sqt term mods
  |psmconv-pep-to-txt1 |results1.pep.xml|tsv       |psm-convert.txt    |psmconv-from-pep1.txt    |
  |psmconv-pep-to-txt2 |results2.pep.xml|tsv       |psm-convert.txt    |psmconv-from-pep2.txt    |
  |psmconv-mzid-to-txt1|results1.mzid   |tsv       |psm-convert.txt    |psmconv-from-mzid1.txt   |
  |psmconv-mzid-to-txt2|results2.mzid   |tsv       |psm-convert.txt    |psmconv-from-mzid2.txt   |

