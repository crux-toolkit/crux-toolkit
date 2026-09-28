Feature: tide-search pepXML output
  tide-search should write the native spectrum identifier of each spectrum
    (spectrumNativeID) and the protein database used (search_database) to
    its pepXML output (issue #535)

  # sliced-mzml.mzML:  Thermo native ids ("controllerType=0 controllerNumber=1 scan=N")
  # sliced-sciex.mzML: Sciex WIFF native ids ("sample=1 period=1 cycle=N experiment=M"),
  #   100 MS2 spectra converted with msconvert from 20230221_Nanxi_bsa.wiff,
  #   PRIDE PXD066231 (CC0)

Scenario Outline: User runs tide-search with pepXML output
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T --seed 7 --num-threads 1 --pepxml-output T <spectra> small-yeast.fasta
  When I run tide-search
  When I ignore lines matching the pattern: /^<msms_pipeline_analysis .*$/
  When I ignore lines matching the pattern: /^<parameter .*$/
  When I ignore lines matching the pattern: /^<search_database local_path=".+small-yeast\.fasta" type="AA" \/>$/
  Then the return value should be 0
  And crux-output/tide-search.target.pep.xml should match good_results/<expected_output>

Examples:
  |test_name          |spectra          |expected_output               |
  |tide-pepxml-thermo |sliced-mzml.mzML |tide-pepxml-thermo.pep.xml    |
  |tide-pepxml-sciex  |sliced-sciex.mzML|tide-pepxml-sciex.pep.xml     |
