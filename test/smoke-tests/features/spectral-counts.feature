Feature: spectral-counts
  spectral-counts should quantify peptides or proteins using one of three
    spectral counting methods

Scenario Outline: User runs spectral-counts
  Given the path to Crux is ../../src/crux
  And I want to run a test named <test_name>
  And I pass the arguments --overwrite T <args> <results>
  When I run spectral-counts
  Then the return value should be 0
  And crux-output/spectral-counts.target.txt should match good_results/<expected_output>

Examples:
  |test_name                     |args                                                                   |results                    |expected_output                   |
  |spectral-counts-raw           |--parameter-file params/raw --protein-database test.fasta              |test.target.txt            |spectral-counts.raw.txt           |
  |spectral-counts-sin           |--parameter-file params/sin --protein-database test.fasta              |test.target.txt            |spectral-counts.sin.txt           |
  |spectral-counts-nsaf          |--parameter-file params/nsaf --protein-database test.fasta             |test.target.txt            |spectral-counts.nsaf.txt          |
  |spectral-counts-empai         |--parameter-file params/empai --protein-database test.fasta            |test.target.txt            |spectral-counts.empai.txt         |
  |spectral-counts-dnsaf         |--parameter-file params/dnsaf --protein-database test.fasta            |test.target.txt            |spectral-counts.dnsaf.txt         |
  |spectral-counts-pepxml        |--parameter-file params/nsaf --protein-database small-yeast.fasta      |sample.target.pep.xml      |spectral-counts.pepxml.txt        |
  |spectral-counts-custom        |--parameter-file params/spc-custom --protein-database small-yeast.fasta|sample.target.pep.xml      |spectral-counts.pepxml.custom.txt |
  |spectral-counts-simple        |--parameter-file params/spc-simple --protein-database small-yeast.fasta|sample.target.pep.xml      |spectral-counts.simple.txt        |
  |spectral-counts-none          |--parameter-file params/spc-none --protein-database small-yeast.fasta  |sample.target.pep.xml      |spectral-counts.none.txt          |
  |spectral-counts-peptideprophet|--parameter-file params/spc-pp --protein-database small-yeast.fasta    |sample.target.pep.xml      |spectral-counts.peptideprophet.txt|
  |spectral-counts-pepxml-nodb   |--parameter-file params/raw                                            |sample.target.pep.xml      |spectral-counts.pepxml.nodb.txt   |
  |spectral-counts-tsv-nodb      |--parameter-file params/raw                                            |sample.target.psms.txt     |spectral-counts.raw.nodb.txt      |
  #|spectral-counts-mzid          |--threshold-type none --measure RAW                                    |Sequest_example_ver1.1.mzid|spectral-counts.mzid.txt          |
  |spectral-counts-barista       |--parameter-file params/nsaf --protein-database small-yeast.fasta      |sample2.target.psms.txt    |spectral-counts.barista.txt       |

