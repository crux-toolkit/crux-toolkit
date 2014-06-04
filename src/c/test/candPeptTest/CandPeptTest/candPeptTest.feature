Feature: candidate peptide selection test
  In order to compare the selected candidate peptides in Comet and Tide
  As a curious tester
  I want to see the difference

Scenario Outline: User runs Tide and Comet
  Given the path to Crux is crux
  And the fasta file is <fasta>
  And the spectra file is <spectra>
  And I have entered the mass_tol_size <mass_tol>
  And I have entered the mass_tol_type <mass_tol_type>
  And I have entered the missed_cleavages <missed_cleavages>
  When I run comet and tide
  Then Candidate peptides should be identical

  Examples:
  
  |fasta           |spectra | mass_tol | mass_tol_type | missed_cleavages |
  |test.fasta      |test.ms2|  30      |mass           | 0  			  |
  