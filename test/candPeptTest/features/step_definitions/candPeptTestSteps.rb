# encoding: utf-8
begin require 'rspec/expectations'; rescue LoadError; require 'spec/expectations'; end 
#require 'cucumber/formatter/unicode'

Before do
  @Tester = CandPeptTest.new();
end

Given /the path to [Cc]rux is (.*)/ do |path|
	@Tester.crux_path(path);
end

Given /the fasta file is (.*)/ do |fasta|
	@Tester.fasta(fasta);
end

Given /the spectra file is (.*)/ do |msms|
	@Tester.msms(msms);
end

Given /I have entered the mass_tol_size (.*)/ do |mass_tol|
	@Tester.mass(mass_tol);
end

Given /I have entered the mass_tol_type (.*)/ do |mass_tol_type|
	@Tester.mass_type(mass_tol_type);
end

Given /I have entered the missed_cleavages (.*)/ do |missed_cleavages|
	@Tester.missed_cleavages(missed_cleavages);
end

When /I run comet and tide/ do
#Specify the other parameters for Tide and Comet
  @Tester.add_args("monoisotopic-precursor=T");
  @Tester.add_args("fragment-mass=mono" );
  @Tester.add_args("decoy-format=none");
  @Tester.add_args("digestion=full-digest");
  @Tester.add_args("min-length=1");
  @Tester.add_args("max-length=100");
  @Tester.add_args("min-mass=100.0");
  @Tester.add_args("max-mass=5000.0");
  @Tester.add_args("use-flanking-peaks=F");
  @Tester.add_args("enzyme=trypsin");
  @Tester.add_args("top-match=500");
  @Tester.add_args("exact-p-value=F");
  @Tester.add_args("min-peaks=1");
 # @Tester.add_args("mods-spec=1MW+15.9949,1QN+0.984");
  @Tester.add_args("max-mods=2");  
  @Tester.add_args("peptide-list=T");

#  Specify the rest of the default parameters for CRUX COMET
  @Tester.add_args("decoy_search = 0");
  @Tester.add_args("mass_type_parent = 1\nmass_type_fragment = 1");
  @Tester.add_args("search_enzyme_number = 1\nnum_enzyme_termini = 2");
  @Tester.add_args("precursor_tolerance_type = 0");
  @Tester.add_args("isotope_error = 0");
  @Tester.add_args("fragment_bin_tol = 1.0005");
  @Tester.add_args("fragment_bin_offset = 0.5");
  @Tester.add_args("theoretical_fragment_ions = 0");
  @Tester.add_args("use_A_ions = 0");
  @Tester.add_args("use_B_ions = 1");
  @Tester.add_args("use_C_ions = 0");
  @Tester.add_args("use_X_ions = 0");
  @Tester.add_args("use_Y_ions = 1");
  @Tester.add_args("use_Z_ions = 0");
  @Tester.add_args("use_NL_ions = 1");
 # @Tester.add_args("variable_mod1 = 15.9949 MW 0 1");
 # @Tester.add_args("variable_mod2 = 0.984 NQ 0 1");
  @Tester.add_args("max_variable_mods_in_peptide = 2");
  @Tester.add_args("use_sparse_matrix = 0");

  @Tester.add_args("minimum_peaks = 1");
  @Tester.add_args("minimum_intensity = 0");
  @Tester.add_args("remove_precursor_peak = 0");
  @Tester.add_args("clip_nterm_methionine = 0");
  @Tester.add_args("max_fragment_charge = 2");
  @Tester.add_args("max_precursor_charge = 5");

  @Tester.add_args("digest_mass_range = 100.0 5000.0");
  @Tester.add_args("num_results = 500");
  @Tester.add_args("output_sqtstream = 0");
  @Tester.add_args("output_sqtfile = 0");
  @Tester.add_args("output_txtfile = 1");
  @Tester.add_args("output_pepxmlfile = 0");
  @Tester.add_args("output_pinxmlfile = 0");
  @Tester.add_args("output_outfiles = 0");
  @Tester.add_args("print_expect_score = 0");
  @Tester.add_args("num_output_lines = 500");
  @Tester.add_args("show_fragment_ions = 0");
  @Tester.add_args("add_C_cysteine = 57.021464");
  
# this must be the last in comet's parameter file.
  @Tester.add_args("[COMET_ENZYME_INFO]");
  @Tester.add_args("0.  No_enzyme              0      -           -]");
  @Tester.add_args("1.  Trypsin                1      KR          P");

  @result = @Tester.runCommand();
#  @result = @Tester.compareCandidatePeptideSets();
#  @result = @Tester.checkResults();
 end
 
Then /Candidate peptides should be identical/ do
	@result.should == true
end


