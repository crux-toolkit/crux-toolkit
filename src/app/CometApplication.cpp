/**
 * \file CometApplication.cpp 
 * \brief Runs hardklor
 *****************************************************************************/
#include "CometSearch/Common.h"
#include "CometSearch/CometSearchManager.h"
#include "model/ModifiedPeptidesIterator.h"
#include "util/CarpStreamBuf.h"
#include "util/StringUtils.h"
#include "CometApplication.h"
#include "io/DelimitedFileWriter.h"
#include "io/DelimitedFile.h"

using namespace std;

/**
 * \returns a blank CometApplication object
 */
CometApplication::CometApplication() {

}

/**
 * Destructor
 */
CometApplication::~CometApplication() {
}

/**
 * main method for CometApplication
 */
int CometApplication::main(int argc, char** argv) {
  initialize(argc, argv);

  /* Re-route stderr to log file */
  CarpStreamBuf buffer;
  streambuf* old = std::cerr.rdbuf();
  std::cerr.rdbuf(&buffer);

  /* set Parmeters */
  vector<InputFileInfo*> pv_input_files;
  CometSearchManager comet_search_mgr;
  setCometParameters(pv_input_files, comet_search_mgr);
  comet_search_mgr.AddInputFiles(pv_input_files);
  
  /* Run search */
  comet_search_mgr.DoSearch();

  /* Recover stderr */
  std::cerr.rdbuf(old);

  return 0;
}

/**
 * Extracts the variable modification info from the string
 */
void CometApplication::calcVarMods(
    const string& var_str, ///< variable modification string
    VarMods& varmods ///< Variable modification structure to set
    ) {
  
  varmods.bBinaryMod = 0;
  varmods.iMaxNumVarModAAPerMod = 0;
  varmods.dVarModMass = 0.0;
  memset(varmods.szVarModChar,0,MAX_VARMOD_AA);

  vector<string> tokens = StringUtils::Split(var_str, ' ');
  from_string<double>(varmods.dVarModMass, tokens[0]);
  strcpy(varmods.szVarModChar, tokens[1].c_str());
  from_string<int>(varmods.bBinaryMod, tokens[2]);
  from_string<int>(varmods.iMaxNumVarModAAPerMod, tokens[3]);
}

/**
 * Sets a double range from a string
 */
void CometApplication::getDoubleRange(
    const string& str, ///< range string to set -in
    DoubleRange& doubleRangeParam ///< DoubleRange parameter -out
    ) {
  vector<string> tokens = StringUtils::Split(str, ' ');
  
  from_string<double>(doubleRangeParam.dStart, tokens[0]);
  from_string<double>(doubleRangeParam.dEnd, tokens[1]);
}

/*
 * sets the integer range from a string
 */
void CometApplication::getIntRange(
    const string& str, ///< range string to set -in
    IntRange& intRangeParam ///< IntRange parameter -out
    ) {
  vector<string> tokens = StringUtils::Split(str, ' ');
  
  from_string<int>(intRangeParam.iStart, tokens[0]);
  from_string<int>(intRangeParam.iEnd, tokens[1]);
}

/**
 * Sets the parameters for the Comet application using the crux parameters
 */
void CometApplication::setCometParameters(
  vector<InputFileInfo*> &pvInputFiles, ///<vector of input spectra files
  CometSearchManager& searchMgr ///< SearchManager to set the parameters
  ) {
  
  VarMods varModsParam;
  IntRange intRangeParam;
  DoubleRange doubleRangeParam;

  
  InputFileInfo *pInputFile = new InputFileInfo();
  pInputFile->iAnalysisType = 0;
  strcpy(pInputFile->szFileName, get_string_parameter("input spectra").c_str());
  
  //VALIDATE
  FILE *fp;
  if ((fp = fopen(pInputFile->szFileName, "r")) == NULL) {
      carp(CARP_FATAL, "Spectra File Not Found:%s", get_string_parameter("input spectra").c_str());
  }
  fclose(fp);
  
  string scan_range_str = get_string_parameter("scan_range");
  if (scan_range_str == "0 0") {
    pInputFile->iAnalysisType = AnalysisType_EntireFile;
  } else {
    pInputFile->iAnalysisType = AnalysisType_SpecificScanRange;
    vector<string> tokens = StringUtils::Split(scan_range_str, ' ');
    from_string<int>(pInputFile->iFirstScan, tokens[0]);
    from_string<int>(pInputFile->iLastScan, tokens[1]);
  }

  pvInputFiles.push_back(pInputFile);
  //TODO - Comet allows multiple spectra to be searched, add this to crux.
  string basename = make_file_path(getName());
  searchMgr.SetOutputFileBaseName(basename.c_str());
  
  searchMgr.SetParam("database_name", get_string_parameter("database name"), get_string_parameter("database name"));
  searchMgr.SetParam("decoy_prefix", get_string_parameter("decoy_prefix"), get_string_parameter("decoy_prefix"));
  searchMgr.SetParam("output_suffix", get_string_parameter("output_suffix"), get_string_parameter("output_suffix"));
  searchMgr.SetParam("nucleotide_reading_frame", get_string_parameter("nucleotide_reading_frame"), get_int_parameter("nucleotide_reading_frame"));
  searchMgr.SetParam("mass_type_parent", get_string_parameter("mass_type_parent"), get_int_parameter("mass_type_parent"));
  searchMgr.SetParam("mass_type_fragment", get_string_parameter("mass_type_fragment"), get_int_parameter("mass_type_fragment"));
  searchMgr.SetParam("show_fragment_ions", get_string_parameter("show_fragment_ions"), get_int_parameter("show_fragment_ions"));
  searchMgr.SetParam("num_threads", get_string_parameter("num_threads"), get_int_parameter("num_threads"));
  searchMgr.SetParam("clip_nterm_methionine", get_string_parameter("clip_nterm_methionine"), get_int_parameter("clip_nterm_methionine"));
  searchMgr.SetParam("theoretical_fragment_ions", get_string_parameter("theoretical_fragment_ions"), get_int_parameter("theoretical_fragment_ions"));
  searchMgr.SetParam("use_A_ions", get_string_parameter("use_A_ions"), get_int_parameter("use_A_ions"));
  searchMgr.SetParam("use_B_ions", get_string_parameter("use_B_ions"), get_int_parameter("use_B_ions"));
  searchMgr.SetParam("use_C_ions", get_string_parameter("use_C_ions"), get_int_parameter("use_C_ions"));
  searchMgr.SetParam("use_X_ions", get_string_parameter("use_X_ions"), get_int_parameter("use_X_ions"));
  searchMgr.SetParam("use_Y_ions", get_string_parameter("use_Y_ions"), get_int_parameter("use_Y_ions"));
  searchMgr.SetParam("use_Z_ions", get_string_parameter("use_Z_ions"), get_int_parameter("use_Z_ions"));
  searchMgr.SetParam("use_NL_ions", get_string_parameter("use_NL_ions"), get_int_parameter("use_NL_ions"));
  searchMgr.SetParam("use_sparse_matrix", get_string_parameter("use_sparse_matrix"), get_int_parameter("use_sparse_matrix"));
  
  if (!get_string_parameter("variable_mod01").empty()) {
    calcVarMods(get_string_parameter("variable_mod01"), varModsParam);
    searchMgr.SetParam("variable_mod01", get_string_parameter("variable_mod01"), varModsParam );
  }
  if (!get_string_parameter("variable_mod02").empty()) {
    calcVarMods(get_string_parameter("variable_mod02"), varModsParam);
    searchMgr.SetParam("variable_mod02", get_string_parameter("variable_mod02"), varModsParam );
  }
  if (!get_string_parameter("variable_mod03").empty()) {
    calcVarMods(get_string_parameter("variable_mod03"), varModsParam);
    searchMgr.SetParam("variable_mod03", get_string_parameter("variable_mod03"), varModsParam );
  }
  if (!get_string_parameter("variable_mod04").empty()) {
    calcVarMods(get_string_parameter("variable_mod04"), varModsParam);
    searchMgr.SetParam("variable_mod04", get_string_parameter("variable_mod04"), varModsParam );
  }
  if (!get_string_parameter("variable_mod05").empty()) {
    calcVarMods(get_string_parameter("variable_mod05"), varModsParam);
    searchMgr.SetParam("variable_mod05", get_string_parameter("variable_mod05"), varModsParam );
  }
  if (!get_string_parameter("variable_mod06").empty()) {
    calcVarMods(get_string_parameter("variable_mod06"), varModsParam);
    searchMgr.SetParam("variable_mod06", get_string_parameter("variable_mod06"), varModsParam );
  }
  if (!get_string_parameter("variable_mod07").empty()) {
    calcVarMods(get_string_parameter("variable_mod07"), varModsParam);
    searchMgr.SetParam("variable_mod07", get_string_parameter("variable_mod07"), varModsParam );
  }
  if (!get_string_parameter("variable_mod08").empty()) {
    calcVarMods(get_string_parameter("variable_mod08"), varModsParam);
    searchMgr.SetParam("variable_mod08", get_string_parameter("variable_mod08"), varModsParam );
  }
  if (!get_string_parameter("variable_mod09").empty()) {
    calcVarMods(get_string_parameter("variable_mod09"), varModsParam);
    searchMgr.SetParam("variable_mod09", get_string_parameter("variable_mod09"), varModsParam );
  }

  searchMgr.SetParam("max_variable_mods_in_peptide", get_string_parameter("max_variable_mods_in_peptide"), get_int_parameter("max_variable_mods_in_peptide"));
  searchMgr.SetParam("fragment_bin_tol", get_string_parameter("fragment_bin_tol"), get_double_parameter("fragment_bin_tol"));
  searchMgr.SetParam("fragment_bin_offset", get_string_parameter("fragment_bin_offset"), get_double_parameter("fragment_bin_offset"));
  searchMgr.SetParam("peptide_mass_tolerance", get_string_parameter("peptide_mass_tolerance"), get_double_parameter("peptide_mass_tolerance"));
  searchMgr.SetParam("peptide_mass_units", get_string_parameter("peptide_mass_units"), get_int_parameter("peptide_mass_units"));
  searchMgr.SetParam("isotope_error", get_string_parameter("isotope_error"), get_int_parameter("isotope_error"));
  searchMgr.SetParam("num_output_lines", get_string_parameter("num_output_lines"), get_int_parameter("num_output_lines"));
  searchMgr.SetParam("num_results", get_string_parameter("num_results"), get_int_parameter("num_results"));
  searchMgr.SetParam("remove_precursor_peak", get_string_parameter("remove_precursor_peak"), get_int_parameter("remove_precursor_peak"));
  searchMgr.SetParam("remove_precursor_tolerance", get_string_parameter("remove_precursor_tolerance"), get_double_parameter("remove_precursor_tolerance"));
  
  getDoubleRange(get_string_parameter("clear_mz_range"), doubleRangeParam );
  searchMgr.SetParam("clear_mz_range", get_string_parameter("clear_mz_range"), doubleRangeParam );

  searchMgr.SetParam("print_expect_score", get_string_parameter("print_expect_score"), get_int_parameter("print_expect_score"));
  searchMgr.SetParam("output_sqtstream", "0", 0);
  
  searchMgr.SetParam("override_charge", get_string_parameter("override_charge"), get_int_parameter("override_charge"));

  searchMgr.SetParam("require_variable_mod", get_string_parameter("require_variable_mod"), get_int_parameter("require_variable_mod"));
  
  searchMgr.SetParam("output_sqtfile", get_string_parameter("output_sqtfile"), get_int_parameter("output_sqtfile"));
  searchMgr.SetParam("output_txtfile", get_string_parameter("output_txtfile"), get_int_parameter("output_txtfile"));
  searchMgr.SetParam("output_pepxmlfile", get_string_parameter("output_pepxmlfile"), get_int_parameter("output_pepxmlfile"));
  searchMgr.SetParam("output_percolatorfile", get_string_parameter("output_percolatorfile"), get_int_parameter("output_percolatorfile"));
  searchMgr.SetParam("output_outfiles", "0", 0);
  searchMgr.SetParam("skip_researching", get_string_parameter("skip_researching"), get_int_parameter("skip_researching"));
  searchMgr.SetParam("add_Cterm_peptide", get_string_parameter("add_Cterm_peptide"), get_double_parameter("add_Cterm_peptide"));
  searchMgr.SetParam("add_Nterm_peptide", get_string_parameter("add_Nterm_peptide"), get_double_parameter("add_Nterm_peptide"));
  searchMgr.SetParam("add_Cterm_protein", get_string_parameter("add_Cterm_protein"), get_double_parameter("add_Cterm_protein"));
  searchMgr.SetParam("add_Nterm_protein", get_string_parameter("add_Nterm_protein"), get_double_parameter("add_Nterm_protein"));
  searchMgr.SetParam("add_G_glycine", get_string_parameter("add_G_glycine"), get_double_parameter("add_G_glycine"));
  searchMgr.SetParam("add_A_alanine", get_string_parameter("add_A_alanine"), get_double_parameter("add_A_alanine"));
  searchMgr.SetParam("add_S_serine", get_string_parameter("add_S_serine"), get_double_parameter("add_S_serine"));
  searchMgr.SetParam("add_P_proline", get_string_parameter("add_P_proline"), get_double_parameter("add_P_proline"));
  searchMgr.SetParam("add_V_valine", get_string_parameter("add_V_valine"), get_double_parameter("add_V_valine"));
  searchMgr.SetParam("add_T_threonine", get_string_parameter("add_T_threonine"), get_double_parameter("add_T_threonine"));
  searchMgr.SetParam("add_C_cysteine", get_string_parameter("add_C_cysteine"), get_double_parameter("add_C_cysteine"));
  searchMgr.SetParam("add_L_leucine", get_string_parameter("add_L_leucine"), get_double_parameter("add_L_leucine"));
  searchMgr.SetParam("add_I_isoleucine", get_string_parameter("add_I_isoleucine"), get_double_parameter("add_I_isoleucine"));
  searchMgr.SetParam("add_N_asparagine", get_string_parameter("add_N_asparagine"), get_double_parameter("add_N_asparagine"));
  searchMgr.SetParam("add_O_ornithine", get_string_parameter("add_O_ornithine"), get_double_parameter("add_O_ornithine"));
  searchMgr.SetParam("add_D_aspartic_acid", get_string_parameter("add_D_aspartic_acid"), get_double_parameter("add_D_aspartic_acid"));
  searchMgr.SetParam("add_Q_glutamine", get_string_parameter("add_Q_glutamine"), get_double_parameter("add_Q_glutamine"));
  searchMgr.SetParam("add_K_lysine", get_string_parameter("add_K_lysine"), get_double_parameter("add_K_lysine"));
  searchMgr.SetParam("add_E_glutamic_acid", get_string_parameter("add_E_glutamic_acid"), get_double_parameter("add_E_glutamic_acid"));
  searchMgr.SetParam("add_M_methionine", get_string_parameter("add_M_methionine"), get_double_parameter("add_M_methionine"));
  searchMgr.SetParam("add_H_histidine", get_string_parameter("add_H_histidine"), get_double_parameter("add_H_histidine"));
  searchMgr.SetParam("add_F_phenylalanine", get_string_parameter("add_F_phenylalanine"), get_double_parameter("add_F_phenylalanine"));
  searchMgr.SetParam("add_R_arginine", get_string_parameter("add_R_arginine"), get_double_parameter("add_R_arginine"));
  searchMgr.SetParam("add_Y_tyrosine", get_string_parameter("add_Y_tyrosine"), get_double_parameter("add_Y_tyrosine"));
  searchMgr.SetParam("add_W_tryptophan", get_string_parameter("add_W_tryptophan"), get_double_parameter("add_W_tryptophan"));
  searchMgr.SetParam("add_B_user_amino_acid", get_string_parameter("add_B_user_amino_acid"), get_double_parameter("add_B_user_amino_acid"));
  searchMgr.SetParam("add_J_user_amino_acid", get_string_parameter("add_J_user_amino_acid"), get_double_parameter("add_J_user_amino_acid"));
  searchMgr.SetParam("add_U_user_amino_acid", get_string_parameter("add_U_user_amino_acid"), get_double_parameter("add_U_user_amino_acid"));
  searchMgr.SetParam("add_X_user_amino_acid", get_string_parameter("add_X_user_amino_acid"), get_double_parameter("add_X_user_amino_acid"));
  searchMgr.SetParam("add_Z_user_amino_acid", get_string_parameter("add_Z_user_amino_acid"), get_double_parameter("add_Z_user_amino_acid"));
  searchMgr.SetParam("search_enzyme_number", get_string_parameter("search_enzyme_number"), get_int_parameter("search_enzyme_number"));
  searchMgr.SetParam("sample_enzyme_number", get_string_parameter("sample_enzyme_number"), get_int_parameter("sample_enzyme_number"));
  searchMgr.SetParam("num_enzyme_termini", get_string_parameter("num_enzyme_termini"), get_int_parameter("num_enzyme_termini"));
  searchMgr.SetParam("allowed_missed_cleavage", get_string_parameter("allowed_missed_cleavage"), get_int_parameter("allowed_missed_cleavage"));

  getIntRange(get_string_parameter("scan_range"), intRangeParam );
  searchMgr.SetParam("scan_range", get_string_parameter("scan_range"), intRangeParam );

  searchMgr.SetParam("spectrum_batch_size", get_string_parameter("spectrum_batch_size"), get_int_parameter("spectrum_batch_size"));
  searchMgr.SetParam("minimum_peaks", get_string_parameter("minimum_peaks"), get_int_parameter("minimum_peaks"));

  getIntRange(get_string_parameter("precursor_charge"), intRangeParam);
  searchMgr.SetParam("precursor_charge", get_string_parameter("precursor_charge"), intRangeParam);
  
  searchMgr.SetParam("max_fragment_charge", get_string_parameter("max_fragment_charge"), get_int_parameter("max_fragment_charge"));
  searchMgr.SetParam("max_precursor_charge", get_string_parameter("max_precursor_charge"), get_int_parameter("max_precursor_charge"));

  getDoubleRange(get_string_parameter("digest_mass_range"), doubleRangeParam);
  searchMgr.SetParam("digest_mass_range", get_string_parameter("digest_mass_range"), doubleRangeParam);
  
  searchMgr.SetParam("ms_level", get_string_parameter("ms_level"), get_int_parameter("ms_level"));
  searchMgr.SetParam("activation_method", get_string_parameter("activation_method"), get_string_parameter("activation_method"));
  searchMgr.SetParam("minimum_intensity", get_string_parameter("minimum_intensity"), get_double_parameter("minimum_intensity"));
  searchMgr.SetParam("decoy_search", get_string_parameter("decoy_search"), get_int_parameter("decoy_search"));
  
  EnzymeInfo enzymeInformation;
  double temp;
  int search_enzyme_number = get_int_parameter("search_enzyme_number");

  if (search_enzyme_number >=0 && search_enzyme_number < get_comet_enzyme_info_lines().size()) {
    const char* szParamBuf = get_comet_enzyme_info_lines()[search_enzyme_number].c_str();
    sscanf(szParamBuf, "%lf %48s %d %20s %20s\n",
      &temp, 
      enzymeInformation.szSearchEnzymeName, 
      &enzymeInformation.iSearchEnzymeOffSet, 
      enzymeInformation.szSearchEnzymeBreakAA, 
      enzymeInformation.szSearchEnzymeNoBreakAA);
  } else {
    carp(CARP_FATAL, "search_enzyme_number=%d out of range (0-%d)", 
      search_enzyme_number, (get_comet_enzyme_info_lines().size()-1));
  }

  int sample_enzyme_number = get_int_parameter("sample_enzyme_number");
  if (sample_enzyme_number >= 0 && sample_enzyme_number < get_comet_enzyme_info_lines().size()) {
    const char* szParamBuf = get_comet_enzyme_info_lines()[sample_enzyme_number].c_str();
    sscanf(szParamBuf, "%lf %48s %d %20s %20s\n",
      &temp, 
      enzymeInformation.szSampleEnzymeName, 
      &enzymeInformation.iSampleEnzymeOffSet, 
      enzymeInformation.szSampleEnzymeBreakAA, 
      enzymeInformation.szSampleEnzymeNoBreakAA);
  } else {
    carp(CARP_FATAL, "sample_enzyme_number=%d out of range (0-%d)", 
      sample_enzyme_number, (get_comet_enzyme_info_lines().size()-1));
  }
  enzymeInformation.iAllowedMissedCleavage = get_int_parameter("allowed_missed_cleavage");
  searchMgr.SetParam("[COMET_ENZYME_INFO]", "TODO", enzymeInformation);
  
}


/**
 * \returns the command name for CometApplication
 */
string CometApplication::getName() const {
  return "comet";
}

/**
 * \returns the description for CometApplication
 */
string CometApplication::getDescription() const {
  return
    "[[nohtml:Search a collection of spectra against a sequence database, "
    "returning a collection of PSMs. This search engine runs directly on a "
    "protein database in FASTA format.]]"
    "[[html:<p>This command searches a protein database with a set of spectra, "
    "assigning peptide sequences to the observed spectra. This search engine "
    "was developed by Jimmy Eng at the University of Washington Proteomics "
    "Resource.</p><p>Although its history goes back two decades, the Comet "
    "search engine was first made publicly available in August 2012 on "
    "SourceForge. Comet is multithreaded and supports multiple input and "
    "output formats.</p><blockquote><a href=\"http://onlinelibrary.wiley.com/"
    "doi/10.1002/pmic.201200439/abstract\">&quot;Comet: an open source tandem "
    "mass spectrometry sequence database search tool.&quot;</a> Eng JK, Jahan "
    "TA, Hoopmann MR. <em>Proteomics</em>. 2012 Nov 12. doi: "
    "10.1002/pmic201200439</blockquote>]]";
}

/**
 * \returns the command arguments
 */
vector<string> CometApplication::getArgs() const {
  string arr[] = {
    "input spectra",
    "database name"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> CometApplication::getOptions() const {
  string arr[] = {
    "fileroot",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity",
    "decoy_search",
    "num_threads",
    "output_suffix",
    "peptide_mass_tolerance",
    "peptide_mass_units",
    "mass_type_parent",
    "mass_type_fragment",
    "isotope_error",
    "search_enzyme_number",
    "num_enzyme_termini",
    "allowed_missed_cleavage",
    "fragment_bin_tol",
    "fragment_bin_offset",
    "theoretical_fragment_ions",
    "use_A_ions",
    "use_B_ions",
    "use_C_ions",
    "use_X_ions",
    "use_Y_ions",
    "use_Z_ions",
    "use_NL_ions",
    "use_sparse_matrix",
    "output_sqtfile",
    "output_pepxmlfile",
    "output_percolatorfile",
    "output_txtfile",
    "output_outfiles",
    "print_expect_score",
    "num_output_lines",
    "show_fragment_ions",
    "sample_enzyme_number",
    "scan_range",
    "precursor_charge",
    "ms_level",
    "activation_method",
    "digest_mass_range",
    "num_results",
    "skip_researching",
    "max_fragment_charge",
    "max_precursor_charge",
    "nucleotide_reading_frame",
    "clip_nterm_methionine",
    "spectrum_batch_size",
    "minimum_peaks",
    "minimum_intensity",
    "remove_precursor_peak",
    "remove_precursor_tolerance",
    "clear_mz_range",
    "variable_mod01",
    "variable_mod02",
    "variable_mod03",
    "variable_mod04",
    "variable_mod05",
    "variable_mod06",
    "variable_mod07",
    "variable_mod08",
    "variable_mod09",
    "require_variable_mod",
    "max_variable_mods_in_peptide",
    "override_charge",
    "add_Cterm_peptide",
    "add_Nterm_peptide",
    "add_Cterm_protein",
    "add_Nterm_protein",
    "add_A_alanine",
    "add_B_user_amino_acid",
    "add_C_cysteine",
    "add_D_aspartic_acid",
    "add_E_glutamic_acid",
    "add_F_phenylalanine",
    "add_G_glycine",
    "add_H_histidine",
    "add_I_isoleucine",
    "add_J_user_amino_acid",
    "add_K_lysine",
    "add_L_leucine",
    "add_M_methionine",
    "add_N_asparagine",
    "add_O_ornithine",
    "add_P_proline",
    "add_Q_glutamine",
    "add_R_arginine",
    "add_S_serine",
    "add_T_threonine",
    "add_U_user_amino_acid",
    "add_V_valine",
    "add_W_tryptophan",
    "add_X_user_amino_acid",
    "add_Y_tyrosine",
    "add_Z_user_amino_acid"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
map<string, string> CometApplication::getOutputs() const {
  map<string, string> outputs;
  outputs["comet.params.txt"] =
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs.";
  outputs["comet.target.txt"] =
    "a tab-delimited text file containing the target PSMs. See <a href=\"txt-format.html\">"
    "txt file format</a> for a list of the fields.";
  outputs["comet.log.txt"] =
    "a log file containing a copy of all messages that were printed to "
    "standard error.";
  return outputs;
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool CometApplication::needsOutputDirectory() const {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
