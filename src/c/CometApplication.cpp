/**
 * \file CometApplication.cpp 
 * \brief Runs hardklor
 *****************************************************************************/
#include "CometSearch/Common.h"
#include "CometSearch/CometSearchManager.h"
#include "ModifiedPeptidesIterator.h"
#include "CarpStreamBuf.h"
#include "CometApplication.h"
#include "DelimitedFileWriter.h"
#include "DelimitedFile.h"

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

   /* Define optional command line arguments */
  const char* option_list[] = {
    "fileroot",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity"
  };

  
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"input spectra","database_name"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */

  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);

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
    const char* var_str, ///< variable modification string
    VarMods& varmods ///< Variable modification structure to set
    ) {
  
  varmods.bBinaryMod = 0;
  varmods.iMaxNumVarModAAPerMod = 0;
  varmods.dVarModMass = 0.0;
  memset(varmods.szVarModChar,0,MAX_VARMOD_AA);


  string temp = var_str;
  vector<string> tokens;
  tokenize(temp, tokens, ' ');
  from_string<double>(varmods.dVarModMass, tokens[0]);
  strcpy(varmods.szVarModChar, tokens[1].c_str());
  from_string<int>(varmods.bBinaryMod, tokens[2]);
  from_string<int>(varmods.iMaxNumVarModAAPerMod, tokens[3]);
}

/**
 * Sets a double range from a string
 */
void CometApplication::getDoubleRange(
    const char* str, ///< range string to set -in
    DoubleRange& doubleRangeParam ///< DoubleRange parameter -out
    ) {
  
  string temp = str;
  vector<string> tokens;
  tokenize(temp, tokens, ' ');
  
  from_string<double>(doubleRangeParam.dStart, tokens[0]);
  from_string<double>(doubleRangeParam.dEnd, tokens[1]);
}

/*
 * sets the integer range from a string
 */
void CometApplication::getIntRange(
    const char* str, ///< range string to set -in
    IntRange& intRangeParam ///< IntRange parameter -out
    ) {
  
  string temp = str;
  vector<string> tokens;
  tokenize(temp, tokens, ' ');
  
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
  strcpy(pInputFile->szFileName, get_string_parameter_pointer("input spectra"));
  
  //VALIDATE
  FILE *fp;
  if ((fp = fopen(pInputFile->szFileName, "r")) == NULL) {
      carp(CARP_FATAL, "Spectra File Not Found:%s", get_string_parameter_pointer("input spectra"));
  }
  fclose(fp);
  
  string scan_range_str = get_string_parameter_pointer("scan_range");
  if (scan_range_str == "0 0") {
    pInputFile->iAnalysisType = AnalysisType_EntireFile;
  } else {
    pInputFile->iAnalysisType = AnalysisType_SpecificScanRange;
    vector<string> tokens;
    tokenize(scan_range_str, tokens, ' ');
    from_string<int>(pInputFile->iFirstScan, tokens[0]);
    from_string<int>(pInputFile->iLastScan, tokens[1]);
  }

  pvInputFiles.push_back(pInputFile);
  //TODO - Comet allows multiple spectra to be searched, add this to crux.
  string basename = make_file_path(getName());
  searchMgr.SetOutputFileBaseName(basename.c_str());
  
  searchMgr.SetParam("database_name", get_string_parameter_pointer("database_name"), get_string_parameter_pointer("database_name"));
  searchMgr.SetParam("decoy_prefix", get_string_parameter_pointer("decoy_prefix"), get_string_parameter_pointer("decoy_prefix"));
  searchMgr.SetParam("output_suffix", get_string_parameter_pointer("output_suffix"), get_string_parameter_pointer("output_suffix"));
  searchMgr.SetParam("nucleotide_reading_frame", get_string_parameter_pointer("nucleotide_reading_frame"), get_int_parameter("nucleotide_reading_frame"));
  searchMgr.SetParam("mass_type_parent", get_string_parameter_pointer("mass_type_parent"), get_int_parameter("mass_type_parent"));
  searchMgr.SetParam("mass_type_fragment", get_string_parameter_pointer("mass_type_fragment"), get_int_parameter("mass_type_fragment"));
  searchMgr.SetParam("show_fragment_ions", get_string_parameter_pointer("show_fragment_ions"), get_int_parameter("show_fragment_ions"));
  searchMgr.SetParam("num_threads", get_string_parameter_pointer("num_threads"), get_int_parameter("num_threads"));
  searchMgr.SetParam("clip_nterm_methionine", get_string_parameter_pointer("clip_nterm_methionine"), get_int_parameter("clip_nterm_methionine"));
  searchMgr.SetParam("theoretical_fragment_ions", get_string_parameter_pointer("theoretical_fragment_ions"), get_int_parameter("theoretical_fragment_ions"));
  searchMgr.SetParam("use_A_ions", get_string_parameter_pointer("use_A_ions"), get_int_parameter("use_A_ions"));
  searchMgr.SetParam("use_B_ions", get_string_parameter_pointer("use_B_ions"), get_int_parameter("use_B_ions"));
  searchMgr.SetParam("use_C_ions", get_string_parameter_pointer("use_C_ions"), get_int_parameter("use_C_ions"));
  searchMgr.SetParam("use_X_ions", get_string_parameter_pointer("use_X_ions"), get_int_parameter("use_X_ions"));
  searchMgr.SetParam("use_Y_ions", get_string_parameter_pointer("use_Y_ions"), get_int_parameter("use_Y_ions"));
  searchMgr.SetParam("use_Z_ions", get_string_parameter_pointer("use_Z_ions"), get_int_parameter("use_Z_ions"));
  searchMgr.SetParam("use_NL_ions", get_string_parameter_pointer("use_NL_ions"), get_int_parameter("use_NL_ions"));
  searchMgr.SetParam("use_sparse_matrix", get_string_parameter_pointer("use_sparse_matrix"), get_int_parameter("use_sparse_matrix"));
  
  if (strncmp(get_string_parameter_pointer("variable_mod1"), "__NULL_STR",10) != 0) {
    calcVarMods(get_string_parameter_pointer("variable_mod1"), varModsParam);
    searchMgr.SetParam("variable_mod1", get_string_parameter_pointer("variable_mod1"), varModsParam );
  }
  if (strncmp(get_string_parameter_pointer("variable_mod2"), "__NULL_STR", 10) != 0) {
    calcVarMods(get_string_parameter_pointer("variable_mod2"), varModsParam);
    searchMgr.SetParam("variable_mod2", get_string_parameter_pointer("variable_mod2"), varModsParam );
  }
  if (strncmp(get_string_parameter_pointer("variable_mod3"), "__NULL_STR", 10) != 0) {
    calcVarMods(get_string_parameter_pointer("variable_mod3"), varModsParam);
    searchMgr.SetParam("variable_mod3", get_string_parameter_pointer("variable_mod3"), varModsParam );
  }
  if (strncmp(get_string_parameter_pointer("variable_mod4"), "__NULL_STR", 10) != 0) {
    calcVarMods(get_string_parameter_pointer("variable_mod4"), varModsParam);
    searchMgr.SetParam("variable_mod4", get_string_parameter_pointer("variable_mod4"), varModsParam );
  }
  if (strncmp(get_string_parameter_pointer("variable_mod5"), "__NULL_STR", 10) != 0) {
    calcVarMods(get_string_parameter_pointer("variable_mod5"), varModsParam);
    searchMgr.SetParam("variable_mod5", get_string_parameter_pointer("variable_mod5"), varModsParam );
  }
  if (strncmp(get_string_parameter_pointer("variable_mod6"), "__NULL_STR", 10) != 0) {
    calcVarMods(get_string_parameter_pointer("variable_mod6"), varModsParam);
    searchMgr.SetParam("variable_mod6", get_string_parameter_pointer("variable_mod6"), varModsParam );
  }

  searchMgr.SetParam("max_variable_mods_in_peptide", get_string_parameter_pointer("max_variable_mods_in_peptide"), get_int_parameter("max_variable_mods_in_peptide"));
  searchMgr.SetParam("fragment_bin_tol", get_string_parameter_pointer("fragment_bin_tol"), get_double_parameter("fragment_bin_tol"));
  searchMgr.SetParam("fragment_bin_offset", get_string_parameter_pointer("fragment_bin_offset"), get_double_parameter("fragment_bin_offset"));
  searchMgr.SetParam("peptide_mass_tolerance", get_string_parameter_pointer("peptide_mass_tolerance"), get_double_parameter("peptide_mass_tolerance"));
  searchMgr.SetParam("precursor_tolerance_type", get_string_parameter_pointer("precursor_tolerance_type"), get_int_parameter("precursor_tolerance_type"));
  searchMgr.SetParam("peptide_mass_units", get_string_parameter_pointer("peptide_mass_units"), get_int_parameter("peptide_mass_units"));
  searchMgr.SetParam("isotope_error", get_string_parameter_pointer("isotope_error"), get_int_parameter("isotope_error"));
  searchMgr.SetParam("num_output_lines", get_string_parameter_pointer("num_output_lines"), get_int_parameter("num_output_lines"));
  searchMgr.SetParam("num_results", get_string_parameter_pointer("num_results"), get_int_parameter("num_results"));
  searchMgr.SetParam("remove_precursor_peak", get_string_parameter_pointer("remove_precursor_peak"), get_int_parameter("remove_precursor_peak"));
  searchMgr.SetParam("remove_precursor_tolerance", get_string_parameter_pointer("remove_precursor_tolerance"), get_double_parameter("remove_precursor_tolerance"));
  
  getDoubleRange(get_string_parameter_pointer("clear_mz_range"), doubleRangeParam );
  searchMgr.SetParam("clear_mz_range", get_string_parameter_pointer("clear_mz_range"), doubleRangeParam );

  searchMgr.SetParam("print_expect_score", get_string_parameter_pointer("print_expect_score"), get_int_parameter("print_expect_score"));
  searchMgr.SetParam("output_sqtstream", "0", 0);
  
  searchMgr.SetParam("output_sqtfile", get_string_parameter_pointer("output_sqtfile"), get_int_parameter("output_sqtfile"));
  searchMgr.SetParam("output_txtfile", get_string_parameter_pointer("output_txtfile"), get_int_parameter("output_txtfile"));
  searchMgr.SetParam("output_pepxmlfile", get_string_parameter_pointer("output_pepxmlfile"), get_int_parameter("output_pepxmlfile"));
  searchMgr.SetParam("output_pinxmlfile", "0", 0); //hardcode to 0
  searchMgr.SetParam("output_outfiles", "0", 0);
  searchMgr.SetParam("skip_researching", get_string_parameter_pointer("skip_researching"), get_int_parameter("skip_researching"));
  searchMgr.SetParam("variable_C_terminus", get_string_parameter_pointer("variable_C_terminus"), get_double_parameter("variable_C_terminus"));
  searchMgr.SetParam("variable_N_terminus", get_string_parameter_pointer("variable_N_terminus"), get_double_parameter("variable_N_terminus"));
  searchMgr.SetParam("variable_C_terminus_distance", get_string_parameter_pointer("variable_C_terminus_distance"), get_int_parameter("variable_C_terminus_distance"));
  searchMgr.SetParam("variable_N_terminus_distance", get_string_parameter_pointer("variable_N_terminus_distance"), get_int_parameter("variable_N_terminus_distance"));
  searchMgr.SetParam("add_Cterm_peptide", get_string_parameter_pointer("add_Cterm_peptide"), get_double_parameter("add_Cterm_peptide"));
  searchMgr.SetParam("add_Nterm_peptide", get_string_parameter_pointer("add_Nterm_peptide"), get_double_parameter("add_Nterm_peptide"));
  searchMgr.SetParam("add_Cterm_protein", get_string_parameter_pointer("add_Cterm_protein"), get_double_parameter("add_Cterm_protein"));
  searchMgr.SetParam("add_Nterm_protein", get_string_parameter_pointer("add_Nterm_protein"), get_double_parameter("add_Nterm_protein"));
  searchMgr.SetParam("add_G_glycine", get_string_parameter_pointer("add_G_glycine"), get_double_parameter("add_G_glycine"));
  searchMgr.SetParam("add_A_alanine", get_string_parameter_pointer("add_A_alanine"), get_double_parameter("add_A_alanine"));
  searchMgr.SetParam("add_S_serine", get_string_parameter_pointer("add_S_serine"), get_double_parameter("add_S_serine"));
  searchMgr.SetParam("add_P_proline", get_string_parameter_pointer("add_P_proline"), get_double_parameter("add_P_proline"));
  searchMgr.SetParam("add_V_valine", get_string_parameter_pointer("add_V_valine"), get_double_parameter("add_V_valine"));
  searchMgr.SetParam("add_T_threonine", get_string_parameter_pointer("add_T_threonine"), get_double_parameter("add_T_threonine"));
  searchMgr.SetParam("add_C_cysteine", get_string_parameter_pointer("add_C_cysteine"), get_double_parameter("add_C_cysteine"));
  searchMgr.SetParam("add_L_leucine", get_string_parameter_pointer("add_L_leucine"), get_double_parameter("add_L_leucine"));
  searchMgr.SetParam("add_I_isoleucine", get_string_parameter_pointer("add_I_isoleucine"), get_double_parameter("add_I_isoleucine"));
  searchMgr.SetParam("add_N_asparagine", get_string_parameter_pointer("add_N_asparagine"), get_double_parameter("add_N_asparagine"));
  searchMgr.SetParam("add_O_ornithine", get_string_parameter_pointer("add_O_ornithine"), get_double_parameter("add_O_ornithine"));
  searchMgr.SetParam("add_D_aspartic_acid", get_string_parameter_pointer("add_D_aspartic_acid"), get_double_parameter("add_D_aspartic_acid"));
  searchMgr.SetParam("add_Q_glutamine", get_string_parameter_pointer("add_Q_glutamine"), get_double_parameter("add_Q_glutamine"));
  searchMgr.SetParam("add_K_lysine", get_string_parameter_pointer("add_K_lysine"), get_double_parameter("add_K_lysine"));
  searchMgr.SetParam("add_E_glutamic_acid", get_string_parameter_pointer("add_E_glutamic_acid"), get_double_parameter("add_E_glutamic_acid"));
  searchMgr.SetParam("add_M_methionine", get_string_parameter_pointer("add_M_methionine"), get_double_parameter("add_M_methionine"));
  searchMgr.SetParam("add_H_histidine", get_string_parameter_pointer("add_H_histidine"), get_double_parameter("add_H_histidine"));
  searchMgr.SetParam("add_F_phenylalanine", get_string_parameter_pointer("add_F_phenylalanine"), get_double_parameter("add_F_phenylalanine"));
  searchMgr.SetParam("add_R_arginine", get_string_parameter_pointer("add_R_arginine"), get_double_parameter("add_R_arginine"));
  searchMgr.SetParam("add_Y_tyrosine", get_string_parameter_pointer("add_Y_tyrosine"), get_double_parameter("add_Y_tyrosine"));
  searchMgr.SetParam("add_W_tryptophan", get_string_parameter_pointer("add_W_tryptophan"), get_double_parameter("add_W_tryptophan"));
  searchMgr.SetParam("add_B_user_amino_acid", get_string_parameter_pointer("add_B_user_amino_acid"), get_double_parameter("add_B_user_amino_acid"));
  searchMgr.SetParam("add_J_user_amino_acid", get_string_parameter_pointer("add_J_user_amino_acid"), get_double_parameter("add_J_user_amino_acid"));
  searchMgr.SetParam("add_U_user_amino_acid", get_string_parameter_pointer("add_U_user_amino_acid"), get_double_parameter("add_U_user_amino_acid"));
  searchMgr.SetParam("add_X_user_amino_acid", get_string_parameter_pointer("add_X_user_amino_acid"), get_double_parameter("add_X_user_amino_acid"));
  searchMgr.SetParam("add_Z_user_amino_acid", get_string_parameter_pointer("add_Z_user_amino_acid"), get_double_parameter("add_Z_user_amino_acid"));
  searchMgr.SetParam("search_enzyme_number", get_string_parameter_pointer("search_enzyme_number"), get_int_parameter("search_enzyme_number"));
  searchMgr.SetParam("sample_enzyme_number", get_string_parameter_pointer("sample_enzyme_number"), get_int_parameter("sample_enzyme_number"));
  searchMgr.SetParam("num_enzyme_termini", get_string_parameter_pointer("num_enzyme_termini"), get_int_parameter("num_enzyme_termini"));
  searchMgr.SetParam("allowed_missed_cleavage", get_string_parameter_pointer("allowed_missed_cleavage"), get_int_parameter("allowed_missed_cleavage"));

  getIntRange(get_string_parameter_pointer("scan_range"), intRangeParam );
  searchMgr.SetParam("scan_range", get_string_parameter_pointer("scan_range"), intRangeParam );

  searchMgr.SetParam("spectrum_batch_size", get_string_parameter_pointer("spectrum_batch_size"), get_int_parameter("spectrum_batch_size"));
  searchMgr.SetParam("minimum_peaks", get_string_parameter_pointer("minimum_peaks"), get_int_parameter("minimum_peaks"));

  getIntRange(get_string_parameter_pointer("precursor_charge"), intRangeParam);
  searchMgr.SetParam("precursor_charge", get_string_parameter_pointer("precursor_charge"), intRangeParam);
  
  searchMgr.SetParam("max_fragment_charge", get_string_parameter_pointer("max_fragment_charge"), get_int_parameter("max_fragment_charge"));
  searchMgr.SetParam("max_precursor_charge", get_string_parameter_pointer("max_precursor_charge"), get_int_parameter("max_precursor_charge"));

  getDoubleRange(get_string_parameter_pointer("digest_mass_range"), doubleRangeParam);
  searchMgr.SetParam("digest_mass_range", get_string_parameter_pointer("digest_mass_range"), doubleRangeParam);
  
  searchMgr.SetParam("ms_level", get_string_parameter_pointer("ms_level"), get_int_parameter("ms_level"));
  searchMgr.SetParam("activation_method", get_string_parameter_pointer("activation_method"), get_string_parameter_pointer("activation_method"));
  searchMgr.SetParam("minimum_intensity", get_string_parameter_pointer("minimum_intensity"), get_double_parameter("minimum_intensity"));
  searchMgr.SetParam("decoy_search", get_string_parameter_pointer("decoy_search"), get_int_parameter("decoy_search"));
  
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
string CometApplication::getName() {
  return "comet";
}

/**
 * \returns the description for CometApplication
 */
string CometApplication::getDescription() {
  return "Search a collection of spectra against a sequence database, "
         "returning a collection of PSMs. This search engine runs directly on "
         "a protein database in FASTA format.";
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool CometApplication::needsOutputDirectory() {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
