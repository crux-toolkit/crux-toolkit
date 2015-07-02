/**
 * \file CometApplication.cpp 
 * \brief Runs comet
 *****************************************************************************/
#include "CometSearch/Common.h"
#include "CometSearch/CometSearchManager.h"
#include "model/ModifiedPeptidesIterator.h"
#include "util/AminoAcidUtil.h"
#include "util/CarpStreamBuf.h"
#include "util/FileUtils.h"
#include "util/Params.h"
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
  return main(Params::GetStrings("input spectra"));
}

int CometApplication::main(const vector<string>& input_files) {
  for (vector<string>::const_iterator i = input_files.begin(); i != input_files.end(); i++) {
    if (!StringUtils::IEndsWith(*i, ".mzXML") && !StringUtils::IEndsWith(*i, ".mzML") &&
        !StringUtils::IEndsWith(*i, ".mz5") && !StringUtils::IEndsWith(*i, ".mzXML.gz") &&
        !StringUtils::IEndsWith(*i, ".mzML.gz") && !StringUtils::IEndsWith(*i, ".raw") &&
        !StringUtils::IEndsWith(*i, ".ms2") && !StringUtils::IEndsWith(*i, ".cms2")) {
      carp(CARP_FATAL, "The format of file '%s' is not supported.", i->c_str());
    }
  }

  /* Re-route stderr to log file */
  CarpStreamBuf buffer;
  streambuf* old = std::cerr.rdbuf();
  std::cerr.rdbuf(&buffer);

  /* set Parmeters */
  vector<InputFileInfo*> pv_input_files;
  CometSearchManager comet_search_mgr;
  setCometParameters(input_files, pv_input_files, comet_search_mgr);
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
  const vector<string>& spec_files,
  vector<InputFileInfo*> &pvInputFiles, ///<vector of input spectra files
  CometSearchManager& searchMgr ///< SearchManager to set the parameters
  ) {
  VarMods varModsParam;
  IntRange intRangeParam;
  DoubleRange doubleRangeParam;

  string scan_range_str = Params::GetString("scan_range");
  int analysis_type = AnalysisType_EntireFile;
  int first_scan, last_scan;
  if (scan_range_str != "0 0") {
    analysis_type = AnalysisType_SpecificScanRange;
    vector<string> tokens = StringUtils::Split(scan_range_str, ' ');
    first_scan = StringUtils::FromString<int>(tokens[0]);
    last_scan = StringUtils::FromString<int>(tokens[1]);
  }

  for (vector<string>::const_iterator i = spec_files.begin(); i != spec_files.end(); i++) {
    if (!FileUtils::Exists(*i)) {
      carp(CARP_FATAL, "Spectra File Not Found:%s", i->c_str());
    }
    InputFileInfo* pInputFile = new InputFileInfo();
    strcpy(pInputFile->szFileName, i->c_str());
    pInputFile->iAnalysisType = analysis_type;
    if (analysis_type == AnalysisType_SpecificScanRange) {
      pInputFile->iFirstScan = first_scan;
      pInputFile->iLastScan = last_scan;
    }
    string basename = make_file_path(getFileStem());
    if (spec_files.size() > 1) {
      basename += "." + FileUtils::Stem(*i);
    }
    strcpy(pInputFile->szBaseName, basename.c_str());
    pvInputFiles.push_back(pInputFile);
  }

  searchMgr.SetParam("database_name", Params::GetString("database name"), Params::GetString("database name"));
  searchMgr.SetParam("decoy_prefix", Params::GetString("decoy_prefix"), Params::GetString("decoy_prefix"));
  searchMgr.SetParam("output_suffix", Params::GetString("output_suffix"), Params::GetString("output_suffix"));
  searchMgr.SetParam("nucleotide_reading_frame", Params::GetString("nucleotide_reading_frame"), Params::GetInt("nucleotide_reading_frame"));
  searchMgr.SetParam("mass_type_parent", Params::GetString("mass_type_parent"), Params::GetInt("mass_type_parent"));
  searchMgr.SetParam("mass_type_fragment", Params::GetString("mass_type_fragment"), Params::GetInt("mass_type_fragment"));
  searchMgr.SetParam("show_fragment_ions", Params::GetString("show_fragment_ions"), Params::GetInt("show_fragment_ions"));
  searchMgr.SetParam("num_threads", Params::GetString("num_threads"), Params::GetInt("num_threads"));
  searchMgr.SetParam("clip_nterm_methionine", Params::GetString("clip_nterm_methionine"), Params::GetInt("clip_nterm_methionine"));
  searchMgr.SetParam("theoretical_fragment_ions", Params::GetString("theoretical_fragment_ions"), Params::GetInt("theoretical_fragment_ions"));
  searchMgr.SetParam("use_A_ions", Params::GetString("use_A_ions"), Params::GetInt("use_A_ions"));
  searchMgr.SetParam("use_B_ions", Params::GetString("use_B_ions"), Params::GetInt("use_B_ions"));
  searchMgr.SetParam("use_C_ions", Params::GetString("use_C_ions"), Params::GetInt("use_C_ions"));
  searchMgr.SetParam("use_X_ions", Params::GetString("use_X_ions"), Params::GetInt("use_X_ions"));
  searchMgr.SetParam("use_Y_ions", Params::GetString("use_Y_ions"), Params::GetInt("use_Y_ions"));
  searchMgr.SetParam("use_Z_ions", Params::GetString("use_Z_ions"), Params::GetInt("use_Z_ions"));
  searchMgr.SetParam("use_NL_ions", Params::GetString("use_NL_ions"), Params::GetInt("use_NL_ions"));
  searchMgr.SetParam("use_sparse_matrix", Params::GetString("use_sparse_matrix"), Params::GetInt("use_sparse_matrix"));

  for (int i = 1; i <= 9; i++) {
    string param = "variable_mod0" + StringUtils::ToString(i);
    string paramValue = Params::GetString(param);
    if (!paramValue.empty()) {
      calcVarMods(paramValue, varModsParam);
      searchMgr.SetParam(param, paramValue, varModsParam);
    }
  }

  searchMgr.SetParam("max_variable_mods_in_peptide", Params::GetString("max_variable_mods_in_peptide"), Params::GetInt("max_variable_mods_in_peptide"));
  searchMgr.SetParam("fragment_bin_tol", Params::GetString("fragment_bin_tol"), Params::GetDouble("fragment_bin_tol"));
  searchMgr.SetParam("fragment_bin_offset", Params::GetString("fragment_bin_offset"), Params::GetDouble("fragment_bin_offset"));
  searchMgr.SetParam("peptide_mass_tolerance", Params::GetString("peptide_mass_tolerance"), Params::GetDouble("peptide_mass_tolerance"));
  searchMgr.SetParam("peptide_mass_units", Params::GetString("peptide_mass_units"), Params::GetInt("peptide_mass_units"));
  searchMgr.SetParam("isotope_error", Params::GetString("isotope_error"), Params::GetInt("isotope_error"));
  searchMgr.SetParam("num_output_lines", Params::GetString("num_output_lines"), Params::GetInt("num_output_lines"));
  searchMgr.SetParam("num_results", Params::GetString("num_results"), Params::GetInt("num_results"));
  searchMgr.SetParam("remove_precursor_peak", Params::GetString("remove_precursor_peak"), Params::GetInt("remove_precursor_peak"));
  searchMgr.SetParam("remove_precursor_tolerance", Params::GetString("remove_precursor_tolerance"), Params::GetDouble("remove_precursor_tolerance"));
  
  getDoubleRange(Params::GetString("clear_mz_range"), doubleRangeParam );
  searchMgr.SetParam("clear_mz_range", Params::GetString("clear_mz_range"), doubleRangeParam );

  searchMgr.SetParam("print_expect_score", Params::GetString("print_expect_score"), Params::GetInt("print_expect_score"));
  searchMgr.SetParam("output_sqtstream", "0", 0);
  
  searchMgr.SetParam("override_charge", Params::GetString("override_charge"), Params::GetInt("override_charge"));

  searchMgr.SetParam("require_variable_mod", Params::GetString("require_variable_mod"), Params::GetInt("require_variable_mod"));
  
  searchMgr.SetParam("output_sqtfile", Params::GetString("output_sqtfile"), Params::GetInt("output_sqtfile"));
  searchMgr.SetParam("output_txtfile", Params::GetString("output_txtfile"), Params::GetInt("output_txtfile"));
  searchMgr.SetParam("output_pepxmlfile", Params::GetString("output_pepxmlfile"), Params::GetInt("output_pepxmlfile"));
  searchMgr.SetParam("output_percolatorfile", Params::GetString("output_percolatorfile"), Params::GetInt("output_percolatorfile"));
  searchMgr.SetParam("output_outfiles", "0", 0);
  searchMgr.SetParam("skip_researching", Params::GetString("skip_researching"), Params::GetInt("skip_researching"));
  searchMgr.SetParam("add_Cterm_peptide", Params::GetString("add_Cterm_peptide"), Params::GetDouble("add_Cterm_peptide"));
  searchMgr.SetParam("add_Nterm_peptide", Params::GetString("add_Nterm_peptide"), Params::GetDouble("add_Nterm_peptide"));
  searchMgr.SetParam("add_Cterm_protein", Params::GetString("add_Cterm_protein"), Params::GetDouble("add_Cterm_protein"));
  searchMgr.SetParam("add_Nterm_protein", Params::GetString("add_Nterm_protein"), Params::GetDouble("add_Nterm_protein"));
  for (char c = 'A'; c <= 'Z'; c++) {
    string aaName = AminoAcidUtil::GetName(c);
    aaName = aaName.empty() ? "user_amino_acid" : StringUtils::Replace(aaName, " ", "_");
    string param = "add_" + string(1, c) + "_" + aaName;
    searchMgr.SetParam(param, Params::GetString(param), Params::GetDouble(param));
  }
  searchMgr.SetParam("search_enzyme_number", Params::GetString("search_enzyme_number"), Params::GetInt("search_enzyme_number"));
  searchMgr.SetParam("sample_enzyme_number", Params::GetString("sample_enzyme_number"), Params::GetInt("sample_enzyme_number"));
  searchMgr.SetParam("num_enzyme_termini", Params::GetString("num_enzyme_termini"), Params::GetInt("num_enzyme_termini"));
  searchMgr.SetParam("allowed_missed_cleavage", Params::GetString("allowed_missed_cleavage"), Params::GetInt("allowed_missed_cleavage"));

  getIntRange(Params::GetString("scan_range"), intRangeParam );
  searchMgr.SetParam("scan_range", Params::GetString("scan_range"), intRangeParam );

  searchMgr.SetParam("spectrum_batch_size", Params::GetString("spectrum_batch_size"), Params::GetInt("spectrum_batch_size"));
  searchMgr.SetParam("minimum_peaks", Params::GetString("minimum_peaks"), Params::GetInt("minimum_peaks"));

  getIntRange(Params::GetString("precursor_charge"), intRangeParam);
  searchMgr.SetParam("precursor_charge", Params::GetString("precursor_charge"), intRangeParam);
  
  searchMgr.SetParam("max_fragment_charge", Params::GetString("max_fragment_charge"), Params::GetInt("max_fragment_charge"));
  searchMgr.SetParam("max_precursor_charge", Params::GetString("max_precursor_charge"), Params::GetInt("max_precursor_charge"));

  getDoubleRange(Params::GetString("digest_mass_range"), doubleRangeParam);
  searchMgr.SetParam("digest_mass_range", Params::GetString("digest_mass_range"), doubleRangeParam);
  
  searchMgr.SetParam("ms_level", Params::GetString("ms_level"), Params::GetInt("ms_level"));
  searchMgr.SetParam("activation_method", Params::GetString("activation_method"), Params::GetString("activation_method"));
  searchMgr.SetParam("minimum_intensity", Params::GetString("minimum_intensity"), Params::GetDouble("minimum_intensity"));
  searchMgr.SetParam("decoy_search", Params::GetString("decoy_search"), Params::GetInt("decoy_search"));
  
  EnzymeInfo enzymeInformation;
  double temp;
  int search_enzyme_number = Params::GetInt("search_enzyme_number");

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

  int sample_enzyme_number = Params::GetInt("sample_enzyme_number");
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
  enzymeInformation.iAllowedMissedCleavage = Params::GetInt("allowed_missed_cleavage");
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
    "input spectra+",
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
vector< pair<string, string> > CometApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("comet.target.txt",
    "a tab-delimited text file containing the target PSMs. See <a href=\"txt-format.html\">"
    "txt file format</a> for a list of the fields."));
  outputs.push_back(make_pair("comet.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  outputs.push_back(make_pair("comet.log.txt",
    "a log file containing a copy of all messages that were printed to "
    "standard error."));
  return outputs;
}

COMMAND_T CometApplication::getCommand() const {
  return COMET_COMMAND;
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
