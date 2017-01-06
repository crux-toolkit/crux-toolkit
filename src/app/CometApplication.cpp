/**
 * \file CometApplication.cpp 
 * \brief Runs comet
 *****************************************************************************/
#include "CometSearch/Common.h"
#include "CometSearch/CometSearchManager.h"
#include "util/AminoAcidUtil.h"
#include "util/CarpStreamBuf.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "CometApplication.h"
#include "ParamMedicApplication.h"
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
        !StringUtils::IEndsWith(*i, ".mzXML.gz") && !StringUtils::IEndsWith(*i, ".mzML.gz") &&
        !StringUtils::IEndsWith(*i, ".raw") && !StringUtils::IEndsWith(*i, ".mgf") &&
        !StringUtils::IEndsWith(*i, ".ms2") && !StringUtils::IEndsWith(*i, ".cms2")) {
      carp(CARP_FATAL, "The format of file '%s' is not supported.", i->c_str());
    }
  }

  /* Re-route stderr to log file */
  CarpStreamBuf buffer;
  streambuf* old = std::cerr.rdbuf();
  std::cerr.rdbuf(&buffer);

  /* set parameters */
  vector<InputFileInfo*> pv_input_files;
  setCometParameters(input_files, pv_input_files);
  searchManager_.AddInputFiles(pv_input_files);
  
  /* Run search */
  bool success = searchManager_.DoSearch();

  /* Recover stderr */
  std::cerr.rdbuf(old);

  return success ? 0 : 1;
}

void CometApplication::setString(const string& param) {
  searchManager_.SetParam(param, Params::GetString(param), Params::GetString(param));
}

void CometApplication::setInt(const string& param) {
  searchManager_.SetParam(param, Params::GetString(param), Params::GetInt(param));
}

void CometApplication::setIntRange(const string& param) {
  vector<int> tokens = StringUtils::Fields<int>(Params::GetString(param));
  IntRange r;
  r.iStart = tokens[0];
  r.iEnd = tokens[1];
  searchManager_.SetParam(param, Params::GetString(param), r);
}

void CometApplication::setDouble(const string& param) {
  searchManager_.SetParam(param, Params::GetString(param), Params::GetDouble(param));
}

void CometApplication::setDoubleRange(const string& param) {
  vector<double> tokens = StringUtils::Fields<double>(Params::GetString(param));
  DoubleRange r;
  r.dStart = tokens[0];
  r.dEnd = tokens[1];
  searchManager_.SetParam(param, Params::GetString(param), r);
}

void CometApplication::setDoubleVector(const string& param) {
  vector<double> v = StringUtils::Fields<double>(Params::GetString(param));
  searchManager_.SetParam(param, Params::GetString(param), v);
}

void CometApplication::setVarMod(const string& param) {
  vector<string> fields = StringUtils::Fields(Params::GetString(param));
  if (fields.empty()) {
    return;
  }
  VarMods m;
  for (size_t i = 0; i < fields.size(); i++) {
    string field = fields[i];
    switch (i) {
      case 0: m.dVarModMass = StringUtils::FromString<double>(field); break;
      case 1: strcpy(m.szVarModChar, field.c_str()); break;
      case 2: m.iBinaryMod = StringUtils::FromString<int>(field); break;
      case 3: m.iMaxNumVarModAAPerMod = StringUtils::FromString<int>(field); break;
      case 4: m.iVarModTermDistance = StringUtils::FromString<int>(field); break;
      case 5: m.iWhichTerm = StringUtils::FromString<int>(field); break;
      case 6: m.bRequireThisMod = StringUtils::FromString<int>(field); break;
    }
  }
  searchManager_.SetParam(param, Params::GetString(param), m);
}

void CometApplication::setEnzyme(
  const string& param,
  const string& searchParam,
  const string& sampleParam,
  const string& missedCleavageParam) {
  EnzymeInfo e;
  double temp;
  int search = Params::GetInt(searchParam);
  if (search >= 0 && (size_t)search < get_comet_enzyme_info_lines().size()) {
    const char* szParamBuf = get_comet_enzyme_info_lines()[search].c_str();
    sscanf(szParamBuf, "%lf %48s %d %20s %20s\n",
      &temp,
      e.szSearchEnzymeName,
      &e.iSearchEnzymeOffSet,
      e.szSearchEnzymeBreakAA,
      e.szSearchEnzymeNoBreakAA);
  } else {
    carp(CARP_FATAL, "search_enzyme_number=%d out of range (0-%d)",
      search, get_comet_enzyme_info_lines().size() - 1);
  }

  int sample = Params::GetInt(sampleParam);
  if (sample >= 0 && (size_t)sample < get_comet_enzyme_info_lines().size()) {
    const char* szParamBuf = get_comet_enzyme_info_lines()[sample].c_str();
    sscanf(szParamBuf, "%lf %48s %d %20s %20s\n",
      &temp,
      e.szSampleEnzymeName,
      &e.iSampleEnzymeOffSet,
      e.szSampleEnzymeBreakAA,
      e.szSampleEnzymeNoBreakAA);
  } else {
    carp(CARP_FATAL, "sample_enzyme_number=%d out of range (0-%d)",
      sample, get_comet_enzyme_info_lines().size() - 1);
  }
  e.iAllowedMissedCleavage = Params::GetInt(missedCleavageParam);
  searchManager_.SetParam(param, "TODO", e);
}

/**
 * Sets the parameters for the Comet application using the crux parameters
 */
void CometApplication::setCometParameters(
  const vector<string>& spec_files,
  vector<InputFileInfo*>& pvInputFiles ///<vector of input spectra files
  ) {
  string scan_range_str = Params::GetString("scan_range");
  int analysis_type = AnalysisType_EntireFile;
  int first_scan, last_scan;
  if (scan_range_str != "0 0") {
    analysis_type = AnalysisType_SpecificScanRange;
    vector<string> tokens = StringUtils::Fields(scan_range_str);
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

  // Database
  setString("database_name");
  setInt("decoy_search");
  // CPU threads
  setInt("num_threads");
  // Masses
  setDouble("peptide_mass_tolerance");
  setInt("peptide_mass_units");
  setInt("mass_type_parent");
  setInt("mass_type_fragment");
  setInt("precursor_tolerance_type");
  setInt("isotope_error");
  // Search enzyme
  setInt("search_enzyme_number");
  setInt("num_enzyme_termini");
  setInt("allowed_missed_cleavage");
  // Fragment ions
  setDouble("fragment_bin_tol");
  setDouble("fragment_bin_offset");
  setInt("theoretical_fragment_ions");
  setInt("use_A_ions");
  setInt("use_B_ions");
  setInt("use_C_ions");
  setInt("use_X_ions");
  setInt("use_Y_ions");
  setInt("use_Z_ions");
  setInt("use_NL_ions");
  // Output
  searchManager_.SetParam("output_sqtstream", "0", 0);
  setInt("output_sqtfile");
  setInt("output_txtfile");
  setInt("output_pepxmlfile");
  setInt("output_percolatorfile");
  setInt("output_outfiles");
  setInt("print_expect_score");
  setInt("num_output_lines");
  setInt("show_fragment_ions");
  setInt("sample_enzyme_number");
  // mzXML/mzML parameters
  setIntRange("scan_range");
  setIntRange("precursor_charge");
  setInt("override_charge");
  setInt("ms_level");
  setString("activation_method");
  // Misc. parameters
  setDoubleRange("digest_mass_range");
  setInt("num_results");
  setInt("skip_researching");
  setInt("max_fragment_charge");
  setInt("max_precursor_charge");
  setInt("nucleotide_reading_frame");
  setInt("clip_nterm_methionine");
  setInt("spectrum_batch_size");
  setString("decoy_prefix");
  setString("output_suffix");
  setDoubleVector("mass_offsets");
  // Spectral processing
  setInt("minimum_peaks");
  setDouble("minimum_intensity");
  setInt("remove_precursor_peak");
  setDouble("remove_precursor_tolerance");
  setDoubleRange("clear_mz_range");
  // Variable modifications
  setVarMod("variable_mod01");
  setVarMod("variable_mod02");
  setVarMod("variable_mod03");
  setVarMod("variable_mod04");
  setVarMod("variable_mod05");
  setVarMod("variable_mod06");
  setVarMod("variable_mod07");
  setVarMod("variable_mod08");
  setVarMod("variable_mod09");
  setInt("max_variable_mods_in_peptide");
  setInt("require_variable_mod");
  // Static modifications
  setDouble("add_Cterm_peptide");
  setDouble("add_Nterm_peptide");
  setDouble("add_Cterm_protein");
  setDouble("add_Nterm_protein");
  for (char c = 'A'; c <= 'Z'; c++) {
    string aaName = AminoAcidUtil::GetName(c);
    aaName = aaName.empty() ? "user_amino_acid" : StringUtils::Replace(aaName, " ", "_");
    string param = "add_" + string(1, c) + "_" + aaName;
    setDouble(param);
  }

  setEnzyme("[COMET_ENZYME_INFO]",
            "search_enzyme_number", "sample_enzyme_number", "allowed_missed_cleavage");
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
    "database_name"
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
    // Database
    "decoy_search",
    // CPU threads
    "num_threads",
    // Masses
    "peptide_mass_tolerance",
    "auto_peptide_mass_tolerance",
    "peptide_mass_units",
    "mass_type_parent",
    "mass_type_fragment",
    "precursor_tolerance_type",
    "isotope_error",
    // Search enzyme
    "search_enzyme_number",
    "num_enzyme_termini",
    "allowed_missed_cleavage",
    // Fragment ions
    "fragment_bin_tol",
    "fragment_bin_offset",
    "auto_fragment_bin_tol",
    "theoretical_fragment_ions",
    "use_A_ions",
    "use_B_ions",
    "use_C_ions",
    "use_X_ions",
    "use_Y_ions",
    "use_Z_ions",
    "use_NL_ions",
    // Output
    "output_sqtfile",
    "output_txtfile",
    "output_pepxmlfile",
    "output_percolatorfile",
    "output_outfiles",
    "print_expect_score",
    "num_output_lines",
    "show_fragment_ions",
    "sample_enzyme_number",
    // mzXML/mzML parameters
    "scan_range",
    "precursor_charge",
    "override_charge",
    "ms_level",
    "activation_method",
    // Misc. parameters
    "digest_mass_range",
    "num_results",
    "skip_researching",
    "max_fragment_charge",
    "max_precursor_charge",
    "nucleotide_reading_frame",
    "clip_nterm_methionine",
    "spectrum_batch_size",
    "decoy_prefix",
    "output_suffix",
    "mass_offsets",
    // Spectral processing
    "minimum_peaks",
    "minimum_intensity",
    "remove_precursor_peak",
    "remove_precursor_tolerance",
    "clear_mz_range",
    // Variable modifications
    "variable_mod01",
    "variable_mod02",
    "variable_mod03",
    "variable_mod04",
    "variable_mod05",
    "variable_mod06",
    "variable_mod07",
    "variable_mod08",
    "variable_mod09",
    "max_variable_mods_in_peptide",
    "require_variable_mod",
    // Static modifications
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
    "add_U_selenocysteine",
    "add_V_valine",
    "add_W_tryptophan",
    "add_X_user_amino_acid",
    "add_Y_tyrosine",
    "add_Z_user_amino_acid",
    // param-medic
    "pm-min-precursor-mz",
    "pm-max-precursor-mz",
    "pm-min-frag-mz",
    "pm-max-frag-mz",
    "pm-min-scan-frag-peaks",
    "pm-max-precursor-delta-ppm",
    "pm-charge",
    "pm-top-n-frag-peaks",
    "pm-pair-top-n-frag-peaks",
    "pm-min-common-frag-peaks",
    "pm-max-scan-separation",
    "pm-min-peak-pairs"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > CometApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("comet.target.txt",
    "a tab-delimited text file containing the target PSMs. See <a href=\""
    "../file-formats/txt-format.html\">"
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

void CometApplication::processParams() {
  const string varModPrefix = "variable_mod0";
  const string varModEmptyDefault = "0.0 null 0 4 -1 0 0";
  for (int i = 1; i <= 9; i++) {
    string mod = varModPrefix + StringUtils::ToString(i);
    if (Params::GetString(mod).empty()) {
      Params::Set(mod, varModEmptyDefault);
    }
  }
  // run param-medic?
  const string autoPrecursor = Params::GetString("auto_peptide_mass_tolerance");
  const string autoFragment = Params::GetString("auto_fragment_bin_tol");
  if (autoPrecursor != "false" || autoFragment != "false") {
    if (autoPrecursor != "false" && Params::GetInt("peptide_mass_units") != 2) {
      carp(CARP_FATAL, "Automatic peptide mass tolerance detection is only supported with ppm "
                       "units. Please rerun with either auto_peptide_mass_tolerance set to 'false' "
                       "or peptide_mass_units set to '2'.");
    }
    ParamMedicErrorCalculator errCalc;
    errCalc.processFiles(Params::GetStrings("input spectra"));
    string precursorFailure, fragmentFailure;
    double precursorSigmaPpm = 0;
    double fragmentSigmaPpm = 0;
    double fragmentSigmaTh = 0;
    double precursorPredictionPpm = 0;
    double fragmentPredictionPpm = 0;
    double fragmentPredictionTh = 0;
    errCalc.calcMassErrorDist(&precursorFailure, &fragmentFailure,
                              &precursorSigmaPpm, &fragmentSigmaPpm,
                              &precursorPredictionPpm, &fragmentPredictionTh);

    if (autoPrecursor != "false") {
      if (precursorFailure.empty()) {
        carp(CARP_INFO, "precursor ppm standard deviation: %f", precursorSigmaPpm);
        carp(CARP_INFO, "Precursor error estimate (ppm): %.2f", precursorPredictionPpm);
        Params::Set("peptide_mass_tolerance", precursorPredictionPpm);
      } else {
        carp(autoPrecursor == "fail" ? CARP_FATAL : CARP_ERROR,
             "failed to calculate precursor error: %s", precursorFailure.c_str());
      }
    }
    if (autoFragment != "false") {
      if (fragmentFailure.empty()) {
        carp(CARP_INFO, "fragment standard deviation (ppm): %f", fragmentSigmaPpm);
        carp(CARP_INFO, "Fragment bin size estimate (Th): %.4f", fragmentPredictionTh);
        Params::Set("fragment_bin_tol", fragmentPredictionTh);
      } else {
        carp(autoFragment == "fail" ? CARP_FATAL : CARP_ERROR,
             "failed to calculate fragment error: %s", fragmentFailure.c_str());
      }
    }
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
