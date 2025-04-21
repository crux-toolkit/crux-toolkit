/**
 * \file CometApplication.cpp 
 * \brief Runs comet
 *****************************************************************************/
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
      case 6: m.iRequireThisMod = StringUtils::FromString<int>(field); break;
    }
  }
  searchManager_.SetParam(param, Params::GetString(param), m);
}

void CometApplication::setEnzyme(
  const string& param,
  const string& searchParam,
  const string& search2Param,
  const string& sampleParam,
  const string& missedCleavageParam) {
  EnzymeInfo e;
  double temp;

  // Search enzyme 1
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
    size_t numLines = get_comet_enzyme_info_lines().size() - 1;
    carp(CARP_FATAL, "%s=%d out of range (0-%d)", searchParam.c_str(), search, numLines);
  }

  // Search enzyme 2
  int search2 = Params::GetInt(search2Param);
  if (search2 >= 0 && (size_t)search2 < get_comet_enzyme_info_lines().size()) {
    const char* szParamBuf = get_comet_enzyme_info_lines()[search2].c_str();
    sscanf(szParamBuf, "%lf %48s %d %20s %20s\n",
      &temp,
      e.szSearchEnzyme2Name,
      &e.iSearchEnzyme2OffSet,
      e.szSearchEnzyme2BreakAA,
      e.szSearchEnzyme2NoBreakAA);
  } else {
    size_t numLines = get_comet_enzyme_info_lines().size() - 1;
    carp(CARP_FATAL, "%s=%d out of range (0-%d)", search2Param.c_str(), search2, numLines);
  }

  // Sample enzyme
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
  setInt("peff_format");
  setString("peff_obo");
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
  setInt("search_enzyme2_number");
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
  setInt("use_Z1_ions");
  setInt("use_NL_ions");
  // Output
  setInt("output_mzidentmlfile");
  setInt("output_pepxmlfile");
  setInt("output_percolatorfile");
  setInt("output_sqtfile");
  searchManager_.SetParam("output_sqtstream", "0", 0);
  setInt("output_sqtstream");
  setInt("output_txtfile");
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
  setInt("clip_nterm_methionine");
  setString("decoy_prefix");
  setDoubleRange("digest_mass_range");
  setInt("equal_I_and_L");
  setDoubleVector("mass_offsets");
  setInt("max_duplicate_proteins");
  setInt("max_fragment_charge");
  setInt("max_index_runtime");
  setInt("max_precursor_charge");
  setInt("num_results");
  setInt("nucleotide_reading_frame");
  setString("output_suffix");
  setInt("peff_verbose_output");
  setIntRange("peptide_length_range");
  setDoubleVector("precursor_NL_ions");
  setInt("skip_researching");
  setInt("spectrum_batch_size");
  setString("text_file_extension");
  setInt("explicit_deltacn");
  setInt( "old_mods_encoding");
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
  setString("require_variable_mod");
  setString("protein_modlist_file");
  // Static modifications
  setDouble("add_Cterm_peptide");
  setDouble("add_Nterm_peptide");
  setDouble("add_Cterm_protein");
  setDouble("add_Nterm_protein");
  for (char c = 'A'; c <= 'Z'; c++) {
    setDouble(staticModParam(c));
  }

  setEnzyme("[COMET_ENZYME_INFO]",
            "search_enzyme_number", "search_enzyme2_number", "sample_enzyme_number",
            "allowed_missed_cleavage");

  // Indexing
  searchManager_.SetParam("create_fragment_index", "0", 0);
  searchManager_.SetParam("create_peptide_index", "0", 0);
  setDouble("fragindex_max_fragmentmass");
  setDouble("fragindex_min_fragmentmass");
  setInt("fragindex_min_ions_report");
  setInt("fragindex_min_ions_score");
  setInt("fragindex_num_spectrumpeaks");
  setInt("fragindex_skipreadprecursors");
}

string CometApplication::staticModParam(char c) {
  if (c < 'A' || c > 'Z') {
    return "";
  }
  string aaName = AminoAcidUtil::GetName(c);
  aaName = aaName.empty() ? "user_amino_acid" : StringUtils::Replace(aaName, " ", "_");
  return "add_" + string(1, c) + "_" + aaName;
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
    "protein database in FASTA format or an index created from a FASTA file.]]"
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
    "peff_format",
    "peff_obo",
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
    "search_enzyme2_number",
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
    "use_Z1_ions",
    "use_NL_ions",
    // Output
    "output_mzidentmlfile",
    "output_pepxmlfile",
    "output_percolatorfile",
    "output_sqtfile",
    "output_sqtstream",
    "output_txtfile",
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
    "clip_nterm_methionine",
    "decoy_prefix",
    "digest_mass_range",
    "equal_I_and_L",
    "mass_offsets",
    "max_duplicate_proteins",
    "max_fragment_charge",
    "max_index_runtime",
    "max_precursor_charge",
    "num_results",
    "nucleotide_reading_frame",
    "output_suffix",
    "peff_verbose_output",
    "peptide_length_range",
    "precursor_NL_ions",
    "skip_researching",
    "spectrum_batch_size",
    "text_file_extension",
    "explicit_deltacn",
    "old_mods_encoding",
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
    "auto_modifications",
    "max_variable_mods_in_peptide",
    "require_variable_mod",
    "protein_modlist_file",
    // Static modifications
    "add_Cterm_peptide",
    "add_Nterm_peptide",
    "add_Cterm_protein",
    "add_Nterm_protein",
    "add_A_alanine",
    "add_C_cysteine",
    "add_D_aspartic_acid",
    "add_E_glutamic_acid",
    "add_F_phenylalanine",
    "add_G_glycine",
    "add_H_histidine",
    "add_I_isoleucine",
    "add_K_lysine",
    "add_L_leucine",
    "add_M_methionine",
    "add_N_asparagine",
    "add_O_pyrrolysine",
    "add_P_proline",
    "add_Q_glutamine",
    "add_R_arginine",
    "add_S_serine",
    "add_T_threonine",
    "add_U_selenocysteine",
    "add_V_valine",
    "add_W_tryptophan",
    "add_Y_tyrosine",
    "add_B_user_amino_acid",
    "add_J_user_amino_acid",
    "add_X_user_amino_acid",
    "add_Z_user_amino_acid",
    // param-medic
    "pm-min-precursor-mz",
    "pm-max-precursor-mz",
    "pm-min-frag-mz",
    "pm-max-frag-mz",
    "pm-min-scan-frag-peaks",
    "pm-max-precursor-delta-ppm",
    "pm-charges",
    "pm-top-n-frag-peaks",
    "pm-pair-top-n-frag-peaks",
    "pm-min-common-frag-peaks",
    "pm-max-scan-separation",
    "pm-min-peak-pairs",
  // Fragment ion fragment indexing
    "fragindex_skipreadprecursors",
    "fragindex_num_spectrumpeaks",
    "fragindex_min_ions_score",
    "fragindex_min_ions_report",
    "fragindex_min_fragmentmass",
    "fragindex_max_fragmentmass"
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
  const bool autoMods = Params::GetBool("auto_modifications");
  if (autoPrecursor != "false" || autoFragment != "false" || autoMods) {
    if (autoPrecursor != "false" && Params::GetInt("peptide_mass_units") != 2) {
      carp(CARP_FATAL, "Automatic peptide mass tolerance detection is only supported with ppm "
                       "units. Please rerun with either auto_peptide_mass_tolerance set to 'false' "
                       "or peptide_mass_units set to '2'.");
    } else if (autoMods) {
      // user cannot specify own mods if param-medic is going to infer them
      for (int i = 1; i <= 9; i++) {
        string mod = varModPrefix + StringUtils::ToString(i);
        if (!Params::IsDefault(mod)) {
          carp(CARP_FATAL, "Automatic modification inference cannot be used with user-specified "
                           "modifications. Please rerun with either auto_modifications set to 'false' "
                           "or with modifications turned off.");
        }
      }
      for (char c = 'A'; c <= 'Z'; c++) {
        string param = staticModParam(c);
        if (!Params::IsDefault(param)) {
          carp(CARP_FATAL, "Automatic modification inference cannot be used with user-specified "
                           "modifications. Please rerun with either auto_modifications set to 'false' "
                           "or with modifications turned off.");
        }
      }
    }
    ParamMedic::RunAttributeResult errorCalcResult;
    vector<ParamMedic::RunAttributeResult> modsResult;
    ParamMedicApplication::processFiles(Params::GetStrings("input spectra"),
      autoPrecursor != "false" || autoFragment != "false", autoMods, &errorCalcResult, &modsResult);

    if (autoPrecursor != "false") {
      string fail = errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_FAILURE);
      if (fail.empty()) {
        double sigma = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_SIGMA));
        double prediction = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_PREDICTION));
        carp(CARP_INFO, "precursor ppm standard deviation: %f", sigma);
        carp(CARP_INFO, "Precursor error estimate (ppm): %.2f", prediction);
        Params::Set("peptide_mass_tolerance", prediction);
      } else {
        carp(autoPrecursor == "fail" ? CARP_FATAL : CARP_ERROR,
             "failed to calculate precursor error: %s", fail.c_str());
      }
    }
    if (autoFragment != "false") {
      string fail = errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_FRAGMENT_FAILURE);
      if (fail.empty()) {
        double sigma = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_FRAGMENT_SIGMA));
        double prediction = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_FRAGMENT_PREDICTION));
        carp(CARP_INFO, "fragment ppm standard deviation: %f", sigma);
        carp(CARP_INFO, "Fragment bin size estimate (Th): %.4f", prediction);
        Params::Set("fragment_bin_tol", prediction);
      } else {
        carp(autoFragment == "fail" ? CARP_FATAL : CARP_ERROR,
             "failed to calculate fragment error: %s", fail.c_str());
      }
    }
    if (autoMods) {
      int modNum = 1;
      vector<ParamMedic::Modification> mods = ParamMedic::Modification::GetFromResults(modsResult);
      for (vector<ParamMedic::Modification>::const_iterator i = mods.begin(); i != mods.end(); i++) {
        string location = i->getLocation();
        const double mass = i->getMassDiff();
        const bool variable = i->getVariable();
        if (variable) {
          string distance = "-1";
          string terminus = "0";
          if (location == ParamMedic::Modification::LOCATION_NTERM) {
            location = "n";
            distance = "0";
            terminus = "2";
          } else if (location == ParamMedic::Modification::LOCATION_CTERM) {
            location = "c";
            distance = "0";
            terminus = "3";
          }
          Params::Set("variable_mod0" + StringUtils::ToString(modNum++),
            StringUtils::ToString(mass) + " " + location + " 0 4 " + distance + " " + terminus + " 0");
        } else {
          vector<string> params;
          if (location == ParamMedic::Modification::LOCATION_NTERM) {
            params.push_back("add_Nterm_peptide");
          } else if (location == ParamMedic::Modification::LOCATION_CTERM) {
            params.push_back("add_Cterm_peptide");
          } else {
            for (string::const_iterator i = location.begin(); i != location.end(); i++) {
              params.push_back(staticModParam(*i));
            }
          }
          for (vector<string>::const_iterator i = params.begin(); i != params.end(); i++) {
            double existing = Params::GetDouble(*i);
            Params::Set(*i, existing + mass);
          }
        }
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
