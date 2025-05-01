/**
 * \file KojakApplication.cpp 
 * \brief Runs Kojak
 *****************************************************************************/
#include "KojakManager.h"
#include "util/AminoAcidUtil.h"
#include "util/CarpStreamBuf.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "KojakApplication.h"
#include "ParamMedicApplication.h"
#include "io/DelimitedFileWriter.h"
#include "io/DelimitedFile.h"
#include "tide/abspath.h"

using namespace std;

/**
 * \returns a blank KojakApplication object
 */
KojakApplication::KojakApplication() {
}

/**
 * Destructor
 */
KojakApplication::~KojakApplication() {
}

/**
 * main method for KojakApplication
 */
int KojakApplication::main(int argc, char** argv) {
  return main(Params::GetStrings("input spectra"));
}

int KojakApplication::main(const vector<string>& input_files) {
  // Vector for parameters that can simply be read with Params::GetString()
  std::vector<std::string> str_params;
  str_params.push_back("threads");
  str_params.push_back("percolator_version");
  str_params.push_back("enrichment");
  str_params.push_back("MS1_resolution");
  str_params.push_back("MS2_resolution");
  str_params.push_back("fixed_modification");
  str_params.push_back("fixed_modification_protC");
  str_params.push_back("fixed_modification_protN");
  str_params.push_back("max_mods_per_peptide");
  str_params.push_back("fragment_bin_size");
  str_params.push_back("decoy_filter");
  str_params.push_back("isotope_error");
  str_params.push_back("max_miscleavages");
  str_params.push_back("max_peptide_mass");
  str_params.push_back("min_peptide_mass");
  str_params.push_back("min_spectrum_peaks");
  str_params.push_back("max_spectrum_peaks");
  str_params.push_back("ppm_tolerance_pre");
  str_params.push_back("prefer_precursor_pred");
  str_params.push_back("top_count");
  str_params.push_back("e_value_depth");
  str_params.push_back("truncate_prot_names");
  str_params.push_back("fragment_bin_offset");
  str_params.push_back("min_peptide_score");
  

  // Vector for boolean parameters, which must be converted to 0/1 strings
  std::vector<std::string> bool_params;
  bool_params.push_back("ion_series_A");
  bool_params.push_back("ion_series_B");
  bool_params.push_back("ion_series_C");
  bool_params.push_back("ion_series_X");
  bool_params.push_back("ion_series_Y");
  bool_params.push_back("ion_series_Z");
  bool_params.push_back("precursor_refinement");
  bool_params.push_back("spectrum_processing");
  bool_params.push_back("export_percolator");
  bool_params.push_back("export_mzID");
  bool_params.push_back("export_pepXML");
  bool_params.push_back("mono_links_on_xl");
  bool_params.push_back("diff_mods_on_xl");
  bool_params.push_back("MS1_centroid");
  bool_params.push_back("MS2_centroid");

  // Vector of list parameters
  // These should be comma-delimited strings of parameters.
  std::vector<std::string> list_params;
  list_params.push_back("modification");
  list_params.push_back("modification_protC");
  list_params.push_back("modification_protN");
  list_params.push_back("cross_link");
  list_params.push_back("mono_link");

  for (vector<string>::const_iterator i = input_files.begin(); i != input_files.end(); i++) {
    if (!StringUtils::IEndsWith(*i, ".mzXML") && !StringUtils::IEndsWith(*i, ".mzML") &&
        !StringUtils::IEndsWith(*i, ".mzXML.gz") && !StringUtils::IEndsWith(*i, ".mzML.gz") &&
        !StringUtils::IEndsWith(*i, ".raw") && !StringUtils::IEndsWith(*i, ".mgf") &&
        !StringUtils::IEndsWith(*i, ".ms2") && !StringUtils::IEndsWith(*i, ".cms2")) {
      carp(CARP_FATAL, "The format of file '%s' is not supported.", i->c_str());
    }
  }

  // Re-route stderr to log file
  CarpStreamBuf buffer;
  streambuf* old_err = cerr.rdbuf();
  cerr.rdbuf(&buffer);

  // set string parameters
  for(size_t i=0; i<str_params.size(); i++) {
    string param = str_params[i] + "=" + Params::GetString(str_params[i]);
    searchManager_.setParam(param);
  }

  // set bool params
  for(size_t i=0; i<bool_params.size(); i++) {
    string param = bool_params[i] + "=" + to_string(Params::GetBool(bool_params[i]));
    searchManager_.setParam(param);
  }

  // set list params
  for(size_t i=0; i<list_params.size(); i++) {
    stringstream param_stream(Params::GetString(list_params[i]));
    while(param_stream.good()) {
      string param_value;
      getline(param_stream, param_value, ',');
      string param = list_params[i] + "=" + param_value;
      searchManager_.setParam(param);
    }
  }

  // Set parameters with different names:
  string instrument = string("instrument=") + Params::GetString("kojak_instrument");
  searchManager_.setParam(instrument);

  stringstream enzyme_stream(Params::GetString("kojak_enzyme"));
  while(enzyme_stream.good()) {
    string enzyme_value;
    getline(enzyme_stream, enzyme_value, ',');
    string enzyme = "enzyme=" + enzyme_value;
    searchManager_.setParam(enzyme);
  }

  string isotope_error = string("isotope_error=") + Params::GetString("kojak_isotope_error");
  searchManager_.setParam(isotope_error);

  string database = string("database=") + Params::GetString("protein database");
  searchManager_.setParam(database);

  // set the output directory
  string output_dir = string("results_path=") + AbsPath(Params::GetString("output-dir"));
  searchManager_.setParam(output_dir);

  // set input files
  int fc=1;
  searchManager_.clearFiles();
  for(size_t i=0;i<input_files.size();i++) fc=searchManager_.setFile(input_files[i].c_str());

  /* Run search */
  bool success = searchManager_.run(); // Returns 0 if successful

  /* Recover stderr */
  cerr.rdbuf(old_err);

  return success ? 1 : 0;
}


void KojakApplication::processParams() {
  const string autoPrecursor = Params::GetString("auto_ppm_tolerance_pre");
  const string autoFragment = Params::GetString("auto_fragment_bin_size");

  if (autoPrecursor != "false" || autoFragment != "false") {
    ParamMedic::RunAttributeResult errorCalcResult;
    vector<ParamMedic::RunAttributeResult> modsResult;
    ParamMedicApplication::processFiles(Params::GetStrings("input spectra"),
      true, false, &errorCalcResult, &modsResult);

    if (autoPrecursor != "false") {
      string fail = errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_FAILURE);
      if (fail.empty()) {
        double sigma = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_SIGMA));
        double prediction = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_PREDICTION));
        carp(CARP_INFO, "precursor ppm standard deviation: %f", sigma);
        carp(CARP_INFO, "Precursor error estimate (ppm): %.2f", prediction);
        Params::Set("ppm_tolerance_pre", prediction);
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
        carp(CARP_INFO, "Fragment bin size estimate (m/z): %.4f", prediction);
        Params::Set("fragment_bin_size", prediction);
      } else {
        carp(autoFragment == "fail" ? CARP_FATAL : CARP_ERROR,
             "failed to calculate fragment error: %s", fail.c_str());
      }
    }
  }
}


/**
 * \returns the command name for KojakApplication
 */
string KojakApplication::getName() const {
  return "kojak";
}

/**
 * \returns the description for KojakApplication
 */
string KojakApplication::getDescription() const {
  return
    "[[nohtml:Identify cross-linked peptide sequences from shotgun mass spectra.]]"
    "[[html:<p>Kojak is a database search tool for the identification of cross-linked "
    "peptides from mass spectra. This search engine was developed and is "
    "maintained by Michael Hoopmann at the Institute for Systems Biology. "
    "Additional information about Kojak can be found at "
    "<a href=\"http://www.kojak-ms.org\">kojak-ms.org</a>.</p>"
    "<p>If you use Kojak, please cite:</p> "
    "<p><a href=\"https://doi.org/10.1021/pr501321h\">\"Kojak: Efficient "
    "analysis of chemically cross-linked protein complexes.\"</a> Hoopmann MR, "
    "Zelter A, Johnson RS, Riffle M, MacCoss MJ, Davis TN, Moritz RL. <em>J "
    "Proteome Res</em> 2015 May 1. doi: 10.1021/pr501321h]]";
}

/**
 * \returns the command arguments
 */
vector<string> KojakApplication::getArgs() const {
  string arr[] = {
    "input spectra+",
    "protein database"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> KojakApplication::getOptions() const {
  string arr[] = {
   // kojak
   "threads",
   "export_percolator",
   "export_mzID",
   "export_pepXML",
   "percolator_version",
   "enrichment",
   "MS1_centroid",
   "MS2_centroid",
   "MS1_resolution",
   "MS2_resolution",
   "cross_link",
   "mono_link",
   "fixed_modification",
   "fixed_modification_protC",
   "fixed_modification_protN",
   "modification",
   "modification_protC",
   "modification_protN",
   "diff_mods_on_xl",
   "max_mods_per_peptide",
   "mono_links_on_xl",
   "fragment_bin_offset",
   "fragment_bin_size",
   "auto_fragment_bin_size",
   "ion_series_A",
   "ion_series_B",
   "ion_series_C",
   "ion_series_X",
   "ion_series_Y",
   "ion_series_Z",
   "decoy_filter",
   "max_miscleavages",
   "max_peptide_mass",
   "min_peptide_mass",
   "min_peptide_score",
   "min_spectrum_peaks",
   "max_spectrum_peaks",
   "ppm_tolerance_pre",
   "auto_ppm_tolerance_pre",
   "precursor_refinement",
   "prefer_precursor_pred",
   "spectrum_processing",
   "top_count",
   "e_value_depth",
   "truncate_prot_names",
   "kojak_instrument",
   "kojak_enzyme",
   "kojak_isotope_error",
   // crux
   "output-dir",
   "overwrite",
   "parameter-file",
   "fileroot",
   "verbosity",
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
   "pm-min-peak-pairs"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > KojakApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("*.kojak.txt",
    "a tab-delimited text file containing the PSMs. See the <a href=\""
    "../file-formats/txt-format.html\">"
    "txt file format</a> for a list of the fields. The \"*\" "
    "indicates the root of the input spectrum file(s)."));
  outputs.push_back(make_pair("*.log", "a log file containing specific "
    "information pertaining to the analysis of the corresponding spectrum "
    "file."));
  outputs.push_back(make_pair("kojak.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  outputs.push_back(make_pair("kojak.log.txt",
    "a log file containing a copy of all messages that were printed to "
    "standard error. These are mainly crux-specific messages."));
  return outputs;
}

COMMAND_T KojakApplication::getCommand() const {
  return KOJAK_COMMAND;
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool KojakApplication::needsOutputDirectory() const {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
