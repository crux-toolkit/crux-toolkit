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

using namespace std;

/**
 * \returns a blank KojakApplication object
 */
KojakApplication::KojakApplication() {
  //populate all the parameters here. Update for parameters added and removed.
  kojak_params.push_back("threads");
  kojak_params.push_back("database");
  kojak_params.push_back("export_percolator");
  kojak_params.push_back("export_pepXML");
  kojak_params.push_back("MS_data_file");
  kojak_params.push_back("percolator_version");
  kojak_params.push_back("enrichment");
  //kojak_params.push_back("kojak_instrument");
  kojak_params.push_back("MS1_centroid");
  kojak_params.push_back("MS2_centroid");
  kojak_params.push_back("MS1_resolution");
  kojak_params.push_back("MS2_resolution");
  kojak_params.push_back("cross_link");
  kojak_params.push_back("mono_link");
  kojak_params.push_back("fixed_modification");
  kojak_params.push_back("fixed_modification_protC");
  kojak_params.push_back("fixed_modification_protN");
  kojak_params.push_back("modification");
  kojak_params.push_back("modification_protC");
  kojak_params.push_back("modification_protN");
  kojak_params.push_back("diff_mods_on_xl");
  kojak_params.push_back("max_mods_per_peptide");
  kojak_params.push_back("mono_links_on_xl");
  kojak_params.push_back("fragment_bin_offset");
  kojak_params.push_back("fragment_bin_size");
  kojak_params.push_back("ion_series_A");
  kojak_params.push_back("ion_series_B");
  kojak_params.push_back("ion_series_C");
  kojak_params.push_back("ion_series_X");
  kojak_params.push_back("ion_series_Y");
  kojak_params.push_back("ion_series_Z");
  kojak_params.push_back("decoy_filter");
  kojak_params.push_back("isotope_error");
  kojak_params.push_back("max_miscleavages");
  kojak_params.push_back("max_peptide_mass");
  kojak_params.push_back("min_peptide_mass");
  kojak_params.push_back("max_spectrum_peaks");
  kojak_params.push_back("ppm_tolerance_pre");
  kojak_params.push_back("prefer_precursor_pred");
  kojak_params.push_back("spectrum_processing");
  kojak_params.push_back("top_count");
  kojak_params.push_back("truncate_prot_names");
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
  streambuf* old = std::cerr.rdbuf();
  std::cerr.rdbuf(&buffer);

  // set parameters
  for(size_t i=0;i<kojak_params.size();i++){
    string param = kojak_params[i] + "=" + Params::GetString(kojak_params[i]);
    searchManager_.setParam(param);
  }

  string instrument = string("instrument") + "=" + Params::GetString("kojak_instrument");
  string enzyme = string("enzyme") + "=" + Params::GetString("kojak_enzyme");
  string fragment_bin_offset = string("fragment_bin_offset") + "=" + Params::GetString("kojak_fragment_bin_offset");
  string isotope_error = string("isotope_error") + "=" + Params::GetString("kojak_isotope_error");
  searchManager_.setParam(instrument);
  searchManager_.setParam(enzyme);
  searchManager_.setParam(fragment_bin_offset);
  searchManager_.setParam(isotope_error);


  // set input files
  int fc=1;
  searchManager_.clearFiles();
  for(size_t i=0;i<input_files.size();i++) fc=searchManager_.setFile(input_files[i].c_str());

  /* Run search */
  bool success = searchManager_.run();

  /* Recover stderr */
  std::cerr.rdbuf(old);

  return success ? 0 : 1;
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
    "[[html:<p>Kojak.</p>]]";
}

/**
 * \returns the command arguments
 */
vector<string> KojakApplication::getArgs() const {
  string arr[] = {
    "input spectra+"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> KojakApplication::getOptions() const {
  return kojak_params;
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > KojakApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("kojak.target.txt",
    "a tab-delimited text file containing the target PSMs. See <a href=\""
    "../file-formats/txt-format.html\">"
    "txt file format</a> for a list of the fields."));
  outputs.push_back(make_pair("kojak.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  outputs.push_back(make_pair("kojak.log.txt",
    "a log file containing a copy of all messages that were printed to "
    "standard error."));
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
