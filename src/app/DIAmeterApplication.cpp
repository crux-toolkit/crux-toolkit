#include <cstdio>
#include "app/tide/abspath.h"
#include "app/tide/records_to_vector-inl.h"

#include "io/carp.h"
#include "parameter.h"
#include "io/SpectrumRecordWriter.h"
#include "TideIndexApplication.h"
#include "TideSearchApplication.h"
#include "DIAmeterApplication.h"
#include "ParamMedicApplication.h"
#include "PSMConvertApplication.h"
#include "tide/mass_constants.h"
#include "TideMatchSet.h"
#include "util/Params.h"
#include "util/FileUtils.h"
#include "util/StringUtils.h"


bool DIAmeterApplication::HAS_DECOYS = false;
bool DIAmeterApplication::PROTEIN_LEVEL_DECOYS = false;

const double DIAmeterApplication::XCORR_SCALING = 100000000.0;
const double DIAmeterApplication::RESCALE_FACTOR = 20.0;

DIAmeterApplication::DIAmeterApplication() { /* do nothing */ }
DIAmeterApplication::~DIAmeterApplication() { /* do nothing */ }


int DIAmeterApplication::main(int argc, char** argv) {
  return main(Params::GetStrings("tide spectra file"));
}
int DIAmeterApplication::main(const vector<string>& input_files) {
  return main(input_files, Params::GetString("tide database"));
}

int DIAmeterApplication::main(const vector<string>& input_files, const string input_index) {
  carp(CARP_INFO, "Running tide-search...");

  NUM_THREADS = Params::GetInt("num-threads");
  if (NUM_THREADS < 1) {
    NUM_THREADS = boost::thread::hardware_concurrency(); // MINIMUM # = 1.
    // (Meaning just main thread) Do not make this value below 1.
  } else if (NUM_THREADS > 64) {
    // make sure that number of threads are reasonable, e.g. user did not specify millions of threads...
    carp(CARP_FATAL, "Requested more than 64 threads.");
  }
  carp(CARP_INFO, "Number of Threads: %d", NUM_THREADS);

  const string index = input_index;
  string peptides_file = FileUtils::Join(index, "pepix");
  string proteins_file = FileUtils::Join(index, "protix");
  string auxlocs_file = FileUtils::Join(index, "auxlocs");

  // Check spectrum-charge parameter
  string charge_string = Params::GetString("spectrum-charge");
  int charge_to_search;
  if (charge_string == "all") {
    carp(CARP_DEBUG, "Searching all charge states");
    charge_to_search = 0;
  } else {
    charge_to_search = atoi(charge_string.c_str());
    if (charge_to_search < 1 || charge_to_search > 6) {
      carp(CARP_FATAL, "Invalid spectrum-charge value %s", charge_string.c_str());
    }
    carp(CARP_INFO, "Searching charge state %d", charge_to_search);
  }

  // Check scan-number parameter
  string scan_range = Params::GetString("scan-number");
  int min_scan, max_scan;
  if (scan_range.empty()) {
    min_scan = 0;
    max_scan = BILLION;
    carp(CARP_DEBUG, "Searching all scans");
  } else if (scan_range.find('-') == string::npos) {
    // Single scan
    min_scan = max_scan = atoi(scan_range.c_str());
    carp(CARP_INFO, "Searching single scan %d", min_scan);
  } else {
    if (!get_range_from_string(scan_range.c_str(), min_scan, max_scan)) {
      carp(CARP_FATAL, "The scan number range '%s' is invalid. "
           "Must be of the form <first>-<last>.", scan_range.c_str());
    }
    if (min_scan > max_scan) {
      carp(CARP_FATAL, "Invalid scan range: %d to %d.", min_scan, max_scan);
    }
    carp(CARP_INFO, "Searching scan range %d to %d.", min_scan, max_scan);
  }

  bin_width_  = Params::GetDouble("mz-bin-width");
  bin_offset_ = Params::GetDouble("mz-bin-offset");
  bool compute_sp = Params::GetBool("compute-sp");

  ProteinVec proteins;
  carp(CARP_INFO, "Reading index %s", index.c_str());
  // Read proteins index file
  pb::Header protein_header;
  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins, proteins_file, &protein_header)) { carp(CARP_FATAL, "Error reading index (%s)", proteins_file.c_str()); }

  // Read auxlocs index file
  vector<const pb::AuxLocation*> locations;
  if (!ReadRecordsToVector<pb::AuxLocation>(&locations, auxlocs_file)) { carp(CARP_FATAL, "Error reading index (%s)", auxlocs_file.c_str()); }
  carp(CARP_DEBUG, "Read %d auxiliary locations.", locations.size());

  // Read peptides index file
  pb::Header peptides_header;
  vector<HeadedRecordReader*> peptide_reader;
  for (int i = 0; i < NUM_THREADS; i++) { peptide_reader.push_back(new HeadedRecordReader(peptides_file, &peptides_header)); }
  if ((peptides_header.file_type() != pb::Header::PEPTIDES) || !peptides_header.has_peptides_header()) { carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str()); }


  const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();
  DECOY_TYPE_T headerDecoyType = (DECOY_TYPE_T)pepHeader.decoys();
  int decoysPerTarget = pepHeader.has_decoys_per_target() ? pepHeader.decoys_per_target() : 0;
  if (headerDecoyType != NO_DECOYS) {
    HAS_DECOYS = true;
    if (headerDecoyType == PROTEIN_REVERSE_DECOYS) { PROTEIN_LEVEL_DECOYS = true; }
  }

  MassConstants::Init(&pepHeader.mods(), 
      &pepHeader.nterm_mods(), &pepHeader.cterm_mods(),
      &pepHeader.nprotterm_mods(), &pepHeader.cprotterm_mods(),
      bin_width_, bin_offset_);
  ModificationDefinition::ClearAll();
  TideMatchSet::initModMap(pepHeader.mods(), ANY);
  TideMatchSet::initModMap(pepHeader.nterm_mods(), PEPTIDE_N);
  TideMatchSet::initModMap(pepHeader.cterm_mods(), PEPTIDE_C);
  TideMatchSet::initModMap(pepHeader.nprotterm_mods(), PROTEIN_N);
  TideMatchSet::initModMap(pepHeader.cprotterm_mods(), PROTEIN_C);

  ofstream* target_file = NULL;
  ofstream* decoy_file = NULL;

  bool overwrite = Params::GetBool("overwrite");
  stringstream ss;
  ss << Params::GetString("enzyme") << '-' << Params::GetString("digestion");
  TideMatchSet::CleavageType = ss.str();
  if (!Params::GetBool("concat")) {
    string target_file_name = make_file_path("tide-search.target.txt");
    target_file = create_stream_in_path(target_file_name.c_str(), NULL, overwrite);
    output_file_name_ = target_file_name;
    if (HAS_DECOYS) {
      string decoy_file_name = make_file_path("tide-search.decoy.txt");
      decoy_file = create_stream_in_path(decoy_file_name.c_str(), NULL, overwrite);
    }
  } else {
    string concat_file_name = make_file_path("tide-search.txt");
    target_file = create_stream_in_path(concat_file_name.c_str(), NULL, overwrite);
    output_file_name_ = concat_file_name;
  }

  if (target_file) {
    TideMatchSet::writeHeaders(target_file, false, decoysPerTarget > 1, compute_sp);
    TideMatchSet::writeHeaders(decoy_file, true, decoysPerTarget > 1, compute_sp);
  }

  vector<InputFile> sr = getInputFiles(input_files);

  return 0;
}

vector<InputFile> DIAmeterApplication::getInputFiles(
  const vector<string>& filepaths
) const {
  // Try to read all spectrum files as spectrumrecords, convert those that fail
  vector<InputFile> input_sr;
//  for (vector<string>::const_iterator f = filepaths.begin(); f != filepaths.end(); f++) {
//    SpectrumCollection spectra;
//    pb::Header spectrum_header;
//    string spectrumrecords = *f;
//    bool keepSpectrumrecords = true;
//    if (!spectra.ReadSpectrumRecords(spectrumrecords, &spectrum_header)) {
//      // Failed, try converting to spectrumrecords file
//      carp(CARP_INFO, "Converting %s to spectrumrecords format", f->c_str());
//      carp(CARP_INFO, "Elapsed time starting conversion: %.3g s", wall_clock() / 1e6);
//      spectrumrecords = Params::GetString("store-spectra");
//      keepSpectrumrecords = !spectrumrecords.empty();
//      if (!keepSpectrumrecords) {
//        spectrumrecords = make_file_path(FileUtils::BaseName(*f) + ".spectrumrecords.tmp");
//      } else if (filepaths.size() > 1) {
//        carp(CARP_FATAL, "Cannot use store-spectra option with multiple input "
//                         "spectrum files");
//      }
//      carp(CARP_DEBUG, "New spectrumrecords filename: %s", spectrumrecords.c_str());
//      if (!SpectrumRecordWriter::convert(*f, spectrumrecords)) {
//        carp(CARP_FATAL, "Error converting %s to spectrumrecords format", f->c_str());
//      }
//      carp(CARP_DEBUG, "Reading converted spectrum file %s", spectrumrecords.c_str());
//      // Re-read converted file as spectrumrecords file
//      if (!spectra.ReadSpectrumRecords(spectrumrecords, &spectrum_header)) {
//        carp(CARP_DEBUG, "Deleting %s", spectrumrecords.c_str());
//        FileUtils::Remove(spectrumrecords);
//        carp(CARP_FATAL, "Error reading spectra file %s", spectrumrecords.c_str());
//      }
//    }
//    input_sr.push_back(InputFile(*f, spectrumrecords, keepSpectrumrecords));
//  }
  return input_sr;
}



string DIAmeterApplication::getName() const {
  return "diameter";
}

string DIAmeterApplication::getDescription() const {
  return "DIAmeter description";
}

vector<string> DIAmeterApplication::getArgs() const {
  string arr[] = {
    "tide spectra file+",
    "tide database"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> DIAmeterApplication::getOptions() const {
  string arr[] = {
    "auto-mz-bin-width",
    "auto-precursor-window",
    "compute-sp",
    "concat",
    "deisotope",
    "elution-window-size",
    "file-column",
    "fileroot",
    "mass-precision",
    "max-precursor-charge",
    "min-peaks",
    "mod-precision",
    "mz-bin-offset",
    "mz-bin-width",
    "mzid-output",
    "num-threads",
    "output-dir",
    "overwrite",
    "parameter-file",
    "fragment-tolerance",
    "evidence-granularity",
    "pepxml-output",
    "pin-output",
    "pm-charges",
    "pm-max-frag-mz",
    "pm-max-precursor-delta-ppm",
    "pm-max-precursor-mz",
    "pm-max-scan-separation",
    "pm-min-common-frag-peaks",
    "pm-min-frag-mz",
    "pm-min-peak-pairs",
    "pm-min-precursor-mz",
    "pm-min-scan-frag-peaks",
    "pm-pair-top-n-frag-peaks",
    "pm-top-n-frag-peaks",
    "precision",
    "precursor-window",
    "precursor-window-type",
    "print-search-progress",
    "remove-precursor-peak",
    "remove-precursor-tolerance",
    "scan-number",
    "skip-preprocessing",
    "spectrum-charge",
    "spectrum-max-mz",
    "spectrum-min-mz",
    "spectrum-parser",
    "store-index",
    "store-spectra",
    "top-match",
    "txt-output",
    "brief-output",
    "use-flanking-peaks",
    "use-neutral-loss-peaks",
    "use-z-line",
    "use-tailor-calibration",    
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector< pair<string, string> > DIAmeterApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("tide-search.target.txt",
    "a tab-delimited text file containing the target PSMs. See <a href=\""
    "../file-formats/txt-format.html\">txt file format</a> for a list of the fields."));
  outputs.push_back(make_pair("tide-search.decoy.txt",
    "a tab-delimited text file containing the decoy PSMs. This file will only "
    "be created if the index was created with decoys."));
  outputs.push_back(make_pair("tide-search.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other Crux programs."));
  outputs.push_back(make_pair("tide-search.log.txt",
    "a log file containing a copy of all messages that were printed to the "
    "screen during execution."));
  return outputs;
}
bool DIAmeterApplication::needsOutputDirectory() const {
  return true;
}

COMMAND_T DIAmeterApplication::getCommand() const {
  return DIAMETER_COMMAND;
}

void DIAmeterApplication::processParams() {

}

