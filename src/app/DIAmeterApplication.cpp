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


bool DIAmeterApplication::HAS_DECOYS = true;
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

  bin_width_  = Params::GetDouble("mz-bin-width");
  bin_offset_ = Params::GetDouble("mz-bin-offset");

  bool use_neutral_loss_peaks_ = Params::GetBool("use-neutral-loss-peaks");
  bool use_flanking_peaks_ = Params::GetBool("use-flanking-peaks");

  // Read proteins index file
  ProteinVec proteins;
  pb::Header protein_header;
  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins, proteins_file, &protein_header)) { carp(CARP_FATAL, "Error reading index (%s)", proteins_file.c_str()); }

  // Read auxlocs index file
  vector<const pb::AuxLocation*> locations;
  if (!ReadRecordsToVector<pb::AuxLocation>(&locations, auxlocs_file)) { carp(CARP_FATAL, "Error reading index (%s)", auxlocs_file.c_str()); }

  // Read peptides index file
  pb::Header peptides_header;
  HeadedRecordReader* peptide_reader = new HeadedRecordReader(peptides_file, &peptides_header);
  if ((peptides_header.file_type() != pb::Header::PEPTIDES) || !peptides_header.has_peptides_header()) { carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str()); }

  const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();
  int decoysPerTarget = pepHeader.has_decoys_per_target() ? pepHeader.decoys_per_target() : 0;

  MassConstants::Init(&pepHeader.mods(), &pepHeader.nterm_mods(), &pepHeader.cterm_mods(), &pepHeader.nprotterm_mods(), &pepHeader.cprotterm_mods(), bin_width_, bin_offset_);
  ModificationDefinition::ClearAll();
  TideMatchSet::initModMap(pepHeader.mods(), ANY);
  TideMatchSet::initModMap(pepHeader.nterm_mods(), PEPTIDE_N);
  TideMatchSet::initModMap(pepHeader.cterm_mods(), PEPTIDE_C);
  TideMatchSet::initModMap(pepHeader.nprotterm_mods(), PROTEIN_N);
  TideMatchSet::initModMap(pepHeader.cprotterm_mods(), PROTEIN_C);


  bool overwrite = Params::GetBool("overwrite");
  ofstream* target_file = NULL;
  ofstream* decoy_file = NULL;

  string target_file_name = make_file_path("tide-search.target.txt");
  target_file = create_stream_in_path(target_file_name.c_str(), NULL, overwrite);
  output_file_name_ = target_file_name;
  if (HAS_DECOYS) {
	 string decoy_file_name = make_file_path("tide-search.decoy.txt");
	 decoy_file = create_stream_in_path(decoy_file_name.c_str(), NULL, overwrite);
  }
  TideMatchSet::writeHeaders(target_file, false, decoysPerTarget > 1, true);
  TideMatchSet::writeHeaders(decoy_file, true, decoysPerTarget > 1, true);

  vector<InputFile> ms1_spectra_files = getInputFiles(input_files, 1);
  vector<InputFile> ms2_spectra_files = getInputFiles(input_files, 2);
  // Loop through spectrum files
  for (vector<InputFile>::const_iterator f = ms2_spectra_files.begin(); f != ms2_spectra_files.end(); f++) {
	  string spectra_file = f->SpectrumRecords;
	  SpectrumCollection* spectra = loadMS2Spectra(spectra_file);

	  double highest_ms2_mz = spectra->FindHighestMZ();
	  carp(CARP_DEBUG, "Maximum observed MS2 m/z = %f.", highest_ms2_mz);
	  MaxBin::SetGlobalMax(highest_ms2_mz);

	  // Active queue to process the indexed peptides
	  ActivePeptideQueue* active_peptide_queue = new ActivePeptideQueue(peptide_reader->Reader(), proteins);
	  active_peptide_queue->SetBinSize(bin_width_, bin_offset_);
	  active_peptide_queue->SetOutputs(NULL, &locations, Params::GetInt("top-match"), true, target_file, decoy_file, highest_ms2_mz);

	  delete spectra;
  }


  return 0;
}

vector<InputFile> DIAmeterApplication::getInputFiles(const vector<string>& filepaths, int ms_level) const {
  vector<InputFile> input_sr;

  if (Params::GetString("spectrum-parser") != "pwiz") { carp(CARP_FATAL, "spectrum-parser must be pwiz instead of %s", Params::GetString("spectrum-parser").c_str() ); }

  for (vector<string>::const_iterator f = filepaths.begin(); f != filepaths.end(); f++) {
	  string spectrum_input_url = *f;
	  string spectrumrecords_url = make_file_path(FileUtils::BaseName(spectrum_input_url) + ".spectrumrecords.ms" + to_string(ms_level));
	  carp(CARP_INFO, "Converting %s to spectrumrecords %s", spectrum_input_url.c_str(), spectrumrecords_url.c_str());
      carp(CARP_DEBUG, "New MS%d spectrumrecords filename: %s", ms_level, spectrumrecords_url.c_str());

      if (!FileUtils::Exists(spectrumrecords_url)) {
    	  if (!SpectrumRecordWriter::convert(spectrum_input_url, spectrumrecords_url, ms_level, true)) {
    		  carp(CARP_FATAL, "Error converting MS2 spectrumrecords from %s", spectrumrecords_url.c_str());
    	  }
      }

      /*SpectrumCollection spectra;
	  pb::Header spectrum_header;
      if (!spectra.ReadSpectrumRecords(spectrumrecords_url, &spectrum_header)) {
    	  FileUtils::Remove(spectrumrecords_url);
    	  carp(CARP_FATAL, "Error reading spectrumrecords: %s", spectrumrecords_url.c_str());
    	  continue;
      }*/
      input_sr.push_back(InputFile(*f, spectrumrecords_url, true));
  }

  return input_sr;
}

SpectrumCollection* DIAmeterApplication::loadMS2Spectra(const std::string& file) {
	SpectrumCollection* spectra = new SpectrumCollection();
	pb::Header spectrum_header;

	if (!spectra->ReadSpectrumRecords(file, &spectrum_header)) {
	    carp(CARP_FATAL, "Error reading spectrum file %s", file.c_str());
	}
	if (string_to_window_type(Params::GetString("precursor-window-type")) == WINDOW_MZ) {
		carp(CARP_FATAL, "Precursor-window-type doesn't support mz in DIAmeter!");
	}
	spectra->Sort();

	return spectra;
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
    "file-column",
    "fileroot",
    "max-precursor-charge",
    "min-peaks",
    "mod-precision",
    "mz-bin-offset",
    "mz-bin-width",
    "output-dir",
    "overwrite",
    "fragment-tolerance",
    "precursor-window",
    "precursor-window-type",
    "spectrum-charge",
//    "spectrum-parser",
    "top-match",
    "use-flanking-peaks",
    "use-neutral-loss-peaks",
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

