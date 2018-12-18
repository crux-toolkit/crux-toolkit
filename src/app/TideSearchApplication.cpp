#include <cstdio>
#include "app/tide/abspath.h"
#include "app/tide/records_to_vector-inl.h"

#include "io/carp.h"
#include "parameter.h"
#include "io/SpectrumRecordWriter.h"
#include "TideIndexApplication.h"
#include "TideSearchApplication.h"
#include "ParamMedicApplication.h"
#include "PSMConvertApplication.h"
#include "tide/mass_constants.h"
#include "TideMatchSet.h"
#include "util/Params.h"
#include "util/FileUtils.h"
#include "util/StringUtils.h"
#include <math.h> //Added by Andy Lin
#include <map> //Added by Andy Lin

bool TideSearchApplication::HAS_DECOYS = false;
bool TideSearchApplication::PROTEIN_LEVEL_DECOYS = false;

/* This constant is the product of the original "magic number" (10000,
 * on line 4622 of search28.c) that was used to rescale the XCorr
 * score, and the integerization constant used by Benjamin Diament in
 * Tide.  In the Tide publication, that constant was reported as 10^7.
 * However, here it appears to be only 10^4.
 *
 * --WSN, 10 March 2015 */
const double TideSearchApplication::XCORR_SCALING = 100000000.0;

/* This constant is used to put the refactored XCorr back into the
 * same range as the original XCorr score.  It is the XCorr "magic
 * number" (10000) divided by the EVIDENCE_SCALE_INT (defined in
 * tide/spectrum_preprocess2.cc). */
const double TideSearchApplication::RESCALE_FACTOR = 20.0;

TideSearchApplication::TideSearchApplication():
  exact_pval_search_(false), remove_index_(""), spectrum_flag_(NULL) {
}

TideSearchApplication::~TideSearchApplication() {
  if (!remove_index_.empty()) {
    carp(CARP_DEBUG, "Removing temp index '%s'", remove_index_.c_str());
    FileUtils::Remove(remove_index_);
  }
}

int TideSearchApplication::main(int argc, char** argv) {
  return main(Params::GetStrings("tide spectra file"));
}
int TideSearchApplication::main(const vector<string>& input_files) {
  return main(input_files, Params::GetString("tide database"));
}

int TideSearchApplication::main(const vector<string>& input_files, const string input_index) {
  carp(CARP_INFO, "Running tide-search...");

  // prevent different output formats from using threading
  if (!Params::GetBool("peptide-centric-search")) {
    NUM_THREADS = Params::GetInt("num-threads");
  } else {
    carp(CARP_INFO, "Threading for peptide-centric formats is not yet supported.");
    NUM_THREADS = 1;
  }
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

  // check to compute exact p-value
  exact_pval_search_ = Params::GetBool("exact-p-value");
  bin_width_  = Params::GetDouble("mz-bin-width");
  bin_offset_ = Params::GetDouble("mz-bin-offset");
  // for now don't allow XCorr p-value or combined p-values
  // searches with variable bin width
  if (exact_pval_search_ && !Params::IsDefault("mz-bin-width")) {
    carp(CARP_FATAL, "Tide-search with variable bin width and "
                     "XCorr p-value or combined p-value "
                     "is not allowed in this version of Crux.");
  }

  // Check that combined p-value is with exact-p-value=T
  SCORE_FUNCTION_T curScoreFunction = string_to_score_function_type(Params::GetString("score-function"));
  if (curScoreFunction == BOTH_SCORE && !exact_pval_search_) {
    carp(CARP_FATAL,"--score-function 'both' must be run with exact-p-value=T");
  }

  // TODO res-ev not implemented with flanking_peak, neutral_loss_peak,
  // and dynamic CTerm and NTerm mass
  bool use_neutral_loss_peaks_ = Params::GetBool("use-neutral-loss-peaks");
  bool use_flanking_peaks_ = Params::GetBool("use-flanking-peaks");
  if (use_flanking_peaks_ && curScoreFunction == RESIDUE_EVIDENCE_MATRIX) {
    carp(CARP_FATAL,"--score-function 'residue-evidence' with "
                    "--use-flanking-peaks true is not implemented yet");
  }
  if (use_neutral_loss_peaks_ && curScoreFunction == RESIDUE_EVIDENCE_MATRIX) {
    carp(CARP_FATAL,"--score-function 'residue-evidence' with "
                    "--use-neutral-loss-peaks 'true' is not implemented yet");
  }

  // Check compute-sp parameter
  bool compute_sp = Params::GetBool("compute-sp");
  if (Params::GetBool("sqt-output") && !compute_sp) {
    compute_sp = true;
    carp(CARP_INFO, "Setting compute-sp=T because SQT output is enabled.");
  }

  vector<int> negative_isotope_errors = getNegativeIsotopeErrors();

  ProteinVec proteins;
  carp(CARP_INFO, "Reading index %s", index.c_str());
  // Read proteins index file
  pb::Header protein_header;
  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins,
      proteins_file, &protein_header)) {
    carp(CARP_FATAL, "Error reading index (%s)", proteins_file.c_str());
  }
  int64_t targetProteinCount = 0;
  for (ProteinVec::const_iterator i = proteins.begin(); i != proteins.end(); i++) {
    if (!(*i)->has_target_pos()) {
      ++targetProteinCount;
    }
  }
  carp(CARP_INFO, "Read %d target proteins", targetProteinCount);

  // Open a copy of peptide buffer for Amino Acid Frequency (AAF) calculation.
  // Is used for curScoreFunction = XCORR_SCORE
  // contains frequency of amino acid masses (masses are binned)
  double* aaFreqN = NULL;
  double* aaFreqI = NULL;
  double* aaFreqC = NULL;
  int* aaMass = NULL;
  int nAA = 0;

  // Is used for curScoreFunction = RESIDUE_EVIDENCE_MATRIX
  // contains frequency of amino acid masses (masses are not binned)
  vector<double> dAAFreqN;
  vector<double> dAAFreqI;
  vector<double> dAAFreqC;
  vector<double> dAAMass;
  int nAARes = 0;

  if (curScoreFunction == RESIDUE_EVIDENCE_MATRIX || curScoreFunction == BOTH_SCORE) {
    pb::Header aaf_peptides_header;
    HeadedRecordReader aaf_peptide_reader(peptides_file, &aaf_peptides_header);

    if (aaf_peptides_header.file_type() != pb::Header::PEPTIDES ||
        !aaf_peptides_header.has_peptides_header()) {
      carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
    }

    MassConstants::Init(&aaf_peptides_header.peptides_header().mods(),
                        &aaf_peptides_header.peptides_header().nterm_mods(),
                        &aaf_peptides_header.peptides_header().cterm_mods(),
                        bin_width_, bin_offset_);

    ActivePeptideQueue* active_peptide_queue =
      new ActivePeptideQueue(aaf_peptide_reader.Reader(), proteins);

    nAARes = active_peptide_queue->CountAAFrequencyRes(bin_width_, bin_offset_,
                                                       dAAFreqN, dAAFreqI, dAAFreqC, dAAMass);
    delete active_peptide_queue;
  }

  // For SCORE_FUNCTION=="XCORR_SCORE" with p-val=T or SCORE_FUNCTION=="BOTH_SCORE"
  if (exact_pval_search_ && (curScoreFunction == XCORR_SCORE || curScoreFunction == BOTH_SCORE)) {
    pb::Header aaf_peptides_header;
    HeadedRecordReader aaf_peptide_reader(peptides_file, &aaf_peptides_header);

    if (( aaf_peptides_header.file_type() != pb::Header::PEPTIDES) ||
         !aaf_peptides_header.has_peptides_header()) {
      carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
    }

    MassConstants::Init(&aaf_peptides_header.peptides_header().mods(),
                        &aaf_peptides_header.peptides_header().nterm_mods(),
                        &aaf_peptides_header.peptides_header().cterm_mods(),
                        bin_width_, bin_offset_);

    ActivePeptideQueue* active_peptide_queue =
      new ActivePeptideQueue(aaf_peptide_reader.Reader(), proteins);

    nAA = active_peptide_queue->CountAAFrequency(bin_width_, bin_offset_,
                                                 &aaFreqN, &aaFreqI, &aaFreqC, &aaMass);
    delete active_peptide_queue;
  } // End calculation of amino acid frequencies.

  // Read auxlocs index file
  vector<const pb::AuxLocation*> locations;
  if (!ReadRecordsToVector<pb::AuxLocation>(&locations, auxlocs_file)) {
    carp(CARP_FATAL, "Error reading index (%s)", auxlocs_file.c_str());
  }
  carp(CARP_DEBUG, "Read %d auxiliary locations.", locations.size());

  // Read peptides index file
  pb::Header peptides_header;

  vector<HeadedRecordReader*> peptide_reader;
  for (int i = 0; i < NUM_THREADS; i++) {
    peptide_reader.push_back(new HeadedRecordReader(peptides_file, &peptides_header));
  }

  if ((peptides_header.file_type() != pb::Header::PEPTIDES) ||
      !peptides_header.has_peptides_header()) {
    carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
  }

  const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();
  DECOY_TYPE_T headerDecoyType = (DECOY_TYPE_T)pepHeader.decoys();
  int decoysPerTarget = pepHeader.has_decoys_per_target() ? pepHeader.decoys_per_target() : 0;
  if (headerDecoyType != NO_DECOYS) {
    HAS_DECOYS = true;
    if (headerDecoyType == PROTEIN_REVERSE_DECOYS) {
      PROTEIN_LEVEL_DECOYS = true;
    }
  }

  MassConstants::Init(&pepHeader.mods(), &pepHeader.nterm_mods(),
    &pepHeader.cterm_mods(), bin_width_, bin_offset_);
  ModificationDefinition::ClearAll();
  TideMatchSet::initModMap(pepHeader.mods(), ANY);
  TideMatchSet::initModMap(pepHeader.nterm_mods(), PEPTIDE_N);
  TideMatchSet::initModMap(pepHeader.cterm_mods(), PEPTIDE_C);

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

  // Loop through spectrum files
  for (vector<InputFile>::const_iterator f = sr.begin(); f != sr.end(); f++) {
    if (!peptide_reader[0]) {
      for (int i = 0; i < NUM_THREADS; i++) {
        peptide_reader[i] = new HeadedRecordReader(peptides_file, &peptides_header);
      }
    }

    vector<ActivePeptideQueue*> active_peptide_queue;
    for (int i = 0; i < NUM_THREADS; i++) {
      active_peptide_queue.push_back(new ActivePeptideQueue(peptide_reader[i]->Reader(), proteins));
      active_peptide_queue[i]->SetBinSize(bin_width_, bin_offset_);
    }

    string spectra_file = f->SpectrumRecords;
    SpectrumCollection* spectra = NULL;
    map<string, SpectrumCollection*>::iterator spectraIter = spectra_.find(spectra_file);
    if (spectraIter == spectra_.end()) {
      carp(CARP_INFO, "Reading spectrum file %s.", spectra_file.c_str());
      spectra = loadSpectra(spectra_file);
      carp(CARP_INFO, "Read %d spectra.", spectra->Size());
    } else {
      spectra = spectraIter->second;
    }

    double highest_mz = spectra->FindHighestMZ();
    unsigned int spectrum_num = spectra->SpecCharges()->size();
    if (spectrum_num > 0 &&
        (exact_pval_search_ || curScoreFunction == RESIDUE_EVIDENCE_MATRIX || curScoreFunction == BOTH_SCORE)) {
      highest_mz = spectra->SpecCharges()->at(spectrum_num - 1).neutral_mass;
    }
    carp(CARP_DEBUG, "Maximum observed m/z = %f.", highest_mz);
    MaxBin::SetGlobalMax(highest_mz);
    // Do the search
    carp(CARP_INFO, "Starting search.");
    if (spectrum_flag_ == NULL) {
      resetMods();
    }
    search(f->OriginalName, spectra->SpecCharges(), active_peptide_queue, proteins,
           locations, Params::GetDouble("precursor-window"),
           string_to_window_type(Params::GetString("precursor-window-type")),
           Params::GetDouble("spectrum-min-mz"), Params::GetDouble("spectrum-max-mz"),
           min_scan, max_scan, Params::GetInt("min-peaks"), charge_to_search,
           Params::GetInt("top-match"), spectra->FindHighestMZ(),
           target_file, decoy_file, compute_sp,
           nAA, aaFreqN, aaFreqI, aaFreqC, aaMass,
           nAARes, dAAFreqN, dAAFreqI, dAAFreqC, dAAMass,
           pepHeader.mods(), pepHeader.nterm_mods(), pepHeader.cterm_mods(),
           decoysPerTarget, &negative_isotope_errors);

    if (spectraIter == spectra_.end()) {
      delete spectra;
    }
    // convert tab delimited to other file formats.
    convertResults();

    // Delete temporary spectrumrecords file
    if (!f->Keep) {
      carp(CARP_DEBUG, "Deleting %s", spectra_file.c_str());
      remove(spectra_file.c_str());
    }

    // Clean up
    for (int i = 0; i < NUM_THREADS; i++) {
      delete active_peptide_queue[i];
      delete peptide_reader[i];
      peptide_reader[i] = NULL;
    }

  } // End of spectrum file loop

  for (ProteinVec::iterator i = proteins.begin(); i != proteins.end(); ++i) {
    delete *i;
  }
  if (target_file) {
    delete target_file;
    if (decoy_file) {
      delete decoy_file;
    }
  }
  delete[] aaFreqN;
  delete[] aaFreqI;
  delete[] aaFreqC;
  delete[] aaMass;

  return 0;
}

vector<int> TideSearchApplication::getNegativeIsotopeErrors() const {
  string isotope_errors_string = Params::GetString("isotope-error");
  if (isotope_errors_string[0] == ',') {
    carp(CARP_FATAL, "Error in isotope_error parameter formatting: (%s)",
         isotope_errors_string.c_str());
  }
  for (string::const_iterator i = isotope_errors_string.begin(); i != isotope_errors_string.end(); i++) {
    if (*i == ',' && (i+1 == isotope_errors_string.end() || *(i+1) == ',')) {
      carp(CARP_FATAL, "Error in isotope_error parameter formatting: (%s) ", isotope_errors_string.c_str());
    }
  }

  vector<int> negative_isotope_errors(1, 0);
  if (!isotope_errors_string.empty()) {
    vector<int> isotope_errors = StringUtils::Split<int>(isotope_errors_string, ',');
    for (vector<int>::iterator i = isotope_errors.begin(); i != isotope_errors.end(); i++) {
      if (*i < 0) {
        carp(CARP_FATAL, "Found a negative isotope error: %d.", *i);
      } else if (find(negative_isotope_errors.begin(), negative_isotope_errors.end(),
                      -(*i)) != negative_isotope_errors.end()) {
        carp(CARP_FATAL, "Found duplicate when parsing isotope_error parameter: %d", *i);
      }
      negative_isotope_errors.push_back(-(*i));
    }
  }
  sort(negative_isotope_errors.begin(), negative_isotope_errors.end());
  return negative_isotope_errors;
}

vector<TideSearchApplication::InputFile> TideSearchApplication::getInputFiles(
  const vector<string>& filepaths
) const {
  // Try to read all spectrum files as spectrumrecords, convert those that fail
  vector<InputFile> input_sr;
  for (vector<string>::const_iterator f = filepaths.begin(); f != filepaths.end(); f++) {
    SpectrumCollection spectra;
    pb::Header spectrum_header;
    string spectrumrecords = *f;
    bool keepSpectrumrecords = true;
    if (!spectra.ReadSpectrumRecords(spectrumrecords, &spectrum_header)) {
      // Failed, try converting to spectrumrecords file
      carp(CARP_INFO, "Converting %s to spectrumrecords format", f->c_str());
      carp(CARP_INFO, "Elapsed time starting conversion: %.3g s", wall_clock() / 1e6);
      spectrumrecords = Params::GetString("store-spectra");
      keepSpectrumrecords = !spectrumrecords.empty();
      if (!keepSpectrumrecords) {
        spectrumrecords = make_file_path(FileUtils::BaseName(*f) + ".spectrumrecords.tmp");
      } else if (filepaths.size() > 1) {
        carp(CARP_FATAL, "Cannot use store-spectra option with multiple input "
                         "spectrum files");
      }
      carp(CARP_DEBUG, "New spectrumrecords filename: %s", spectrumrecords.c_str());
      if (!SpectrumRecordWriter::convert(*f, spectrumrecords)) {
        carp(CARP_FATAL, "Error converting %s to spectrumrecords format", f->c_str());
      }
      carp(CARP_DEBUG, "Reading converted spectrum file %s", spectrumrecords.c_str());
      // Re-read converted file as spectrumrecords file
      if (!spectra.ReadSpectrumRecords(spectrumrecords, &spectrum_header)) {
        carp(CARP_DEBUG, "Deleting %s", spectrumrecords.c_str());
        FileUtils::Remove(spectrumrecords);
        carp(CARP_FATAL, "Error reading spectra file %s", spectrumrecords.c_str());
      }
    }
    input_sr.push_back(InputFile(*f, spectrumrecords, keepSpectrumrecords));
  }
  return input_sr;
}

SpectrumCollection* TideSearchApplication::loadSpectra(const string& file) {
  SpectrumCollection* spectra = new SpectrumCollection();
  pb::Header header;
  if (!spectra->ReadSpectrumRecords(file, &header)) {
    carp(CARP_FATAL, "Error reading spectrum file %s", file.c_str());
  }
  if (string_to_window_type(Params::GetString("precursor-window-type")) != WINDOW_MZ) {
    spectra->Sort();
  } else {
    spectra->Sort<ScSortByMz>(ScSortByMz(Params::GetDouble("precursor-window")));
  }
  return spectra;
}

void TideSearchApplication::search(void* threadarg) {
  struct thread_data *my_data = (struct thread_data *) threadarg;

  const string& spectrum_filename = my_data->spectrum_filename;
  const vector<SpectrumCollection::SpecCharge>* spec_charges = my_data->spec_charges;
  ActivePeptideQueue* active_peptide_queue = my_data->active_peptide_queue;
  ProteinVec& proteins = my_data->proteins;
  vector<const pb::AuxLocation*>& locations = my_data->locations;
  double precursor_window = my_data->precursor_window;
  WINDOW_TYPE_T window_type = my_data->window_type;
  double spectrum_min_mz = my_data->spectrum_min_mz;
  double spectrum_max_mz = my_data->spectrum_max_mz;
  int min_scan = my_data->min_scan;
  int max_scan = my_data->max_scan;
  int min_peaks = my_data->min_peaks;
  int search_charge = my_data->search_charge;
  int top_matches = my_data->top_matches;
  double highest_mz = my_data->highest_mz;
  ofstream* target_file = my_data->target_file;
  ofstream* decoy_file = my_data->decoy_file;
  bool compute_sp = my_data->compute_sp;
  int64_t thread_num = my_data->thread_num;
  int64_t num_threads = my_data->num_threads;
  int nAA = my_data->nAA;
  double* aaFreqN = my_data->aaFreqN;
  double* aaFreqI = my_data->aaFreqI;
  double* aaFreqC = my_data->aaFreqC;
  int* aaMass = my_data->aaMass;
  int nAARes = my_data->nAARes;
  const vector<double> dAAFreqN = *(my_data->dAAFreqN);
  const vector<double> dAAFreqI = *(my_data->dAAFreqI);
  const vector<double> dAAFreqC = *(my_data->dAAFreqC);
  const vector<double> aaMassDouble = *(my_data->dAAMass);
  const pb::ModTable mod_table = *(my_data->mod_table);
  const pb::ModTable nterm_mod_table = *(my_data->nterm_mod_table);
  const pb::ModTable cterm_mod_table = *(my_data->cterm_mod_table);
  const int numDecoys = my_data->decoysPerTarget;

  vector<boost::mutex*> locks_array = my_data->locks_array;
  vector<int>* negative_isotope_errors = my_data->negative_isotope_errors;

  double bin_width = my_data->bin_width;
  double bin_offset = my_data->bin_offset;
  bool exact_pval_search = my_data->exact_pval_search;
  map<pair<string, unsigned int>, bool>* spectrum_flag = my_data->spectrum_flag;

  int* sc_index = my_data->sc_index;
  int* total_candidate_peptides = my_data->total_candidate_peptides;

  // params
  bool peptide_centric = Params::GetBool("peptide-centric-search");
  bool use_neutral_loss_peaks = Params::GetBool("use-neutral-loss-peaks");
  bool use_flanking_peaks = Params::GetBool("use-flanking-peaks");
  int max_charge = Params::GetInt("max-precursor-charge");
  // Added by Andy Lin on 2/9/2016
  // Determines which score function to use for scoring PSMs and store in SCORE_FUNCTION enum
  SCORE_FUNCTION_T curScoreFunction = string_to_score_function_type(Params::GetString("score-function"));

  // This is the main search loop.
  ObservedPeakSet observed(bin_width, bin_offset,
                           use_neutral_loss_peaks,
                           use_flanking_peaks);

  // Keep track of observed peaks that get filtered out in various ways.
  long int num_range_skipped = 0;
  long int num_precursors_skipped = 0;
  long int num_isotopes_skipped = 0;
  long int num_retained = 0;

  // cycle through spectrum-charge pairs, sorted by neutral mass
  FLOAT_T sc_total = (FLOAT_T)spec_charges->size();
  int print_interval = Params::GetInt("print-search-progress");

  for (vector<SpectrumCollection::SpecCharge>::const_iterator sc = spec_charges->begin()+thread_num;
       sc < spec_charges->begin() + (spec_charges->size());
       sc = sc + num_threads) {
    locks_array[LOCK_REPORTING]->lock();
    ++(*sc_index);
    if (print_interval > 0 && *sc_index > 0 && *sc_index % print_interval == 0) {
      carp(CARP_INFO, "%d spectrum-charge combinations searched, %.0f%% complete",
           *sc_index, *sc_index / sc_total * 100);
    }
    locks_array[LOCK_REPORTING]->unlock();

    Spectrum* spectrum = sc->spectrum;
    double precursor_mz = spectrum->PrecursorMZ();
    double precursorMass = sc->neutral_mass;  //Added by Andy Lin (needed for residue evidence)
    int charge = sc->charge;
    int scan_num = spectrum->SpectrumNumber();
    if (spectrum_flag != NULL) {
      locks_array[LOCK_CASCADE]->lock();
      map<pair<string, unsigned int>, bool>::iterator spectrum_id;
      spectrum_id = spectrum_flag->find(pair<string, unsigned int>(
        spectrum_filename, scan_num * 10 + charge));
      if (spectrum_id != spectrum_flag->end()) {
        locks_array[LOCK_CASCADE]->unlock();
        continue;
      }
      locks_array[LOCK_CASCADE]->unlock();
    }

    if (precursor_mz < spectrum_min_mz || precursor_mz > spectrum_max_mz ||
        scan_num < min_scan || scan_num > max_scan ||
        spectrum->Size() < min_peaks ||
        (search_charge != 0 && charge != search_charge) || charge > max_charge) {
      continue;
    }
    // The active peptide queue holds the candidate peptides for spectrum.
    // Calculate and set the window, depending on the window type.
    vector<double>* min_mass = new vector<double>();
    vector<double>* max_mass = new vector<double>();
    vector<bool>* candidatePeptideStatus = new vector<bool>();
    double min_range, max_range;
    computeWindow(*sc, window_type, precursor_window, max_charge,
                  negative_isotope_errors, min_mass, max_mass, &min_range, &max_range);

    //TODO throw error when fragment-tolerance and evidence-granularity parameters are defined

    if (curScoreFunction == XCORR_SCORE && !exact_pval_search_) {  //execute original tide-search program
      // Normalize the observed spectrum and compute the cache of
      // frequently-needed values for taking dot products with theoretical
      // spectra.
      observed.PreprocessSpectrum(*spectrum, charge, &num_range_skipped,
                                  &num_precursors_skipped,
                                  &num_isotopes_skipped, &num_retained);
      int nCandPeptide = active_peptide_queue->SetActiveRange(
        min_mass, max_mass, min_range, max_range, candidatePeptideStatus);
      if (nCandPeptide == 0) {
        continue;
      }
      locks_array[LOCK_CANDIDATES]->lock();
      *total_candidate_peptides += nCandPeptide;
      locks_array[LOCK_CANDIDATES]->unlock();

      int candidatePeptideStatusSize = candidatePeptideStatus->size();
      TideMatchSet::Arr2 match_arr2(candidatePeptideStatusSize); // Scored peptides will go here.

      // Programs for taking the dot-product with the observed spectrum are laid
      // out in memory managed by the active_peptide_queue, one program for each
      // candidate peptide. The programs will store the results directly into
      // match_arr. We now pass control to those programs.
      collectScoresCompiled(active_peptide_queue, spectrum, observed, &match_arr2,
                            candidatePeptideStatusSize, charge);

      // matches will arrange the results in a heap by score, return the top
      // few, and recover the association between counter and peptide. We output
      // the top matches.
      if (peptide_centric) {
        deque<Peptide*>::const_iterator iter_ = active_peptide_queue->iter_;
        TideMatchSet::Arr2::iterator it = match_arr2.begin();
        for (; it != match_arr2.end(); ++iter_, ++it) {
          int peptide_idx = candidatePeptideStatusSize - (it->second);
          if ((*candidatePeptideStatus)[peptide_idx]) {
            (*iter_)->AddHit(spectrum, it->first, 0.0, it->second, charge);
          }
        }
      } else {  //spectrum centric match report.
        TideMatchSet::Arr match_arr(nCandPeptide);
        for (TideMatchSet::Arr2::iterator it = match_arr2.begin();
             it != match_arr2.end();
             ++it) {
          int peptide_idx = candidatePeptideStatusSize - (it->second);
          if ((*candidatePeptideStatus)[peptide_idx]) {
            TideMatchSet::Scores curScore;
            curScore.xcorr_score = (double)(it->first / XCORR_SCALING);
            curScore.rank = it->second;
            match_arr.push_back(curScore);
          }
        }

        TideMatchSet matches(&match_arr, highest_mz);
        matches.exact_pval_search_ = exact_pval_search;
        matches.cur_score_function_ = curScoreFunction;

        matches.report(target_file, decoy_file, top_matches, numDecoys, spectrum_filename,
                       spectrum, charge, active_peptide_queue, proteins,
                       locations, compute_sp, true, locks_array[LOCK_RESULTS]);
      }  //end peptide_centric == false
    } else { //This runs curScoreFunction=BOTH_SCORE, curScoreFunction=RESIUDUE_EVIDENCE_MATRIX, and xcorr p-val

      int nCandPeptide = active_peptide_queue->SetActiveRangeBIons(min_mass, max_mass, min_range, max_range, candidatePeptideStatus);
      int candidatePeptideStatusSize = candidatePeptideStatus->size();
      if (nCandPeptide == 0) {
        continue;
      }

      locks_array[LOCK_CANDIDATES]->lock();
      *total_candidate_peptides += nCandPeptide;
      locks_array[LOCK_CANDIDATES]->unlock();

      //TODO so this includes ALL amino acids seen (including modified, NTerm mod, CTerm Mod)
      //as a result -- we will look for NTerm mod amino acids throughout spectrum instead of
      //just amino acids without NTerm mod
      vector<int> aaMassInt;
      if (curScoreFunction != XCORR_SCORE) {
        for (int i = 0; i < aaMassDouble.size(); i++) {
          int tmpMass = MassConstants::mass2bin(aaMassDouble[i]);
          aaMassInt.push_back(tmpMass);
        }
      }
      int maxPrecurMassBin = floor(MaxBin::Global().CacheBinEnd() + 50.0);
      double fragTol = Params::GetDouble("fragment-tolerance");
      int granularityScale = Params::GetInt("evidence-granularity");

      //TODO look at this
      int minDeltaMass;
      int maxDeltaMass;
      if (curScoreFunction == XCORR_SCORE) {
        minDeltaMass = aaMass[0];
        maxDeltaMass = aaMass[nAA - 1];
      } else {
        minDeltaMass = aaMassInt[0];
        maxDeltaMass = aaMassInt[nAARes - 1];
      }

      TideMatchSet::Arr match_arr(nCandPeptide); // scored peptides will go here.

      // iterators needed at multiple places in following code
      deque<Peptide*>::const_iterator iter_ = active_peptide_queue->iter_;
      deque<TheoreticalPeakSetBIons>::const_iterator iter1_ = active_peptide_queue->iter1_;

      //************************************************************************
      /* For one observed spectrum, calculates:
       *  - vector of cleavage evidence
       *  - score count vectors for a range of integer masses
       *  - p-values of XCorr match scores between spectrum and all selected candidate target and decoy peptides
       * Written by Jeff Howbert, October 2013.
       * Ported to and integrated with Tide by Jeff Howbert, November 2013.
       *
       * In addition calculates:
       *   - a residue evidence matrix for a range of int masses
       *   - score count vectors for range of int masses (using res-ev matrix)
       *   - p-values of residue-evidence match between spectrum and all selected target and decoy peptides
       * Written by Jeff Howbert
       * Ported to and integrated with Tide by Andy Lin, Nov 2016
       */
      int peidx, pe, ma;
      vector<int> pepMassInt;
      pepMassInt.reserve(nCandPeptide);
      vector<int> pepMassIntUnique;
      pepMassIntUnique.reserve(nCandPeptide);

      //For each candidate peptide, determine which discretized mass bin it is in
      //pepMassInt contains the corresponding mass bin for each candidate peptide
      //pepMassIntUnique contains the unique set of mass bins that candidate peptides fall in
      getMassBin(pepMassInt, pepMassIntUnique, active_peptide_queue, candidatePeptideStatus);
      int nPepMassIntUniq = (int)pepMassIntUnique.size();

      //XCORR
      vector< vector<int> > evidenceObs(nPepMassIntUniq, vector<int>(maxPrecurMassBin, 0));
      int* scoreOffsetObs = new int[nPepMassIntUniq];
      double** pValueScoreObs = new double*[nPepMassIntUniq];
      int* intensArrayTheor = new int [maxPrecurMassBin]; // initialized later in loop
      //END XCORR

      //RES-EV
      //Creates a 3D vector representing 3D matrix
      //nPepMassIntUniq: number of mass bins candidate are in
      //nAARes: number of amino acids
      //maxPrecurMassBin: max number of mass bins
      vector<vector<vector<double> > > residueEvidenceMatrix(nPepMassIntUniq,
	      vector<vector<double> >(nAARes, vector<double>(maxPrecurMassBin, 0)));

      //Stores the score offset needed calculating res-ev p-values
      vector<int> scoreResidueOffsetObs(maxPrecurMassBin, -1);

      //For each mass bin, a vector hold the p-values for each corresponding res-ev score
      vector<vector<double> > pValuesResidueObs(maxPrecurMassBin);

      //TODO assumption is that there is one nterm mod per peptide
      int nTermMassBin;
      double nTermMass;
      if (nterm_mod_table.static_mod_size() > 0) {
        nTermMassBin = MassConstants::mass2bin(
                         MassConstants::mono_h + nterm_mod_table.static_mod(0).delta());
        nTermMass = MassConstants::mono_h + nterm_mod_table.static_mod(0).delta();
      } else {
        nTermMassBin = MassConstants::mass2bin(MassConstants::mono_h);
        nTermMass = MassConstants::mono_h;
      }

      //TODO assumption is that there is one cterm mod per peptide
      int cTermMassBin;
      double cTermMass;
      if (cterm_mod_table.static_mod_size() > 0) {
        cTermMassBin = MassConstants::mass2bin(
		      MassConstants::mono_oh + cterm_mod_table.static_mod(0).delta());
        cTermMass = MassConstants::mono_oh + cterm_mod_table.static_mod(0).delta();
      } else {
        cTermMassBin = MassConstants::mass2bin(MassConstants::mono_oh);
        cTermMass = MassConstants::mono_oh;
      }

      map<int, bool> calcDPMatrix; //for each precursor mass bin, bool determines whether to calc DP matrix
      //END RES-EV

      //Create a residue evidence matrix and evidence vector
      //for each mass bin candidate peptides are in
      for (pe = 0; pe < nPepMassIntUniq; pe++) {
        //XCORR
        if (curScoreFunction != RESIDUE_EVIDENCE_MATRIX) {
          scoreOffsetObs[pe] = 0;
          int pepMaInt = pepMassIntUnique[pe]; // TODO should be accessed with an iterator

          //preprocess to create one integerized evidence vector for each cluster of masses among selected peptides
          double pepMassMonoMean = (pepMaInt - 0.5 + bin_offset_) * bin_width_;
          evidenceObs[pe] = spectrum->CreateEvidenceVectorDiscretized(
            bin_width, bin_offset, charge, pepMassMonoMean, maxPrecurMassBin,
            &num_range_skipped, &num_precursors_skipped, &num_isotopes_skipped, &num_retained);
        }
        //END XCORR

        //RES-EV
        if (curScoreFunction != XCORR_SCORE) {
          // note: aaMassDouble differs from aaMass
          // aaMassDouble contains amino acids masses in float form
          // aaMass contains amino acid asses in integer form
          // precursorMass is the neutral mass
          observed.CreateResidueEvidenceMatrix(*spectrum, charge, maxPrecurMassBin, precursorMass,
                                               nAARes, aaMassDouble, fragTol, granularityScale,
                                               nTermMass, cTermMass,&num_range_skipped, 
                                               &num_precursors_skipped, &num_isotopes_skipped, &num_retained,
                                               residueEvidenceMatrix[pe]);
          vector<vector<double> > curResidueEvidenceMatrix = residueEvidenceMatrix[pe];

          //Get rid of values larger than curPepMassInt
          int curPepMassInt = pepMassIntUnique[pe];
          for (int i = 0; i < curResidueEvidenceMatrix.size(); i++) {
            curResidueEvidenceMatrix[i].resize(curPepMassInt);
          }
          residueEvidenceMatrix[pe] = curResidueEvidenceMatrix;
          calcDPMatrix[curPepMassInt] = false;
        }
        //END RES-Ev
      }

      //Calculates a residue evidence score and a xcorr score
      //between a spectrum and all possible peptide candidates
      //based upon the residue evidence matrix and the theoretical spectrum
      int scoreResidueEvidence;
      int scoreRefactInt;
      vector<int> resEvScores;
      vector<int> xcorrScores;
      pe = 0;
      for (peidx = 0; peidx < candidatePeptideStatusSize; peidx++) {
        if ((*candidatePeptideStatus)[peidx]) {
          int pepMassIntIdx = 0;
          int curPepMassInt;
          for (ma = 0; ma < nPepMassIntUniq; ma++ ) { //TODO should probably use iterator instead
            if (pepMassIntUnique[ma] == pepMassInt[pe]) { //TODO pepMassIntUnique should be accessed with an interator
              pepMassIntIdx = ma;
              curPepMassInt = pepMassIntUnique[ma];
              break;
            }
          }

          //XCORR
          // score XCorr for target peptide with integerized evidenceObs array
          if (curScoreFunction != RESIDUE_EVIDENCE_MATRIX) {
            for (ma = 0; ma < maxPrecurMassBin; ma++) {
              intensArrayTheor[ma] = 0;
            }

            for (vector<unsigned int>::const_iterator iter_uint = iter1_->unordered_peak_list_.begin();
                 iter_uint != iter1_->unordered_peak_list_.end();
                 iter_uint++) {
              intensArrayTheor[*iter_uint] = 1;
            }

            scoreRefactInt = 0;
            for (ma = 0; ma < maxPrecurMassBin; ma++) {
              scoreRefactInt += evidenceObs[pepMassIntIdx][ma] * intensArrayTheor[ma];
            }
            xcorrScores.push_back(scoreRefactInt);
          }
          //END XCORR

          //RES-EV
          if (curScoreFunction != XCORR_SCORE) {
            vector<vector<double> > curResidueEvidenceMatrix = residueEvidenceMatrix[pepMassIntIdx];
            Peptide* curPeptide = (*iter_);

            vector<unsigned int> intensArrayTheorResEv;
            for (vector<unsigned int>::const_iterator iter_uint = iter1_->unordered_peak_list_.begin();
                 iter_uint != iter1_->unordered_peak_list_.end();
                 iter_uint++) {
              intensArrayTheorResEv.push_back(*iter_uint);
            }

            scoreResidueEvidence = calcResEvScore(curResidueEvidenceMatrix,intensArrayTheorResEv,aaMassDouble,curPeptide);
            resEvScores.push_back(scoreResidueEvidence);

            if (scoreResidueEvidence > 0) { // if > 0, set bool to true to create DP matrix
              calcDPMatrix[curPepMassInt] = true;
            }
          }
          //END RES-EV
          pe++;
        }
        ++iter_;
        ++iter1_;
      }

      if (curScoreFunction == RESIDUE_EVIDENCE_MATRIX || curScoreFunction == BOTH_SCORE) {
        assert(resEvScores.size() == nCandPeptide);
      }
      if (curScoreFunction == XCORR_SCORE || curScoreFunction == BOTH_SCORE) {
        assert(xcorrScores.size() == nCandPeptide);
      }

      //XCORR
      //Create a dynamic programming vector is there is a xcorr
      //and if user specified as a score function either 'xcorr' or 'both'
      if (curScoreFunction != RESIDUE_EVIDENCE_MATRIX) {
        for (pe = 0; pe < nPepMassIntUniq; pe++) { // TODO should probably instead use iterator over pepMassIntUnique
          int pepMaInt = pepMassIntUnique[pe]; // TODO should be accessed with an iterator

          // NOTE: will have to go back to separate dynamic programming for
          //       target and decoy if they have different probNI and probC
          int maxEvidence = *std::max_element(evidenceObs[pe].begin(), evidenceObs[pe].end());
          int minEvidence = *std::min_element(evidenceObs[pe].begin(), evidenceObs[pe].end());

          // estimate maxScore and minScore
          int maxNResidue = (int)floor((double)pepMaInt / (double)minDeltaMass);
          vector<int> sortEvidenceObs(evidenceObs[pe].begin(), evidenceObs[pe].end());
          std::sort(sortEvidenceObs.begin(), sortEvidenceObs.end(), greater<int>());
          int maxScore = 0;
          int minScore = 0;
          for (int sc = 0; sc < maxNResidue; sc++) {
            maxScore += sortEvidenceObs[sc];
          }

          for (int sc = maxPrecurMassBin - maxNResidue; sc < maxPrecurMassBin; sc++) {
            minScore += sortEvidenceObs[sc];
          }

          int bottomRowBuffer = maxEvidence + 1;
          int topRowBuffer = -minEvidence;
          int nRowDynProg = bottomRowBuffer - minScore + 1 + maxScore + topRowBuffer;
          pValueScoreObs[pe] = new double[nRowDynProg];

          scoreOffsetObs[pe] = calcScoreCount(maxPrecurMassBin, &evidenceObs[pe][0], pepMaInt,
                               maxEvidence, minEvidence, maxScore, minScore,
                               nAA, aaFreqN, aaFreqI, aaFreqC, aaMass,
                               pValueScoreObs[pe]);
        }
      }
      //END XCORR

      //RES-EV
      //Create dyanamic programming matrix if there is a res-ev score greater than 0
      //and if user specified as a score function either 'residue-evidence matrix' or 'both'
      if (curScoreFunction != XCORR_SCORE) {
        for (pe=0 ; pe<nPepMassIntUniq ; pe++) {
          int curPepMassInt = pepMassIntUnique[pe];
          if (calcDPMatrix[curPepMassInt] == false) {
            continue;
          }

          vector<vector<double> > curResidueEvidenceMatrix = residueEvidenceMatrix[pe];
          vector<int> maxColEvidence(curPepMassInt,0);

          //maxColEvidence is edited by reference
          int maxEvidence = getMaxColEvidence(curResidueEvidenceMatrix,maxColEvidence,curPepMassInt);
          int maxNResidue = floor((double)curPepMassInt / 57.0);

          std::sort(maxColEvidence.begin(),maxColEvidence.end(),greater<int>());
          int maxScore = 0;
          for(int i = 0; i < maxNResidue; i++) { //maxColEvidence has been sorted
            maxScore += maxColEvidence[i];
          }

          int scoreOffset;
          vector<double> scoreResidueCount;

          calcResidueScoreCount(nAARes,curPepMassInt,curResidueEvidenceMatrix,aaMassInt,
                                dAAFreqN, dAAFreqI, dAAFreqC,nTermMassBin,cTermMassBin,
                                minDeltaMass,maxDeltaMass,maxEvidence,maxScore,
                                scoreResidueCount,scoreOffset);
          scoreResidueOffsetObs[curPepMassInt] = scoreOffset;

          double totalCount = 0;
          for (int i=scoreOffset ; i<scoreResidueCount.size() ; i++) {
            totalCount += scoreResidueCount[i];
          }
          for (int i=scoreResidueCount.size()-2 ; i>-1; i--) {
            scoreResidueCount[i] = scoreResidueCount[i] + scoreResidueCount[i+1];
          }
          for (int i = 0; i < scoreResidueCount.size(); i++) {
            //Avoid potential underflow
            scoreResidueCount[i] = exp(log(scoreResidueCount[i]) - log(totalCount));
          }
          pValuesResidueObs[curPepMassInt] = scoreResidueCount;
        }
      }
      //END RES-EV

      /************ calculate p-values for PSMs using residue evidence matrix ****************/
      iter_ = active_peptide_queue->iter_;
      iter1_ = active_peptide_queue->iter1_;
      int curPepMassInt;
      double pValue_xcorr;
      double pValue_resEv;
      double pValue_both;
      pe = 0;
      for (peidx = 0; peidx < candidatePeptideStatusSize; peidx++) {
        if ((*candidatePeptideStatus)[peidx]) {
          int pepMassIntIdx = 0;

          for (ma = 0; ma < nPepMassIntUniq; ma++ ) { //TODO should probably use iterator instead
            if (pepMassIntUnique[ma] == pepMassInt[pe]) { //TODO pepMassIntUnique should be accessed with an interator
              pepMassIntIdx = ma;
              curPepMassInt = pepMassIntUnique[ma];
              break;
            }
          }

          int scoreCountIdx;
          //XCORR
          if (curScoreFunction != RESIDUE_EVIDENCE_MATRIX) {
            scoreRefactInt = xcorrScores[pe];
            scoreCountIdx = scoreRefactInt + scoreOffsetObs[pepMassIntIdx];
            pValue_xcorr = pValueScoreObs[pepMassIntIdx][scoreCountIdx];
          }
          //END XCORR

          //RES-EV
          if (curScoreFunction != XCORR_SCORE) {
            scoreResidueEvidence = resEvScores[pe];
            if (calcDPMatrix[curPepMassInt]) {
              scoreCountIdx = scoreResidueEvidence + scoreResidueOffsetObs[curPepMassInt];
              pValue_resEv = pValuesResidueObs[curPepMassInt][scoreCountIdx];
            } else {
              pValue_resEv = 1.0;
            }
          }
          //END RES-EV

          //BOTH SCORE
          if (curScoreFunction == BOTH_SCORE) {
            double cPval = pValue_xcorr * pValue_resEv;

            double m = 1.2; // This value has been empircally determined
            pValue_both = calcCombinedPval(m,cPval,2); //2 is the # of p-values that are combined
          }
          //END BOTH_SCORE

          if (curScoreFunction == XCORR_SCORE && pValue_xcorr == 0.0) {
            carp(CARP_FATAL, "PSM p-value should not be equal to 0.0");
          } else if (curScoreFunction == RESIDUE_EVIDENCE_MATRIX && pValue_resEv == 0.0) {
            carp(CARP_FATAL, "PSM p-value should not be equal to 0.0");
          } else if (curScoreFunction == BOTH_SCORE && pValue_both == 0.0) {
            carp(CARP_FATAL, "PSM p-value should not be equal to 0.0");
          }

          if (peptide_centric) {
            carp(CARP_FATAL, "residue-evidence has not been implemented with 'peptide-centric-search T' yet.");
          } else {
            TideMatchSet::Scores curScore;
            curScore.xcorr_score = (double)scoreRefactInt / RESCALE_FACTOR;
            curScore.xcorr_pval = pValue_xcorr;
            curScore.resEv_pval = pValue_resEv;
            curScore.resEv_score = scoreResidueEvidence;
            curScore.combinedPval = pValue_both;
            //TODO ugly hack to conform with the way these indices are generated in standard tide-search
            curScore.rank = candidatePeptideStatusSize - peidx;
            match_arr.push_back(curScore);
          }
          pe++;
        }
        ++iter_;
        ++iter1_;
      }

      //clean up
      if (curScoreFunction != RESIDUE_EVIDENCE_MATRIX) {
        for (pe = 0; pe < nPepMassIntUniq; pe++) {
          delete [] pValueScoreObs[pe];
        }
      }
      delete [] scoreOffsetObs;
      delete [] pValueScoreObs;
      delete [] intensArrayTheor;

      if (!peptide_centric) {
        // below text is copied from text above in the exact-p-value XCORR case
        // matches will arrange the results in a heap by score, return the top
        // few, and recover the association between counter and peptide. We output
        // the top matches.
        TideMatchSet matches(&match_arr, highest_mz);
        matches.exact_pval_search_ = exact_pval_search_;
        matches.cur_score_function_ = curScoreFunction;

        if (curScoreFunction == RESIDUE_EVIDENCE_MATRIX && exact_pval_search_ == false) {
          matches.report(target_file, decoy_file, top_matches, numDecoys, spectrum_filename,
                         spectrum, charge, active_peptide_queue, proteins,
                         locations, compute_sp, true, locks_array[LOCK_RESULTS]);
        } else {
          matches.report(target_file, decoy_file, top_matches, numDecoys, spectrum_filename,
                         spectrum, charge, active_peptide_queue, proteins,
                         locations, compute_sp, false, locks_array[LOCK_RESULTS]);
        }
      } //end peptide_centric == false
    }
    delete min_mass;
    delete max_mass;
    delete candidatePeptideStatus;
  }

  if (!Params::GetBool("skip-preprocessing")) {
    locks_array[LOCK_REPORTING]->lock();
    if (curScoreFunction == BOTH_SCORE) {
      num_precursors_skipped = num_precursors_skipped / 2;
      num_isotopes_skipped = num_isotopes_skipped / 2;
      num_range_skipped = num_range_skipped / 2;
      num_retained = num_retained / 2;
    }

    long int total_peaks = num_precursors_skipped + num_isotopes_skipped + num_range_skipped + num_retained;
    if (total_peaks == 0) {
      carp(CARP_INFO, "[Thread %d]: Warning: no peaks found.", thread_num);
    } else {
      carp(CARP_INFO,
           "[Thread %d]: Deleted %d precursor, %d isotope and %d out-of-range peaks.",
           thread_num, num_precursors_skipped, num_isotopes_skipped, num_range_skipped);
    }
    if (num_retained == 0) {
      carp(CARP_INFO, "[Thread %d]: Warning: no peaks retained.", thread_num);
    } else {
      carp(CARP_INFO, "[Thread %d]: Retained %g%% of peaks.",
           thread_num, (100.0 * num_retained) / total_peaks);
    }
    locks_array[LOCK_REPORTING]->unlock();
  }
}

void TideSearchApplication::search(
  const string& spectrum_filename,
  const vector<SpectrumCollection::SpecCharge>* spec_charges,
  vector<ActivePeptideQueue*> active_peptide_queue,
  ProteinVec& proteins,
  vector<const pb::AuxLocation*>& locations,
  double precursor_window,
  WINDOW_TYPE_T window_type,
  double spectrum_min_mz,
  double spectrum_max_mz,
  int min_scan,
  int max_scan,
  int min_peaks,
  int search_charge,
  int top_matches,
  double highest_mz,
  ofstream* target_file,
  ofstream* decoy_file,
  bool compute_sp,
  int nAA,
  double* aaFreqN,
  double* aaFreqI,
  double* aaFreqC,
  int* aaMass,
  int nAARes,
  const vector<double>& dAAFreqN,
  const vector<double>& dAAFreqI,
  const vector<double>& dAAFreqC,
  const vector<double>& dAAMass,
  const pb::ModTable& mod_table,
  const pb::ModTable& nterm_mod_table,
  const pb::ModTable& cterm_mod_table,
  int numDecoys,
  vector<int>* negative_isotope_errors
) {
  // Create an array of locks.
  vector<boost::mutex *> locks_array;
  for (int i = 0; i < NUMBER_LOCK_TYPES; i++) {
    locks_array.push_back(new boost::mutex());
  }

  int elution_window = Params::GetInt("elution-window-size");
  bool peptide_centric = Params::GetBool("peptide-centric-search");

  // initialize fields required for output
  int* sc_index = new int(-1);
  int* total_candidate_peptides = new int(0);
  FLOAT_T sc_total = (FLOAT_T)spec_charges->size();

  if (peptide_centric == false) {
    elution_window = 0;
  }

  for (int i = 0; i < NUM_THREADS; i++) {
    active_peptide_queue[i]->setElutionWindow(elution_window);
    active_peptide_queue[i]->setPeptideCentric(peptide_centric);
  }

  if (elution_window > 0 && elution_window % 2 == 0) {
    for (int i = 0; i < NUM_THREADS; i++) {
      active_peptide_queue[i]->setElutionWindow(elution_window+1);
    }
  }

  if (!peptide_centric || !exact_pval_search_) {
    for (int i = 0; i < NUM_THREADS; i++) {
      active_peptide_queue[i]->setElutionWindow(0);
    }
  }

  for (int i = 0; i < NUM_THREADS; i++) {
    active_peptide_queue[i]->SetOutputs(
      NULL, &locations, top_matches, compute_sp, target_file, decoy_file, highest_mz);
  }

  // Creating structs to hold information required for each thread to search through
  // a spec charge

  vector<thread_data> thread_data_array;
  for (int i= 0; i < NUM_THREADS; i++) {
      thread_data_array.push_back(thread_data(spectrum_filename, spec_charges, active_peptide_queue[i],
      proteins, locations, precursor_window, window_type, spectrum_min_mz,
      spectrum_max_mz, min_scan, max_scan, min_peaks, search_charge, top_matches,
      highest_mz, target_file, decoy_file, compute_sp,
      i, NUM_THREADS, nAA, aaFreqN, aaFreqI, aaFreqC, aaMass,
      nAARes, &dAAFreqN, &dAAFreqI, &dAAFreqC, &dAAMass,
      &mod_table, &nterm_mod_table, &cterm_mod_table, numDecoys, locks_array, //TODO do I need to delete pointer somewhere?
      bin_width_, bin_offset_, exact_pval_search_, spectrum_flag_, sc_index, total_candidate_peptides, negative_isotope_errors));
  }

  boost::thread_group threadgroup;

  // Launch threads
  for (int64_t t = 1; t < NUM_THREADS; t++) {
    boost::thread * currthread = new boost::thread(boost::bind(&TideSearchApplication::search, this, (void *) &(thread_data_array[t])));
    threadgroup.add_thread(currthread);
  }

  // Searches through part of the spec charge vector while waiting for threads are busy
  search( (void *) &(thread_data_array[0]) );

  // Join threads
  threadgroup.join_all();

  carp(CARP_INFO, "Time per spectrum-charge combination: %lf s.", wall_clock() / (1e6*sc_total));
  carp(CARP_INFO, "Average number of candidates per spectrum-charge combination: %lf ",
                  (*total_candidate_peptides) / sc_total);
  for (int i = 0; i < NUMBER_LOCK_TYPES; i++) {
    delete locks_array[i];
  }
  delete sc_index;
  delete total_candidate_peptides;

}

#ifdef _WIN64
#pragma optimize( "g", off )
#endif
void TideSearchApplication::collectScoresCompiled(
  ActivePeptideQueue* active_peptide_queue,
  const Spectrum* spectrum,
  const ObservedPeakSet& observed,
  TideMatchSet::Arr2* match_arr,
  int queue_size,
  int charge
) {
  if (!active_peptide_queue->HasNext()) {
    return;
  }
  // prog gets the address of the dot-product program for the first peptide
  // in the active queue.
  const void* prog = active_peptide_queue->NextPeptide()->Prog(charge);
  const int* cache = observed.GetCache();
  // results will get (score, counter) pairs, where score is the dot product
  // of the observed peak set with a candidate peptide. The candidate
  // peptide is given by counter, which refers to the index within the
  // ActivePeptideQueue, counting from the back. This complication
  // simplifies the generated programs, which now simply dump the counter.
  pair<int, int>* results = match_arr->data();

  // See compiler.h for a description of the programs beginning at prog and
  // how they are generated. Here we initialize certain registers to the
  // values expected by the programs and call the first one (*prog).
  //
  // See gnu assembler format for more on this format. We tell the compiler
  // to set these registers:
  // edx/rdx points to the cache.
  // eax/rax points to the first program.
  // ecx/rcx is the counter and gets the size of the active queue.
  // edi/rdi points to the results buffer.
  //
  // The push and pop operations are a workaround for a compiler that
  // doesn't understand that %ecx and %edi (or %rcx and %rdi) get
  // clobbered. Since they're already input registers, they can't be
  // included in the clobber list.

#ifdef _MSC_VER
#ifdef _WIN64
  DWORD64 rcx;
  DWORD64 rdi;
  volatile bool restored = false;
  CONTEXT context;
  RtlCaptureContext(&context);
  if (!restored) {
    rcx = context.Rcx;
    rdi = context.Rdi;

    context.Rdx = (DWORD64)cache;
    context.Rax = (DWORD64)prog;
    context.Rcx = (DWORD64)queue_size;
    context.Rdi = (DWORD64)results;

    restored = true;
    RtlRestoreContext(&context, NULL);
  } else {
    ((void(*)(void))prog)();
  }

  restored = false;
  RtlCaptureContext(&context);
  if (!restored) {
    context.Rcx = rcx;
    context.Rdi = rdi;

    restored = true;
    RtlRestoreContext(&context, NULL);
  }
#else
  __asm {
    cld
    push ecx
    push edi
    mov edx, cache
    mov eax, prog
    mov ecx, queue_size
    mov edi, results
    call eax
    pop edi
    pop ecx
  }
#endif
#else
  __asm__ __volatile__("cld\n" // stos operations increment edi
#ifdef __x86_64__
                       "push %%rcx\n"
                       "push %%rdi\n"
                       "call *%%rax\n"
                       "pop %%rdi\n"
                       "pop %%rcx\n"
#else
                       "push %%ecx\n"
                       "push %%edi\n"
                       "call *%%eax\n"
                       "pop %%edi\n"
                       "pop %%ecx\n"
#endif
                       : // no outputs
                       : "d" (cache),
                         "a" (prog),
                         "c" (queue_size),
                         "D" (results)
  );
#endif

  // match_arr is filled by the compiled programs, not by calls to
  // push_back(). We have to set the final size explicitly.
  match_arr->set_size(queue_size);
}
#ifdef _WIN64
#pragma optimize( "g", on )
#endif

void TideSearchApplication::convertResults() const {
  PSMConvertApplication converter;
  if (!Params::GetBool("concat")) {
    string target_file_name = make_file_path("tide-search.target.txt");
    if (Params::GetBool("pin-output")) {
      converter.convertFile("tsv", "pin", target_file_name, "tide-search.target.", Params::GetString("protein-database"), true);
    }
    if (Params::GetBool("pepxml-output")) {
      converter.convertFile("tsv", "pepxml", target_file_name, "tide-search.target.", Params::GetString("protein-database"), true);
    }
    if (Params::GetBool("mzid-output")) {
      converter.convertFile("tsv", "mzidentml", target_file_name, "tide-search.target.", Params::GetString("protein-database"), true);
    }
    if (Params::GetBool("sqt-output")) {
      converter.convertFile("tsv", "sqt", target_file_name, "tide-search.target.", Params::GetString("protein-database"), true);
    }

    if (HAS_DECOYS) {
      string decoy_file_name = make_file_path("tide-search.decoy.txt");
      if (Params::GetBool("pin-output")) {
        converter.convertFile("tsv", "pin", decoy_file_name, "tide-search.decoy.", Params::GetString("protein-database"), true);
      }
      if (Params::GetBool("pepxml-output")) {
        converter.convertFile("tsv", "pepxml", decoy_file_name, "tide-search.decoy.", Params::GetString("protein-database"), true);
      }
      if (Params::GetBool("mzid-output")) {
        converter.convertFile("tsv", "mzidentml", decoy_file_name, "tide-search.decoy.", Params::GetString("protein-database"), true);
      }
      if (Params::GetBool("sqt-output")) {
        converter.convertFile("tsv", "sqt", decoy_file_name, "tide-search.decoy.", Params::GetString("protein-database"), true);
      }
    }
  } else {
    string concat_file_name = make_file_path("tide-search.txt");
    if (Params::GetBool("pin-output")) {
      converter.convertFile("tsv", "pin", concat_file_name, "tide-search.", Params::GetString("protein-database"), true);
    }
    if (Params::GetBool("pepxml-output")) {
      converter.convertFile("tsv", "pepxml", concat_file_name, "tide-search.", Params::GetString("protein-database"), true);
    }
    if (Params::GetBool("mzid-output")) {
      converter.convertFile("tsv", "mzidentml", concat_file_name, "tide-search.", Params::GetString("protein-database"), true);
    }
    if (Params::GetBool("sqt-output")) {
      converter.convertFile("tsv", "sqt", concat_file_name, "tide-search.", Params::GetString("protein-database"), true);
    }
  }
}

void TideSearchApplication::computeWindow(
  const SpectrumCollection::SpecCharge& sc,
  WINDOW_TYPE_T window_type,
  double precursor_window,
  int max_charge,
  vector<int>* negative_isotope_errors,
  vector<double>* out_min,
  vector<double>* out_max,
  double* min_range,
  double* max_range
) {
  double unit_dalton = BIN_WIDTH;
  switch (window_type) {
  case WINDOW_MASS:
    for (vector<int>::const_iterator ie = negative_isotope_errors->begin(); ie != negative_isotope_errors->end(); ++ie) {
      out_min->push_back(sc.neutral_mass + (*ie * unit_dalton) - precursor_window);
      out_max->push_back(sc.neutral_mass + (*ie * unit_dalton) + precursor_window);
    }
    *min_range = (sc.neutral_mass + (negative_isotope_errors->front() * unit_dalton)) - precursor_window;
    *max_range = (sc.neutral_mass + (negative_isotope_errors->back() * unit_dalton)) + precursor_window;
    break;
  case WINDOW_MZ: {
    double mz_minus_proton = sc.spectrum->PrecursorMZ() - MASS_PROTON;
    for (vector<int>::const_iterator ie = negative_isotope_errors->begin(); ie != negative_isotope_errors->end(); ++ie) {
      out_min->push_back((mz_minus_proton - precursor_window) * sc.charge + (*ie * unit_dalton));
      out_max->push_back((mz_minus_proton + precursor_window) * sc.charge + (*ie * unit_dalton));
    }
    *min_range = (mz_minus_proton*sc.charge + (negative_isotope_errors->front() * unit_dalton)) - precursor_window*max_charge;
    *max_range = (mz_minus_proton*sc.charge + (negative_isotope_errors->back() * unit_dalton)) + precursor_window*max_charge;
    break;
  }
  case WINDOW_PPM: {
    double tiny_precursor = precursor_window * 1e-6;
    for (vector<int>::const_iterator ie = negative_isotope_errors->begin(); ie != negative_isotope_errors->end(); ++ie) {
      out_min->push_back((sc.neutral_mass + (*ie * unit_dalton)) * (1.0 - tiny_precursor));
      out_max->push_back((sc.neutral_mass + (*ie * unit_dalton)) * (1.0 + tiny_precursor));
    }
    *min_range = (sc.neutral_mass + (negative_isotope_errors->front() * unit_dalton)) * (1.0 - tiny_precursor);
    *max_range = (sc.neutral_mass + (negative_isotope_errors->back() * unit_dalton)) * (1.0 + tiny_precursor);
    break;
  }
  default:
    carp(CARP_FATAL, "Invalid window type");
  }
  carp(CARP_DETAILED_DEBUG, "Scan=%d Charge=%d Mass window=[%f, %f]",
       sc.spectrum->SpectrumNumber(), sc.charge, (*out_min)[0], (*out_max)[0]);
}

bool TideSearchApplication::proteinLevelDecoys() {
  return PROTEIN_LEVEL_DECOYS;
}

string TideSearchApplication::getName() const {
  return "tide-search";
}

string TideSearchApplication::getDescription() const {
  return
    "[[nohtml:Search a collection of spectra against a sequence database, "
    "returning a collection of peptide-spectrum matches (PSMs). This is a "
    "fast search engine but requires that you first build an index with "
    "tide-index.]]"
    "[[html:<p>Tide is a tool for identifying peptides from tandem mass "
    "spectra. It is an independent reimplementation of the SEQUEST<sup>&reg;"
    "</sup> algorithm, which assigns peptides to spectra by comparing the "
    "observed spectra to a catalog of theoretical spectra derived from a "
    "database of known proteins. Tide's primary advantage is its speed. Our "
    "published paper provides more detail on how Tide works. If you use Tide "
    "in your research, please cite:</p><blockquote>Benjamin J. Diament and "
    "William Stafford Noble. <a href=\"http://dx.doi.org/10.1021/pr101196n\">"
    "&quot;Faster SEQUEST Searching for Peptide Identification from Tandem "
    "Mass Spectra&quot;</a>. <em>Journal of Proteome Research</em>. "
    "10(9):3871-9, 2011.</blockquote> "
    "<p>When <code>tide-search</code> runs, it performs "
    "several intermediate steps, as follows:</p><ol>"
    "<li>If a FASTA file was provided, convert it to an index using "
    "<code>tide-index</code>.</li>"
    "<li>Convert the given "
    "fragmentation spectra to a binary format.</li><li>Search the spectra "
    "against the database and store the results in binary format.</li><li>"
    "Convert the results to one or more requested output formats.</li></ol><p>"
    "By default, the intermediate binary files are stored in the output "
    "directory and deleted when Tide finishes execution. If you plan to search "
    "against given database more than once or search a given set of spectra "
    "more than once, then you can direct Tide to save the binary spectrum "
    "files using the <code>--store-index</code> and "
    "<code>--store-spectra</code> options. "
    "Subsequent runs of the program will go faster "
    "if provided with inputs in binary format.</p>]]";
}

vector<string> TideSearchApplication::getArgs() const {
  string arr[] = {
    "tide spectra file+",
    "tide database"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> TideSearchApplication::getOptions() const {
  string arr[] = {
    "auto-mz-bin-width",
    "auto-precursor-window",
    "compute-sp",
    "concat",
    "deisotope",
    "elution-window-size",
    "exact-p-value",
    "file-column",
    "fileroot",
    "isotope-error",
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
    "peptide-centric-search",
    "score-function",
    "fragment-tolerance",
    "evidence-granularity",
    "pepxml-output",
    "pin-output",
    "pm-charge",
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
    "sqt-output",
    "store-index",
    "store-spectra",
    "top-match",
    "txt-output",
    "use-flanking-peaks",
    "use-neutral-loss-peaks",
    "use-z-line",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector< pair<string, string> > TideSearchApplication::getOutputs() const {
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
bool TideSearchApplication::needsOutputDirectory() const {
  return true;
}

COMMAND_T TideSearchApplication::getCommand() const {
  return TIDE_SEARCH_COMMAND;
}

/* Calculates counts of peptides with various XCorr scores, given a preprocessed
 * MS2 spectrum, using dynamic programming.
 * Written by Jeff Howbert, October, 2012 (as function calcScoreCount).
 * Ported to and integrated with Tide by Jeff Howbert, November, 2013.
 */
int TideSearchApplication::calcScoreCount(
  int numelEvidenceObs,
  int* evidenceObs,
  int pepMassInt,
  int maxEvidence,
  int minEvidence,
  int maxScore,
  int minScore,
  int nAA,
  double* aaFreqN,
  double* aaFreqI,
  double* aaFreqC,
  int* aaMass,
  double* pValueScoreObs
) {
  const int nDeltaMass = nAA;
  int minDeltaMass = aaMass[0];
  int maxDeltaMass = aaMass[nDeltaMass - 1];

  // internal variables
  int row;
  int col;
  int ma;
  int evidence;
  int de;
  int evidenceRow;
  double sumScore;

  int bottomRowBuffer = maxEvidence + 1;
  int topRowBuffer = -minEvidence;
  int colBuffer = maxDeltaMass;
  int colStart = MassConstants::mass2bin(MassConstants::mono_h);
  int scoreOffsetObs = bottomRowBuffer - minScore;

  int nRow = bottomRowBuffer - minScore + 1 + maxScore + topRowBuffer;
  int nCol = colBuffer + pepMassInt;
  int rowFirst = bottomRowBuffer;
  int rowLast = rowFirst - minScore + maxScore;
  int colFirst = colStart + MassConstants::mass2bin(MassConstants::mono_h);
  int colLast = MassConstants::mass2bin(MassConstants::bin2mass(pepMassInt)
    - MassConstants::mono_oh);
  int initCountRow = bottomRowBuffer - minScore;
  int initCountCol = maxDeltaMass + colStart;

  double** dynProgArray = new double*[nRow];
  for (row = 0; row < nRow; row++) {
    dynProgArray[row] = new double[nCol];
    for (col = 0; col < nCol; col++) {
      dynProgArray[row][col] = 0.0;
    }
  }
  double* scoreCountBinAdjust = 0;
  scoreCountBinAdjust = new double[nRow];
  for (row = 0; row < nRow; row++) {
    scoreCountBinAdjust[row] = 0.0;
  }

  dynProgArray[initCountRow][initCountCol] = 1.0; // initial count of peptides with mass = 1
  vector<int> deltaMassCol(nDeltaMass);
  // populate matrix with scores for first (i.e. N-terminal) amino acid in sequence
  for (de = 0; de < nDeltaMass; de++) {
    ma = aaMass[de];
    row = initCountRow + evidenceObs[ma + colStart];
    col = initCountCol + ma;
    if (col <= maxDeltaMass + colLast) {
      dynProgArray[row][col] += dynProgArray[initCountRow][initCountCol] * aaFreqN[de];
    }
  }
  // set to zero now that score counts for first amino acid are in matrix
  dynProgArray[initCountRow][initCountCol] = 0.0;
  // populate matrix with score counts for non-terminal amino acids in sequence
  for (ma = colFirst; ma < colLast; ma++) {
    col = maxDeltaMass + ma;
    evidence = evidenceObs[ma];
    for (de = 0; de < nDeltaMass; de++) {
      deltaMassCol[de] = col - aaMass[de];
    }
    for (row = rowFirst; row <= rowLast; row++) {
      evidenceRow = row - evidence;
      sumScore = dynProgArray[row][col];
      for (de = 0; de < nDeltaMass; de++) {
        sumScore += dynProgArray[evidenceRow][deltaMassCol[de]] * aaFreqI[de];
      }
      dynProgArray[row][col] = sumScore;
    }
  }
  // populate matrix with score counts for last (i.e. C-terminal) amino acid in sequence
  ma = colLast;
  col = maxDeltaMass + ma;
  evidence = 0; // no evidence should be added for last amino acid in sequence
  for (de = 0; de < nDeltaMass; de++) {
    deltaMassCol[de] = col - aaMass[de];
  }
  for (row = rowFirst; row <= rowLast; row++) {
    evidenceRow = row - evidence;
    sumScore = 0.0;
    for (de = 0; de < nDeltaMass; de++) {
      sumScore += dynProgArray[evidenceRow][deltaMassCol[de]] * aaFreqC[de];  // C-terminal residue
    }
    dynProgArray[row][col] = sumScore;
  }

  int colScoreCount = maxDeltaMass + colLast;
  double totalCount = 0.0;
  for (row = 0; row < nRow; row++) {
    // at this point pValueScoreObs just holds counts from last column of dynamic programming array
    pValueScoreObs[row] = dynProgArray[row][colScoreCount];
    totalCount += pValueScoreObs[row];
    scoreCountBinAdjust[row] = pValueScoreObs[row] / 2.0;
  }
  // convert from counts to cumulative sum of counts
  for (row = nRow - 2; row >= 0; row--) {
    pValueScoreObs[row] += pValueScoreObs[row + 1];
  }
  double logTotalCount = log(totalCount);
  for (row = 0; row < nRow; row++) {
    // adjust counts to reflect center of bin, not edge
    pValueScoreObs[row] -= scoreCountBinAdjust[row];
    // normalize distribution; use exp( log ) to avoid potential underflow
    pValueScoreObs[row] = exp(log(pValueScoreObs[row]) - logTotalCount);
  }

  // clean up
  for (row = 0; row < nRow; row++) {
    delete [] dynProgArray[row];
  }
  delete [] dynProgArray;
  delete [] scoreCountBinAdjust;

  return scoreOffsetObs;
}

/*
 * Calculates counts of peptides with various residue evidence scores, given
 * a preprocessed residue evidence matrix, using dynamic programming
 *
 * This version:
 *  - calculates most of dimension and inexing variables required or
 *    dynamic programming inside of function, instead of externally
 *    in MATLAB
 *  - uses uniform amino acid probabilities for all positions in peptide
 *
 * This used to be a MEX-file to be called in MATLAB
 *  - has been incorporated into Tide/Crux
 *
 * Written by Jeff Howbert, August, 2015.
 *
 * Added by Andy Lin, March 2-16
 * Edited to work within Crux code instead of with original MATLAB code
 */
void TideSearchApplication::calcResidueScoreCount (
  int nAa,
  int pepMassInt,
  vector<vector<double> >& residueEvidenceMatrix,
  vector<int>& aaMass,
  const vector<double>& aaFreqN,
  const vector<double>& aaFreqI,
  const vector<double>& aaFreqC,
  int nTermMass, //this is nTermMassBin
  int cTermMass, //this is cTermMassBin
  int minAaMass,
  int maxAaMass,
  int maxEvidence,
  int maxScore,
  vector<double>& scoreCount, //this is returned for later use
  int& scoreOffset //this is returned for later use
) {
  int minEvidence  = 0;
  int minScore     = 0;

  int row;
  int col;
  int ma;
  int evid;
  int de;
  int evidRow;
  double sumScore;

  int bottomRowBuffer = maxEvidence;
  int topRowBuffer = -minEvidence;
  int colBuffer = maxAaMass;
  int colStart = nTermMass;
  int nRow = bottomRowBuffer - minScore + 1 + maxScore + topRowBuffer;
  int nCol = colBuffer + pepMassInt;
  int rowFirst = bottomRowBuffer + 1;
  int rowLast = rowFirst - minScore + maxScore;
  int colFirst = colStart + 1;
  int colLast = pepMassInt - cTermMass;
  int initCountRow = bottomRowBuffer - minScore + 1;
  int initCountCol = maxAaMass + colStart;

  // convert to zero-based indexing
  rowFirst = rowFirst - 1;
  rowLast = rowLast - 1;
  colFirst = colFirst - 1;
  colLast = colLast - 1;
  initCountRow = initCountRow - 1;
  initCountCol = initCountCol - 1;

  double** dynProgArray = 0;
  dynProgArray = new double*[nRow];
  for (row = 0; row < nRow; row++) {
    dynProgArray[row] = new double[nCol];
    for (col = 0; col < nCol; col++) {
      dynProgArray[row][col] = 0.0;
    }
  }

  // initial count of peptides with mass = nTermMass
  dynProgArray[initCountRow][initCountCol] = 1.0;

  int* aaMassCol = new int[nAa];
  // populate matrix with scores for first (i.e. N-terminal) amino acid in sequence
  for (de = 0; de < nAa; de++) {
    ma = aaMass[de];

    //&& -1 is to account for zero-based indexing in evidence vector
    //row = initCountRow + residueEvidueMatrix[ de ][ ma + nTermMass - 1 ]; //original
    row = initCountRow + residueEvidenceMatrix[de][ma + 1 - 1]; //+1 for N-Term H and -1 for 0 indexing

    //TODO need to change this to based off bool
    if (nTermMass == 1) { //N-Term not modified
      col = initCountCol + ma;
    } else { //N-Term is modified
      col = initCountCol + ma - nTermMass + 1;
    }

//    if ( col <= maxAaMass + colLast ) { //original
    if (col <= maxAaMass + colLast && col >= initCountCol) { //TODO not sure if below or above is correct
      //dynProgArray[ row ][ col ] += dynProgArray[ initCountRow ][ initCountCol ];
      dynProgArray[row][col] += dynProgArray[initCountRow][initCountCol] * aaFreqN[de];
    }
  }

  //set to zero now that score counts for first amino acid are in matrix
  dynProgArray[initCountRow][initCountCol] = 0.0;

  // populate matrix with score counts for non-terminal amino acids in sequence
  for (ma = colFirst; ma < colLast; ma++) {
    col = maxAaMass + ma;

    for (de = 0; de < nAa; de++) {
      aaMassCol[de] = col - aaMass[de];
    }
    for (row = rowFirst; row <= rowLast; row++) {
      sumScore = dynProgArray[row][col];
      for (de = 0; de < nAa; de++) {
        evidRow = row - residueEvidenceMatrix[de][ma];
        //sumScore += dynProgArray[ evidRow ][ aaMassCol[ de ] ];
        sumScore += dynProgArray[evidRow][aaMassCol[de]] * aaFreqI[de];
      }
      dynProgArray[row][col] = sumScore;
    }
  }

  // populate matrix with score counts for last (i.e. C-terminal) amino acid in sequence
  ma = colLast;
  col = maxAaMass + ma;

  //no evidence should be added for last amino acid in sequence
  evid = 0;
  for (de = 0; de < nAa; de++) {
    aaMassCol[de] = col - aaMass[de];
  }
  for (row = rowFirst; row <= rowLast; row++) {
    evidRow = row - evid;
    sumScore = 0.0;
    for (de = 0; de < nAa; de++) {
      //sumScore += dynProgArray[ evidRow ][ aaMassCol[ de ] ];
      sumScore += dynProgArray[evidRow][aaMassCol[de]] * aaFreqC[de];
    }
    dynProgArray[row][col] = sumScore;
  }

  int colScoreCount = maxAaMass + colLast;
  scoreCount.resize(nRow);
  for (int row = 0; row < nRow; row++) {
    scoreCount[row] = dynProgArray[row][colScoreCount];
  }
  scoreOffset = initCountRow;

  // clean up
  for (int row = 0; row < nRow; row++) {
    delete [] dynProgArray[row];
  }
  delete [] dynProgArray;
  delete [] aaMassCol;
}

void TideSearchApplication::processParams() {
  const string index = Params::GetString("tide database");
  if (!FileUtils::Exists(index)) {
    carp(CARP_FATAL, "'%s' does not exist", index.c_str());
  } else if (FileUtils::IsRegularFile(index)) {
    // Index is FASTA file
    carp(CARP_INFO, "Creating index from '%s'", index.c_str());
    string targetIndexName = Params::GetString("store-index");
    if (targetIndexName.empty()) {
      targetIndexName = FileUtils::Join(Params::GetString("output-dir"),
                                        "tide-search.tempindex");
      remove_index_ = targetIndexName;
    }
    TideIndexApplication indexApp;
    indexApp.processParams();
    if (indexApp.main(index, targetIndexName) != 0) {
      carp(CARP_FATAL, "tide-index failed.");
    }
    Params::Set("tide database", targetIndexName);
  } else {
    // Index is Tide index directory
    pb::Header peptides_header;
    string peptides_file = FileUtils::Join(index, "pepix");
    HeadedRecordReader peptide_reader(peptides_file, &peptides_header);
    if ((peptides_header.file_type() != pb::Header::PEPTIDES) ||
        !peptides_header.has_peptides_header()) {
      carp(CARP_FATAL, "Error reading index (%s).", peptides_file.c_str());
    }

    const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();

    Params::Set("enzyme", pepHeader.enzyme());
    const char* digestString =
      digest_type_to_string(pepHeader.full_digestion() ? FULL_DIGEST : PARTIAL_DIGEST);
    Params::Set("digestion", digestString);
    Params::Set("isotopic-mass", pepHeader.monoisotopic_precursor() ? "mono" : "average");
  }
  // run param-medic?
  const string autoPrecursor = Params::GetString("auto-precursor-window");
  const string autoFragment = Params::GetString("auto-mz-bin-width");
  if (autoPrecursor != "false" || autoFragment != "false") {
    if (autoPrecursor != "false" && Params::GetString("precursor-window-type") != "ppm") {
      carp(CARP_FATAL, "Automatic peptide mass tolerance detection is only supported with ppm "
                       "units. Please re-run with auto-precursor-window set to 'false' or "
                       "precursor-window-type set to 'ppm'.");
    }
    ParamMedicErrorCalculator errCalc;
    errCalc.processFiles(Params::GetStrings("tide spectra file"));
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
        carp(CARP_INFO, "Precursor ppm standard deviation: %f", precursorSigmaPpm);
        carp(CARP_INFO, "Precursor error estimate (ppm): %.2f", precursorPredictionPpm);
        Params::Set("precursor-window", precursorPredictionPpm);
      } else {
        carp(autoPrecursor == "fail" ? CARP_FATAL : CARP_ERROR,
             "failed to calculate precursor error: %s", precursorFailure.c_str());
      }
    }
    if (autoFragment != "false") {
      if (fragmentFailure.empty()) {
        carp(CARP_INFO, "Fragment standard deviation (ppm): %f", fragmentSigmaPpm);
        carp(CARP_INFO, "Fragment bin size estimate (Th): %.4f", fragmentPredictionTh);
        Params::Set("mz-bin-width", fragmentPredictionTh);
      } else {
        carp(autoFragment == "fail" ? CARP_FATAL : CARP_ERROR,
             "failed to calculate fragment error: %s", fragmentFailure.c_str());
      }
    }
  }
}

void TideSearchApplication::setSpectrumFlag(map<pair<string, unsigned int>, bool>* spectrum_flag) {
  spectrum_flag_ = spectrum_flag;
}

string TideSearchApplication::getOutputFileName() {
  return output_file_name_;
}

//Added by Andy Lin in Feb 2016
//Determines the mass bin each peptide candidate (active_peptide_queue) is in
//pepMassInt will contain the a mass bin for each peptide candidate
//pepMassIntUnique will contain the unique set of mass bins
//pepMassInt and pepMassIntUnique are initalized right before call
//The length of pepMassInt is the number of canddiate peptides
void TideSearchApplication::getMassBin(
  vector<int>& pepMassInt,
  vector<int>& pepMassIntUnique,
  ActivePeptideQueue* active_peptide_queue,
  vector<bool>* candidatePeptideStatus
) {
  int pe = 0;
  int peidx = 0;

  deque<Peptide*>::const_iterator iter_ = active_peptide_queue->iter_;

  for (iter_ = active_peptide_queue->iter_; iter_ != active_peptide_queue->end_; ++iter_) {
    if ((*candidatePeptideStatus)[peidx]) {
      double pepMass = (*iter_)->Mass();
      int pepMaInt = MassConstants::mass2bin(pepMass);
      pepMassInt[pe] = pepMaInt;
      pepMassIntUnique.push_back(pepMaInt);
      pe++;
    }
    peidx++;
  }

  //For pepMassIntUnique vector
  //Sort vector, take unique of vector, get rid of extra space in vector
  std::sort(pepMassIntUnique.begin(), pepMassIntUnique.end());
  vector<int>::iterator last = std::unique(pepMassIntUnique.begin(),
                                           pepMassIntUnique.end());
  pepMassIntUnique.erase(last, pepMassIntUnique.end());
}

//Added by Andy Lin in March 2016
//Functions returns max value in curResidueEvidenceMatrix
//Function assumes that all values in curResidueEvidenceMatrix have been rounded to int
//Once function runs, maxColEvidence will contain the max evidence in
//each column of curResidueEvidenceMatrix
int TideSearchApplication::getMaxColEvidence(
  const vector<vector<double> >& curResidueEvidenceMatrix,
  vector<int>& maxColEvidence,
  int pepMassInt
) {
  assert(maxColEvidence.size() == curResidueEvidenceMatrix[0].size());

  int maxEvidence = -1;

  for (int curAA = 0; curAA < curResidueEvidenceMatrix.size(); curAA++) {
    for (int curMassBin = 0; curMassBin < pepMassInt; curMassBin++) {
      if (curResidueEvidenceMatrix[curAA][curMassBin] > maxColEvidence[curMassBin]) {
        maxColEvidence[curMassBin] = curResidueEvidenceMatrix[curAA][curMassBin];
      }
      if (curResidueEvidenceMatrix[curAA][curMassBin] > maxEvidence) {
        maxEvidence = curResidueEvidenceMatrix[curAA][curMassBin];
      }
    }
  }
  assert(maxEvidence >= 0);
  return maxEvidence;
}

//Added by Andy Lin in Nov 2016
//Calculates residue evidence score given a
//residue evidence matrix and a theoretical spectrum
int TideSearchApplication::calcResEvScore(
  const vector<vector<double> >& curResidueEvidenceMatrix,
  const vector<unsigned int>& intensArrayTheor,
  const vector<double>& aaMassDouble,
  Peptide* curPeptide
) {
  //Make sure the number of theoretical peaks match pepLen
  int pepLen = curPeptide->Len();
  assert(intensArrayTheor.size() == pepLen - 1);

  int scoreResidueEvidence = 0;
  double* residueMasses = curPeptide->getAAMasses(); //retrieves the amino acid masses, modifications included
  for (int res = 0; res < pepLen - 1; res++) {
    double tmpAAMass = residueMasses[res];
    int tmpAA = find(aaMassDouble.begin(),aaMassDouble.end(),tmpAAMass) - aaMassDouble.begin();
    scoreResidueEvidence += curResidueEvidenceMatrix[tmpAA][intensArrayTheor[res]-1];
  }
  delete residueMasses;
  return scoreResidueEvidence;
}

//Added by Andy Lin in Dec 2016
//Function takes a value, which results from multiplying various p-values together,
//and computes a new p-value from the distribution of correlated p-values
//Use eqn 3 from Tim Baily and Bill Noble Grundy RECOMB99 paper
double TideSearchApplication::calcCombinedPval(
  double m, //parameter
  double p, //value is the multiplication of p-values that will be combined,
  int numPval //number of p-values to combine
) {
  //compute useful quantities
  int intPartofM = int(m);
  double realPartofM = m - intPartofM;
  double y = m / numPval;
  double lnpy = -1.0 * log(pow(p, y));

  //compute first term (p-values are completely dependent)
  double firstTerm = 0.0;
  for (int i = 0; i < intPartofM; i++) {
    firstTerm += pow(lnpy, i) / double(factorial(i));
  }
  firstTerm = pow(p, y) * firstTerm;

  //compute second term (p-values are completely independent)
  double secondTerm = pow(p, y) * realPartofM * pow(lnpy, intPartofM) / (double)factorial(intPartofM);

  return firstTerm + secondTerm;
}

int TideSearchApplication::factorial(int n) {
  int product = 1;
  for (int i = 1; i <= n; i++) {
    product *= i;
  }
  return product;
}
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
