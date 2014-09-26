#include <cstdio>
#include "tide/abspath.h"
#include "tide/records_to_vector-inl.h"

#include "carp.h"
#include "parameter.h"
#include "SpectrumRecordWriter.h"
#include "TideIndexApplication.h"
#include "TideSearchApplication.h"

extern HASH_T* parameters;
extern AA_MOD_T* list_of_mods[MAX_AA_MODS];
extern int num_mods;

bool TideSearchApplication::HAS_DECOYS = false;

TideSearchApplication::TideSearchApplication() {
  exact_pval_search_ = false;
}

TideSearchApplication::~TideSearchApplication() {
}

int TideSearchApplication::main(int argc, char** argv) {

  const char* option_list[] = {
    "precursor-window",
    "precursor-window-type",
    "spectrum-min-mz",
    "spectrum-max-mz",
    "min-peaks",
    "spectrum-charge",
    "scan-number",
    "top-match",
    "store-spectra",
    "concat",
    "compute-sp",
    "remove-precursor-peak",
    "remove-precursor-tolerance",
    "print-search-progress",
    "txt-output",
    "sqt-output",
    "pepxml-output",
    "mzid-output",
    "pin-output",
    "fileroot",
    "output-dir",
    "overwrite",
    "parameter-file",
    "exact-p-value",
    "use-neutral-loss-peaks",
    "use-flanking-peaks",
    "mz-bin-width",
    "mz-bin-offset",
    "verbosity"
  };
  int num_options = sizeof(option_list) / sizeof(char*);
  const char* arg_list[] = {
    "tide spectra file",
    "tide database index"
  };
  int num_args = sizeof(arg_list) / sizeof(char*);
  initialize(arg_list, num_args, option_list, num_options, argc, argv);

  carp(CARP_INFO, "Running tide-search...");

  string cmd_line = "crux tide-search";
  for (int i = 1; i < argc; ++i) {
    cmd_line += " ";
    cmd_line += argv[i];
  }

  string index_dir = get_string_parameter_pointer("tide database index");
  string peptides_file = index_dir + "/pepix";
  string proteins_file = index_dir + "/protix";
  string auxlocs_file = index_dir + "/auxlocs";
  string spectra_file = get_string_parameter_pointer("tide spectra file");

  double window = get_double_parameter("precursor-window");
  WINDOW_TYPE_T window_type = get_window_type_parameter("precursor-window-type");

  // Check spectrum-charge parameter
  string charge_string = get_string_parameter_pointer("spectrum-charge");
  int charge_to_search;
  if (charge_string == "all") {
    carp(CARP_DEBUG, "Searching all charge states");
    charge_to_search = 0;
  } else {
    charge_to_search = atoi(charge_string.c_str());
    if (charge_to_search < 1 || charge_to_search > 6) {
      carp(CARP_FATAL, "Invalid spectrum-charge value %s",
           charge_string.c_str());
    }
    carp(CARP_DEBUG, "Searching charge state %d", charge_to_search);
  }

  // Check scan-number parameter
  string scan_range = get_string_parameter_pointer("scan-number");
  int min_scan, max_scan;
  if (scan_range == "__NULL_STR") {
    min_scan = 0;
    max_scan = BILLION;
    carp(CARP_DEBUG, "Searching all scans");
  }
  else if (scan_range.find('-') == string::npos) {
    // Single scan
    min_scan = max_scan = atoi(scan_range.c_str());
    carp(CARP_DEBUG, "Searching single scan %d", min_scan);
  } else {
    if (!get_range_from_string(scan_range.c_str(), min_scan, max_scan)) {
      carp(CARP_FATAL, "The scan number range '%s' is invalid. "
           "Must be of the form <first>-<last>", scan_range.c_str());
    } else {
      if (min_scan > max_scan) {
        int tmp_scan = min_scan;
        min_scan = max_scan;
        max_scan = tmp_scan;
        carp(CARP_DEBUG, "Switched scan range min and max");
      }
      carp(CARP_DEBUG, "Searching scan range %d-%d", min_scan, max_scan);
    }
  }
  //check to compute exact p-value 
  exact_pval_search_ = get_boolean_parameter("exact-p-value");
  bin_width_  = get_double_parameter("mz-bin-width");
  bin_offset_ = get_double_parameter("mz-bin-offset");
  // for now don't allow XCorr p-value searches with variable bin width
  if (exact_pval_search_ && abs(bin_width_ - BIN_WIDTH_MONO) > 0.000001) {
    carp(CARP_FATAL, "tide-search with XCorr p-values and variable bin width "
                     "is not allowed in this version of Crux.");
  }

  // Check concat parameter
  bool concat = get_boolean_parameter("concat");
  if (concat) {
    OutputFiles::setConcat();
  }

  // Check compute-sp parameter
  bool compute_sp = get_boolean_parameter("compute-sp");
  if (get_boolean_parameter("sqt-output") && !compute_sp){
    compute_sp = true;
    carp(CARP_INFO, "Enabling parameter compute-sp since SQT output is enabled "
                    " (this will increase runtime).");
  }

  const double max_mz = 0.0;

  carp(CARP_INFO, "Reading index %s", index_dir.c_str());
  // Read proteins index file
  ProteinVec proteins;
  ActivePeptideQueue* active_peptide_queue = NULL;
  pb::Header protein_header;
  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins,
      proteins_file, &protein_header)) {
    carp(CARP_FATAL, "Error reading index (%s)", proteins_file.c_str());
  }
  carp(CARP_DEBUG, "Read %d proteins", proteins.size());

  //open a copy of peptide buffer for Amino Acid Frequency (AAF) calculation.
  double* aaFreqN = NULL;
  double* aaFreqI = NULL;
  double* aaFreqC = NULL;
  int* aaMass = NULL;
  int nAA;
  
  if (exact_pval_search_) {
    pb::Header aaf_peptides_header;
    HeadedRecordReader aaf_peptide_reader(peptides_file, &aaf_peptides_header);

    if (!aaf_peptides_header.file_type() == pb::Header::PEPTIDES ||
        !aaf_peptides_header.has_peptides_header()) {
      carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
    }
    MassConstants::Init(&aaf_peptides_header.peptides_header().mods(),
                        bin_width_, bin_offset_);
    active_peptide_queue = new ActivePeptideQueue(aaf_peptide_reader.Reader(), proteins);
    nAA = active_peptide_queue->CountAAFrequency(bin_width_, bin_offset_,
                                                 &aaFreqN, &aaFreqI, &aaFreqC, &aaMass);
    delete active_peptide_queue;
  } // End calculation AA frequencies

  // Read auxlocs index file
  vector<const pb::AuxLocation*> locations;
  if (!ReadRecordsToVector<pb::AuxLocation>(&locations, auxlocs_file)) {
    carp(CARP_FATAL, "Error reading index (%s)", auxlocs_file.c_str());
  }
  carp(CARP_DEBUG, "Read %d auxlocs", locations.size());

  // Read peptides index file
  pb::Header peptides_header;
  HeadedRecordReader peptide_reader(peptides_file, &peptides_header);
  if (!peptides_header.file_type() == pb::Header::PEPTIDES ||
      !peptides_header.has_peptides_header()) {
    carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
  }
  
  const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();
  DECOY_TYPE_T headerDecoyType = (DECOY_TYPE_T)pepHeader.decoys();
  if (headerDecoyType != NO_DECOYS) {
    HAS_DECOYS = true;
    if (headerDecoyType == PROTEIN_REVERSE_DECOYS) {
      OutputFiles::setProteinLevelDecoys();
    }
  }

  MassConstants::Init(&pepHeader.mods(), bin_width_, bin_offset_);
  TideMatchSet::initModMap(pepHeader.mods());

  active_peptide_queue = new ActivePeptideQueue(peptide_reader.Reader(), proteins);

  active_peptide_queue->SetBinSize(bin_width_, bin_offset_);

  carp(CARP_INFO, "Reading spectra file %s", spectra_file.c_str());
  // Try to read file as spectrumrecords file
  SpectrumCollection spectra;
  pb::Header spectrum_header;
  string delete_spectra_file = "";
  if (!spectra.ReadSpectrumRecords(spectra_file, &spectrum_header)) {
    // Failed, try converting to spectrumrecords file
    carp(CARP_INFO, "Converting %s to spectrumrecords format",
                    spectra_file.c_str());
    string converted_spectra_file = get_string_parameter_pointer("store-spectra");
    if (converted_spectra_file.empty()) {
      delete_spectra_file = converted_spectra_file =
      make_file_path("spectrumrecords.tmp");
    }
    carp(CARP_DEBUG, "New spectrumrecords filename: %s",
                     converted_spectra_file.c_str());
    if (!SpectrumRecordWriter::convert(spectra_file, converted_spectra_file)) {
      carp(CARP_FATAL, "Error converting %s to spectrumrecords format",
                       spectra_file.c_str());
    }
    carp(CARP_DEBUG, "Reading converted spectra file %s",
                     spectra_file.c_str());
    // Re-read converted file as spectrumrecords file
    if (!spectra.ReadSpectrumRecords(converted_spectra_file, &spectrum_header)) {
      carp(CARP_DEBUG, "Deleting %s", converted_spectra_file.c_str());
      remove(converted_spectra_file.c_str());
      carp(CARP_FATAL, "Error reading spectra file %s",
                       converted_spectra_file.c_str());
    }
  } else {
    // Successfully read file as spectrumrecords format
    if (get_int_parameter("remove_precursor_peak") != 0) {
      carp(CARP_FATAL, "remove_precursor_peak can only be used during conversion "
                       "to spectrumrecords format.");
    }
  }

  carp(CARP_INFO, "Sorting spectra");
  if (window_type != WINDOW_MZ) {
    spectra.Sort();
  } else {
    spectra.Sort<ScSortByMz>(ScSortByMz(window));
  }
  if (max_mz == 0) {
    double highest_mz = spectra.FindHighestMZ();
    unsigned int spectrum_num = spectra.SpecCharges()->size();
    if (spectrum_num > 0 && exact_pval_search_) {
      highest_mz = spectra.SpecCharges()->at(spectrum_num - 1).neutral_mass;
    }

    carp(CARP_DEBUG, "Max m/z %f", highest_mz);
    MaxBin::SetGlobalMax(highest_mz);
  } else {
    MaxBin::SetGlobalMaxFromFlag();
  }

  char* digestString =
    digest_type_to_string(pepHeader.full_digestion() ? FULL_DIGEST : PARTIAL_DIGEST);

  bool txt_only = !get_boolean_parameter("sqt-output") &&
                  !get_boolean_parameter("pepxml-output") &&
                  !get_boolean_parameter("mzid-output") &&
                  !get_boolean_parameter("pin-output");
  OutputFiles* output_files = NULL;
  ofstream* target_file = NULL;
  ofstream* decoy_file = NULL;
  if (!txt_only) {
    carp(CARP_DEBUG, "Using OutputFiles to write matches");
    // Overwrite enzyme/digestion parameters in the hash
    // TODO Find a better way to do this?
    add_or_update_hash(parameters, "enzyme", pepHeader.enzyme().c_str());
    add_or_update_hash(parameters, "digestion", digestString);
    add_or_update_hash(parameters, "monoisotopic-precursor", pepHeader.monoisotopic_precursor()?"T":"F");
    free(digestString);
    output_files = new OutputFiles(this);
  } else {
    carp(CARP_DEBUG, "Using TideMatchSet to write matches");
    bool overwrite = get_boolean_parameter("overwrite");
    stringstream ss;
    ss << pepHeader.enzyme() << '-' << digestString;
    free(digestString);
    TideMatchSet::setCleavageType(ss.str());
    if (!concat) {
      string target_file_name = make_file_path("tide-search.target.txt");
      target_file = create_stream_in_path(target_file_name.c_str(), NULL, overwrite);
      if (HAS_DECOYS) {
        string decoy_file_name = make_file_path("tide-search.decoy.txt");
        decoy_file = create_stream_in_path(decoy_file_name.c_str(), NULL, overwrite);
      }
    } else {
      string concat_file_name = make_file_path("tide-search.txt");
      target_file = create_stream_in_path(concat_file_name.c_str(), NULL, overwrite);
    }
  }

  // Do the search
  carp(CARP_INFO, "Running search");
  cleanMods();
  search(spectra.SpecCharges(), active_peptide_queue, proteins, locations, window,
         window_type, get_double_parameter("spectrum-min-mz"),
         get_double_parameter("spectrum-max-mz"), min_scan, max_scan,
         get_int_parameter("min-peaks"), charge_to_search,
         get_int_parameter("top-match"), spectra.FindHighestMZ(),
         output_files, target_file, decoy_file, compute_sp,
         nAA, aaFreqN, aaFreqI, aaFreqC, aaMass);
  // Delete temporary spectrumrecords file
  if (!delete_spectra_file.empty()) {
    carp(CARP_DEBUG, "Deleting %s", delete_spectra_file.c_str());
    remove(delete_spectra_file.c_str());
  }

  // Clean up
  delete active_peptide_queue;
  for (ProteinVec::const_iterator i = proteins.begin(); i != proteins.end(); ++i) {
    delete const_cast<pb::Protein*>(*i);
  }
  if (output_files) {
    delete output_files;
  }
  if (target_file) {
    target_file->close();
    delete target_file;
    if (decoy_file) {
      decoy_file->close();
      delete decoy_file;
    }
  }
  delete aaFreqN;
  delete aaFreqI;
  delete aaFreqC;
  delete aaMass;

  return 0;
}

/**
 * Free all existing mods
 */
void TideSearchApplication::cleanMods() {
  for (int i = 0; i < MAX_AA_MODS; ++i) {
    free_aa_mod(list_of_mods[i]);
    list_of_mods[i] = NULL;
  }
  num_mods = 0;
}

void TideSearchApplication::search(
  const vector<SpectrumCollection::SpecCharge>* spec_charges,
  ActivePeptideQueue* active_peptide_queue,
  const ProteinVec& proteins,
  const vector<const pb::AuxLocation*>& locations,
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
  OutputFiles* output_files,
  ofstream* target_file,
  ofstream* decoy_file,
  bool compute_sp,
  int nAA, 
  double* aaFreqN,
  double* aaFreqI,
  double* aaFreqC,
  int* aaMass
) {
  if (output_files) {
    output_files->exact_pval_search_ = exact_pval_search_;
    output_files->writeHeaders();
  } else if (target_file) {
    TideMatchSet::writeHeaders(target_file, false, compute_sp, exact_pval_search_);
    TideMatchSet::writeHeaders(decoy_file, true, compute_sp, exact_pval_search_);
  }

  // This is the main search loop.
  ObservedPeakSet observed(bin_width_, bin_offset_,
  get_boolean_parameter("use-neutral-loss-peaks"), 
  get_boolean_parameter("use-flanking-peaks"));

  // cycle through spectrum-charge pairs, sorted by neutral mass
  unsigned sc_index = 0;
  FLOAT_T sc_total = (FLOAT_T)spec_charges->size();
  int print_interval = get_int_parameter("print-search-progress");
  for (vector<SpectrumCollection::SpecCharge>::const_iterator sc = spec_charges->begin();
       sc != spec_charges->end();
       ++sc) {
    const Spectrum* spectrum = sc->spectrum;

    double precursor_mz = spectrum->PrecursorMZ();
    int charge = sc->charge;
    int scan_num = spectrum->SpectrumNumber();
    if (precursor_mz < spectrum_min_mz || precursor_mz > spectrum_max_mz ||
        scan_num < min_scan || scan_num > max_scan ||
        spectrum->Size() < min_peaks ||
        (search_charge != 0 && charge != search_charge)) {
      continue;
    }

    // The active peptide queue holds the candidate peptides for spectrum.
    // Calculate and set the window, depending on the window type.
    double min_mass, max_mass;
    computeWindow(*sc, window_type, precursor_window, &min_mass, &max_mass);
    if (!exact_pval_search_) {  //execute original tide-search program

      // Normalize the observed spectrum and compute the cache of
      // frequently-needed values for taking dot products with theoretical
      // spectra.
      observed.PreprocessSpectrum(*spectrum, charge);
      int nCandPeptide = active_peptide_queue->SetActiveRange(min_mass, max_mass);
      TideMatchSet::Arr2 match_arr2(nCandPeptide); // Scored peptides will go here.

      // Programs for taking the dot-product with the observed spectrum are laid
      // out in memory managed by the active_peptide_queue, one program for each
      // candidate peptide. The programs will store the results directly into
      // match_arr. We now pass control to those programs.
      collectScoresCompiled(active_peptide_queue, spectrum, observed, &match_arr2,
                            nCandPeptide, charge);

      // matches will arrange the results in a heap by score, return the top
      // few, and recover the association between counter and peptide. We output
      // the top matches.
      TideMatchSet::Arr match_arr(nCandPeptide);
      for (TideMatchSet::Arr2::iterator it = match_arr2.begin();
           it != match_arr2.end();
           ++it) {
        TideMatchSet::Pair pair;
        pair.first.first = (double)(it->first / 100000000.0);
        pair.first.second = 0.0;
        pair.second = it->second;
        match_arr.push_back(pair);
      }

      TideMatchSet matches(&match_arr, highest_mz);
      matches.exact_pval_search_ = exact_pval_search_;
      if (output_files) {
        matches.report(output_files, top_matches, spectrum, charge,
                       active_peptide_queue, proteins, locations, compute_sp, true);
      } else {
        matches.report(target_file, decoy_file, top_matches, spectrum, charge,
                       active_peptide_queue, proteins, locations, compute_sp, true);
      }
    } else {  // execute exact-pval-search
      const int minDeltaMass = aaMass[0];
      const int maxDeltaMass = aaMass[nAA - 1];

      int maxPrecurMass = floor(MaxBin::Global().CacheBinEnd() + 50.0); // TODO works, but is this the best way to get?
      int nCandPeptide = active_peptide_queue->SetActiveRangeBIons(min_mass, max_mass);
      TideMatchSet::Arr match_arr(nCandPeptide); // scored peptides will go here.
  
      // iterators needed at multiple places in following code
      deque<Peptide*>::const_iterator iter_ = active_peptide_queue->iter_;
      deque<TheoreticalPeakSetBIons>::const_iterator iter1_ = active_peptide_queue->iter1_;
      vector<int>::const_iterator iter_int;
      vector<unsigned int>::const_iterator iter_uint;

      //************************************************************************
      /* For one observed spectrum, calculates:
       *  - vector of cleavage evidence
       *  - score count vectors for a range of integer masses
       *  - p-values of XCorr match scores between spectrum and all selected candidate target and decoy peptides
       * Written by Jeff Howbert, October, 2013.
       * Ported to and integrated with Tide by Jeff Howbert, November, 2013.
       */
      int pe;
      int ma;
      int pepMaInt;
      int* pepMassInt = new int[nCandPeptide];
      vector<int> pepMassIntUnique;
      pepMassIntUnique.reserve(nCandPeptide);
      pe = 0;
      for (iter_ = active_peptide_queue->iter_;
           iter_ != active_peptide_queue->end_;
           ++iter_) {
        double pepMass = (*iter_)->Mass();
        pepMaInt = (int)floor(pepMass / bin_width_ + 1.0 - bin_offset_);
        pepMassInt[pe] = pepMaInt;
        pepMassIntUnique.push_back(pepMaInt);
        pe++;
      }
      std::sort(pepMassIntUnique.begin(), pepMassIntUnique.end());
      vector<int>::iterator last = std::unique(pepMassIntUnique.begin(),
                                               pepMassIntUnique.end());
      pepMassIntUnique.erase(last, pepMassIntUnique.end());
      int nPepMassIntUniq = (int)pepMassIntUnique.size();

      int** evidenceObs = new int*[nPepMassIntUniq];
      int* scoreOffsetObs = new int[nPepMassIntUniq];
      double** pValueScoreObs = new double*[nPepMassIntUniq];
      int* intensArrayTheor = new int [maxPrecurMass]; // initialized later in loop
      for (pe = 0; pe < nPepMassIntUniq; pe++) { // TODO should probably instead use iterator over pepMassIntUnique
        evidenceObs[pe] = new int[maxPrecurMass];
        for (ma = 0; ma < maxPrecurMass; ma++) {
          evidenceObs[pe][ma] = 0;
        }
        scoreOffsetObs[pe] = 0;
        pepMaInt = pepMassIntUnique[pe]; // TODO should be accessed with an iterator
        // preprocess to create one integerized evidence vector for each cluster of masses among selected peptides
        double pepMassMonoMean = (pepMaInt - 1.0 + bin_offset_) * bin_width_ + 0.5;
        observed.CreateEvidenceVector(*spectrum, bin_width_, bin_offset_, charge,
                                      pepMassMonoMean, maxPrecurMass, evidenceObs[pe]);
        // NOTE: will have to go back to separate dynamic programming for
        //       target and decoy if they have different probNI and probC
        int maxEvidence = *std::max_element(evidenceObs[pe], evidenceObs[pe] + maxPrecurMass);
        int minEvidence = *std::min_element(evidenceObs[pe], evidenceObs[pe] + maxPrecurMass);
        // estimate maxScore and minScore
        int maxNResidue = (int)floor((double)pepMaInt / (double)minDeltaMass);
        vector<int> sortEvidenceObs (evidenceObs[pe], evidenceObs[pe] + maxPrecurMass);
        std::sort(sortEvidenceObs.begin(), sortEvidenceObs.end(), greater<int>());
        int maxScore = 0;
        int minScore = 0;
        for (int sc = 0; sc < maxNResidue; sc++) {
          maxScore += sortEvidenceObs[sc];
        }
        for (int sc = maxPrecurMass - maxNResidue; sc < maxPrecurMass; sc++) {
          minScore += sortEvidenceObs[sc];
        }
        int bottomRowBuffer = maxEvidence + 1;
        int topRowBuffer = -minEvidence;
        int nRowDynProg = bottomRowBuffer - minScore + 1 + maxScore + topRowBuffer;
        pValueScoreObs[pe] = new double[nRowDynProg];
        scoreOffsetObs[pe] = calcScoreCount(maxPrecurMass, evidenceObs[pe], pepMaInt,
                                            maxEvidence, minEvidence, maxScore, minScore, 
                                            nAA, aaFreqN, aaFreqI, aaFreqC, aaMass,
                                            pValueScoreObs[pe]);
      }

      // ***** calculate p-values for peptide-spectrum matches ***********************************
      iter_ = active_peptide_queue -> iter_;
      iter1_ = active_peptide_queue -> iter1_;
      for (pe = 0; pe < nCandPeptide; pe++) { // TODO should probably use iterator instead
        int pepMassIntIdx = 0;
        for (ma = 0; ma < nPepMassIntUniq; ma++) { // TODO should probably use iterator instead
          if (pepMassIntUnique[ma] == pepMassInt[pe]) { // TODO pepMassIntUnique should be accessed with an iterator
            pepMassIntIdx = ma;
            break;
          }
        }
        // score XCorr for target peptide with integerized evidenceObs array
        for (ma = 0; ma < maxPrecurMass; ma++) {
          intensArrayTheor[ma] = 0;
        }
        for (iter_uint = iter1_->unordered_peak_list_.begin();
             iter_uint != iter1_->unordered_peak_list_.end();
             iter_uint++) {
          intensArrayTheor[*iter_uint] = 1;
        }

        int scoreRefactInt = 0;
        for (ma = 0; ma < maxPrecurMass; ma++) {
          scoreRefactInt += evidenceObs[pepMassIntIdx][ma] * intensArrayTheor[ma];
        }
        int scoreCountIdx = scoreRefactInt + scoreOffsetObs[pepMassIntIdx];
        double pValue = pValueScoreObs[pepMassIntIdx][scoreCountIdx];

        TideMatchSet::Pair pair;
        pair.first.first = pValue;
        pair.first.second = (double)scoreRefactInt;
        pair.second = nCandPeptide - pe; // TODO ugly hack to conform with the way these indices are generated in standard tide-search
        match_arr.push_back(pair);

        // move to next peptide and b ion queue
        ++iter_; // TODO need to add test to make sure haven't gone past available peptides
        ++iter1_; // TODO need to add test to make sure haven't gone past available b ion queues
      }

      // clean up
      delete [] pepMassInt;
      delete [] scoreOffsetObs;
      for (pe = 0; pe < nPepMassIntUniq; pe++) {
        delete [] evidenceObs[pe];
        delete [] pValueScoreObs[pe];
      }
      delete [] evidenceObs;
      delete [] pValueScoreObs;
      delete [] intensArrayTheor;

      // matches will arrange the results in a heap by score, return the top
      // few, and recover the association between counter and peptide. We output
      // the top matches.
      TideMatchSet matches(&match_arr, highest_mz);
      matches.exact_pval_search_ = exact_pval_search_;
      if (output_files) {
        matches.report(output_files, top_matches, spectrum, charge,
                     active_peptide_queue, proteins, locations, compute_sp, false);
      } else {
        matches.report(target_file, decoy_file, top_matches, spectrum, charge,
                     active_peptide_queue, proteins, locations, compute_sp, false);
      }
    }

    ++sc_index;
    if (print_interval > 0 && sc_index % print_interval == 0) {
      carp(CARP_INFO, "%d spectra searched, %.0f%% complete",
           sc_index, sc_index / sc_total * 100);
    }
  }
  if (output_files) {
    output_files->writeFooters();
  }
}

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

void TideSearchApplication::computeWindow(
  const SpectrumCollection::SpecCharge& sc,
  WINDOW_TYPE_T window_type,
  double precursor_window,
  double* out_min,
  double* out_max
) {
  switch (window_type) {
  case WINDOW_MASS:
    *out_min = sc.neutral_mass - precursor_window;
    *out_max = sc.neutral_mass + precursor_window;
    break;
  case WINDOW_MZ: {
    double mz_minus_proton = sc.spectrum->PrecursorMZ() - MASS_PROTON;
    *out_min = (mz_minus_proton - precursor_window) * sc.charge;
    *out_max = (mz_minus_proton + precursor_window) * sc.charge;
    break;
    }
  case WINDOW_PPM: {
    double tiny_precursor = precursor_window * 1e-6;
    *out_min = sc.neutral_mass * (1.0 - tiny_precursor);
    *out_max = sc.neutral_mass * (1.0 + tiny_precursor);
    break;
    }
  default:
    carp(CARP_FATAL, "Invalid window type");
  }
  carp(CARP_DETAILED_DEBUG, "Scan %d.%d mass window is [%f, %f]",
       sc.spectrum->SpectrumNumber(), sc.charge, *out_min, *out_max);
}

bool TideSearchApplication::hasDecoys() {
  return HAS_DECOYS;
}

string TideSearchApplication::getName() {
  return "tide-search";
}

string TideSearchApplication::getDescription() {
  return
  "Search a collection of spectra against a sequence database, returning a "
  "collection of peptide-spectrum matches (PSMs). This is a fast search engine "
  "but requires that you first build an index with tide-index.";
}

bool TideSearchApplication::needsOutputDirectory() {
  return true;
}

COMMAND_T TideSearchApplication::getCommand() {
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
  int colStart = (int)floor(1.0 / bin_width_ + 1.0 - bin_offset_);
  int scoreOffsetObs = bottomRowBuffer - minScore;

  int nRow = bottomRowBuffer - minScore + 1 + maxScore + topRowBuffer;
  int nCol = colBuffer + pepMassInt;
  int rowFirst = bottomRowBuffer;
  int rowLast = rowFirst - minScore + maxScore;
  int colFirst = colStart + 1;
  int colLast = pepMassInt - 17;
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
  scoreCountBinAdjust = new double [nRow];
  for (row = 0; row < nRow; row++) {
    scoreCountBinAdjust[row] = 0.0;
  }

  dynProgArray[initCountRow][initCountCol] = 1.0; // initial count of peptides with mass = 1
  int deltaMassCol[nDeltaMass];
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
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
