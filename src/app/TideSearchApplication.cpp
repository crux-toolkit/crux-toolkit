#include <cstddef>
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
#include "tide/spectrum_collection.h"
#include "util/Params.h"
#include "util/FileUtils.h"
#include "util/StringUtils.h"
#include <math.h> 
#include <map>
#include "tide/ActivePeptideQueue.h"
#include "residue_stats.pb.h"
#include "crux_version.h"

//here :O
//#include <iostream>
//#include <set>
/*
#include <vector>
#include <queue>
#include <sstream>
#include <cstdlib>
*/
/* This constant is the product of the original "magic number" (10000,
 * on line 4622 of search28.c) that was used to rescale the XCorr
 * score, and the integerization constant used by Benjamin Diament in
 * Tide.  In the Tide publication, that constant was reported as 10^7.
 * --WSN, 10 March 2015 */
const double TideSearchApplication::XCORR_SCALING = 100000000.0;

/* This constant is used to put the refactored XCorr back into the
 * same range as the original XCorr score.  It is the XCorr "magic
 * number" (10000) divided by the EVIDENCE_SCALE_INT (defined in
 * tide/spectrum_preprocess2.cc). */
const double TideSearchApplication::RESCALE_FACTOR = 20.0;

// Things done:
// -- handle terminal mod structure in the peptides internally. GH issue: #639
// -- missing charge states, 0 charege states,, override charge states handled. --> need to update doc. GH issue: #557, #607

#define CHECK(x) GOOGLE_CHECK(x)

TideSearchApplication::TideSearchApplication() {
  remove_index_ = "";
  spectrum_flag_ = NULL;
  decoy_num_ = 0;
  num_range_skipped_ = 0;
  num_precursors_skipped_ = 0;
  num_isotopes_skipped_ = 0;
  num_retained_ = 0;
  num_spectra_ = 0;
  num_spectra_searched_ = 0;
  total_candidate_peptides_ = 0;
  precursor_window_ = 0;
  spectrum_min_mz_ = 0; 
  spectrum_max_mz_ = 0;
  min_scan_ = 0;
  max_scan_ = 0;
  min_peaks_ = 0;
  min_precursor_charge_ = 0;
  max_precursor_charge_ = 0;
  out_tsv_target_  = NULL; // original tide-search output format in tab-delimited text files (txt)
  out_tsv_decoy_ = NULL;  // original tide-search output format in tab-delimited text files (txt) for the decoy psms only
  out_mztab_target_ = NULL;      // mzTAB output format
  out_mztab_decoy_ = NULL;      // mzTAB output format for the decoy psms only
  out_pin_target_ = NULL;        // pin output format for percolator
  out_pin_decoy_ = NULL;        // pin output format for percolator for the decoy psms only
  total_spectra_num_ = 0;       // The total number of spectra searched. This is counted during the spectrum conversion

  for (int i = 0; i < NUMBER_LOCK_TYPES; i++) {  // LOCK_TYPES are defined in model/objects.h
    locks_array_.push_back(new boost::mutex());
  }
}

TideSearchApplication::~TideSearchApplication() {
  for (int i = 0; i < NUMBER_LOCK_TYPES; i++) {
    delete locks_array_[i];
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

  print_interval_ = Params::GetInt("print-search-progress");  
  
  bin_width_  = Params::GetDouble("mz-bin-width");
  bin_offset_ = Params::GetDouble("mz-bin-offset");

  use_neutral_loss_peaks_ = Params::GetBool("use-neutral-loss-peaks");
  use_flanking_peaks_ = Params::GetBool("use-flanking-peaks");

  negative_isotope_errors_ = getNegativeIsotopeErrors();  

  window_type_= string_to_window_type(Params::GetString("precursor-window-type"));
  precursor_window_ = Params::GetDouble("precursor-window");

  spectrum_min_mz_ =  Params::GetDouble("spectrum-min-mz") ;
  spectrum_max_mz_ = Params::GetDouble("spectrum-max-mz") ;
  min_peaks_ = Params::GetInt("min-peaks");
  min_precursor_charge_ = Params::GetInt("min-precursor-charge");
  max_precursor_charge_ = Params::GetInt("max-precursor-charge");  

  fragTol_ = Params::GetDouble("fragment-tolerance");
  granularityScale_ = Params::GetInt("evidence-granularity");  

  top_matches_ = Params::GetInt("top-match");  
  num_threads_ = Params::GetInt("num-threads");
  if (num_threads_ < 1) {
    num_threads_ = boost::thread::hardware_concurrency(); // MINIMUM # = 1.
    // (Meaning just main thread) Do not make this value below 1.
  } else if (num_threads_ > 64) {
    // make sure that number of threads are reasonable, e.g. user did not specify millions of threads...
    carp(CARP_FATAL, "Requested more than 64 threads.");
  }
  carp(CARP_INFO, "Number of Threads: %d", num_threads_);


  // Check scan-number parameter
  string scan_range = Params::GetString("scan-number");
  if (scan_range.empty()) {
    min_scan_ = 0;
    max_scan_ = BILLION;
    carp(CARP_DEBUG, "Searching all scans");
  } else if (scan_range.find('-') == string::npos) {
    // Single scan
    min_scan_ = max_scan_ = atoi(scan_range.c_str());
    carp(CARP_INFO, "Searching single scan %d", min_scan_);
  } else {
    if (!get_range_from_string(scan_range.c_str(), min_scan_, max_scan_)) {
      carp(CARP_FATAL, "The scan number range '%s' is invalid. "
           "Must be of the form <first>-<last>.", scan_range.c_str());
    }
    if (min_scan_ > max_scan_) {
      carp(CARP_FATAL, "Invalid scan range: %d to %d.", min_scan_, max_scan_);
    }
    carp(CARP_INFO, "Searching scan range %d to %d.", min_scan_, max_scan_);
  }

  // Determine which score function to use for scoring PSMs and store in SCORE_FUNCTION enum. Change this in ./src/util/crux-utils.cpp: 68
  curScoreFunction_ = string_to_score_function_type(Params::GetString("score-function"));  

  if (curScoreFunction_ == PVALUES && bin_width_ < 1.0) {
    carp(CARP_FATAL, "Tide-search does not support P-value calculation with bin-width less than 1.0 Da.");
  }
  if (curScoreFunction_ >= NUMBER_SCORE_FUNCTIONS) {
    carp(CARP_FATAL, "Invalid score function.");
  }
 
  // Get a peptide reader to the peptide index datasets along with proteins, auxlocs. 
  vector<const pb::Protein*> proteins;
  vector<const pb::AuxLocation*> locations;
  pb::Header peptides_header;
  string peptides_file = FileUtils::Join(input_index, "pepix");  
  HeadedRecordReader peptide_reader = HeadedRecordReader(peptides_file, &peptides_header);
  getPeptideIndexData(input_index, proteins, locations, peptides_header);
  tide_index_mzTab_file_path_ = FileUtils::Join(input_index, TideIndexApplication::tide_index_mzTab_filename_);

  TideMatchSet::curScoreFunction_ = curScoreFunction_;
  TideMatchSet::top_matches_ = top_matches_;
  TideMatchSet::decoy_num_ = decoy_num_;
  TideMatchSet::mass_precision_ =  Params::GetInt("mass-precision");
  TideMatchSet::score_precision_ = Params::GetInt("precision");
  TideMatchSet::mod_precision_ = Params::GetInt("mod-precision");
  TideMatchSet::concat_ = Params::GetBool("concat");  

  // Create the output files, print headers
  createOutputFiles(); 

  // Convert the original file names into spectrum records if needed 
  // Update the file names in the variable inputFiles_ locally.
  // Run spectrum file convertion in parallel.
  for (vector<string>::const_iterator original_file_name = input_files.begin(); original_file_name != input_files.end(); ++original_file_name) {
    inputFiles_.push_back(TideSearchApplication::InputFile(*original_file_name, *original_file_name, false));
  }
  // Launch threads to convert files
  boost::thread_group threadgroup_input_files;
  for (int t = 1; t < num_threads_; ++t) {
    boost::thread * currthread = new boost::thread(boost::bind(&TideSearchApplication::getInputFiles, this, t));
    threadgroup_input_files.add_thread(currthread);
  }
  getInputFiles(0);
  // Join threads
  threadgroup_input_files.join_all();

  if (total_spectra_num_ > 0) {
    carp(CARP_INFO, "There were a total of %d spectrum conversions from %d input spectrum files.",
         total_spectra_num_, inputFiles_.size());
  }
  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);

  // Create the active_peptide_queues and peptide_readers for each threads
  vector<HeadedRecordReader*> peptide_reader_threads;
  vector<ActivePeptideQueue*> APQ;
  for (int i = 0; i < num_threads_; i++) {
    peptide_reader_threads.push_back(new HeadedRecordReader(peptides_file, &peptides_header));
    APQ.push_back(new ActivePeptideQueue(peptide_reader_threads.back()->Reader(), proteins, &locations));
  }

  carp(CARP_INFO, "Starting search.");
  // Read the first spectrum records from each input files 
  int file_cnt = 0;
  pb::Spectrum pb_spectrum; // spectrum message struct
  for (vector<InputFile>::iterator spectrum_file = inputFiles_.begin(); spectrum_file != inputFiles_.end(); ++spectrum_file, ++file_cnt) {

    string spectrum_records_file = spectrum_file->SpectrumRecords;
    spectrum_reader_.push_back(new HeadedRecordReader(spectrum_records_file));

    if (!spectrum_reader_.back()->Done() ) {
      spectrum_reader_.back()->Read(&pb_spectrum);
      spectrum_heap_.push_back(make_pair(pb_spectrum, file_cnt));
    }
    if ( !spectrum_reader_.back()->OK() ){
      carp(CARP_FATAL, "Spectrum records file %s is corrupt.", spectrum_records_file.c_str());
    }

  }
  // make a heap with the spectrum records 
  make_heap(spectrum_heap_.begin(), spectrum_heap_.end(), compare_spectrum());

  // Create thread data
  vector<thread_data> thread_data_array;
  for (int t = 0; t < num_threads_; ++t) {
      thread_data_array.push_back(thread_data(APQ[t], t));
  }

  // Launch threads
  boost::thread_group threadgroup;
  for (int t = 1; t < num_threads_; ++t) {
    boost::thread * currthread = new boost::thread(boost::bind(&TideSearchApplication::spectrum_search, this, (void *) &(thread_data_array[t])));
    threadgroup.add_thread(currthread);
  }

  // Searches through part of the spec charge vector while waiting for threads are busy
  spectrum_search( (void *)  &(thread_data_array[0]) );

  // Join threads
  threadgroup.join_all();

  // Print statistics
  long int total_peaks = num_precursors_skipped_ + num_isotopes_skipped_ + num_range_skipped_ + num_retained_;
  if (total_peaks == 0) {
    carp(CARP_INFO, "Warning: no peaks found.");
  } else {
    carp(CARP_INFO,
          "Deleted %d precursor, %d isotope and %d out-of-range peaks.",
          num_precursors_skipped_, num_isotopes_skipped_, num_range_skipped_);
  }
  if (num_retained_ == 0) {
    carp(CARP_INFO, "Warning: no peaks retained.");
  } else {
    carp(CARP_INFO, "Retained %g%% of peaks.", (100.0 * num_retained_) / total_peaks);
  }
  
  carp(CARP_INFO, "Time per spectrum-charge combination: %lf s.", wall_clock() / (1e6*num_spectra_searched_));
  carp(CARP_INFO, "Average number of candidates per spectrum-charge combination: %lf ",
    ((double)total_candidate_peptides_) /  (double)num_spectra_searched_ );
  carp(CARP_INFO, "%d spectrum-charge combinations loaded, %d spectrum-charge combinations searched. ", num_spectra_, num_spectra_searched_);
  

  // Delete temporary spectrumrecords file
 for (vector<TideSearchApplication::InputFile>::iterator original_file_name = inputFiles_.begin(); original_file_name != inputFiles_.end(); ++original_file_name) {
    carp(CARP_DEBUG, "Deleting %s", (*original_file_name).SpectrumRecords.c_str());
    remove((*original_file_name).SpectrumRecords.c_str());
  }

  // Delete stuffs
  if (out_tsv_target_ != NULL)
    delete out_tsv_target_;
  if (out_tsv_decoy_ != NULL)
    delete out_tsv_decoy_;
  if (out_mztab_target_ != NULL)
    delete out_mztab_target_;
  if (out_mztab_decoy_ != NULL)
    delete out_mztab_decoy_;
  if (out_pin_target_ != NULL)
    delete out_pin_target_;
  if (out_pin_decoy_ != NULL)
    delete out_pin_decoy_;

  return 0;
}

void TideSearchApplication::spectrum_search(void *threadarg) {  
  struct thread_data *my_data = (struct thread_data *) threadarg;

  ActivePeptideQueue* active_peptide_queue = my_data->active_peptide_queue_;
  int thread_id = my_data->thread_id_;

  int input_file_source;
  pb::Spectrum pb_spectrum;  
  while (true){

    // Get the next spectrum records with the smallest neutral mass from the heap and load the next spectrum records from the input files.
    locks_array_[LOCK_SPECTRUM_READING]->lock();
    if (spectrum_heap_.size() == 0) {
      locks_array_[LOCK_SPECTRUM_READING]->unlock();
      return;
    }
    // access the lightest spectra in the heap
    auto spectrum_pair = spectrum_heap_.front();   
    input_file_source = spectrum_pair.second;
    string spectrum_file_name = inputFiles_[input_file_source].OriginalName;

    // remove the lightest spectra from the heap.
    pop_heap(spectrum_heap_.begin(), spectrum_heap_.end(), compare_spectrum());
    spectrum_heap_.pop_back();   

    // read the next spectra from the input files and put it in the heap
    if (!spectrum_reader_[input_file_source]->Done()) {
      spectrum_reader_[input_file_source]->Read(&pb_spectrum);
      spectrum_heap_.push_back(make_pair(pb_spectrum, input_file_source));
      push_heap(spectrum_heap_.begin(), spectrum_heap_.end(), compare_spectrum());
    }
    if ( !spectrum_reader_[input_file_source]->OK() ){
      carp(CARP_FATAL, "Spectrum records file %s is corrupt.", spectrum_file_name.c_str());
    }
    ++num_spectra_;
    if (print_interval_ > 0 && num_spectra_ > 0 && num_spectra_ % print_interval_ == 0) {
      carp(CARP_INFO, "%d spectrum-charge combinations searched.", num_spectra_);
    }

    pb_spectrum = spectrum_pair.first;
    locks_array_[LOCK_SPECTRUM_READING]->unlock();
    
    Spectrum* spectrum = new Spectrum(pb_spectrum); 
    
    int charge = spectrum->ChargeState(0);
    double neutral_mass = pb_spectrum.neutral_mass();
    SpectrumCollection::SpecCharge* sc = new SpectrumCollection::SpecCharge(neutral_mass, charge, spectrum, 0);
   
    // Search one spectrum against its candidate peptides

    double precursor_mz = spectrum->PrecursorMZ();
    int scan_num = spectrum->SpectrumNumber();

    if (precursor_mz < spectrum_min_mz_|| 
        precursor_mz > spectrum_max_mz_ || 
        scan_num < min_scan_ || 
        scan_num > max_scan_ ||
        spectrum->Size() < min_peaks_  ||
        charge < min_precursor_charge_ || 
        charge >max_precursor_charge_ ) {
        delete spectrum;
        delete sc;
      
      continue; 
   }

    if (spectrum_flag_ != NULL) {  // TODO: Do something, possibly for cascade search
    }

    double min_range, max_range;
    vector<double>* min_mass = new vector<double>();
    vector<double>* max_mass = new vector<double>();
    
    computeWindow(*sc, min_mass, max_mass, &min_range, &max_range);
    active_peptide_queue->SetActiveRange(min_mass, max_mass, min_range, max_range);
    delete min_mass;
    delete max_mass;

    if (active_peptide_queue->nCandPeptides_ == 0) { // No peptides to score.
      delete spectrum;
      delete sc;  
      continue; 
    }
    long num_range_skipped = 0;
    long num_precursors_skipped = 0;
    long num_isotopes_skipped = 0;
    long num_retained = 0;

    ObservedPeakSet observed(use_neutral_loss_peaks_, use_flanking_peaks_);
    observed.PreprocessSpectrum(*(sc->spectrum), charge, &num_range_skipped,
      &num_precursors_skipped,
      &num_isotopes_skipped, &num_retained);

    locks_array_[LOCK_CANDIDATES]->lock();
    total_candidate_peptides_ += active_peptide_queue->nCandPeptides_;
    ++num_spectra_searched_;    
    num_range_skipped_ += num_range_skipped;
    num_precursors_skipped_ += num_precursors_skipped;
    num_isotopes_skipped_ += num_isotopes_skipped;
    num_retained_ += num_retained;
    locks_array_[LOCK_CANDIDATES]->unlock();  
 
    // allocate PSMscores for N scores
    TideMatchSet psm_scores(active_peptide_queue, &observed);  //nPeptides_ includes acitve and inacitve peptides

    // Calculate the scores needed
    switch (curScoreFunction_) {
      case PVALUES:
        PValueScoring(sc, active_peptide_queue, psm_scores);
        //break; // Run standard xcorr scoring in case of combined p-value calculations
      case XCORR_SCORE:
        // Spectrum preprocessing for xcorr scoring
        XCorrScoring(sc->charge, observed, active_peptide_queue, psm_scores);
        break;
      // case HYPERSCOR: TODO add new scoring functions here
    } 
    // Print the top-N results to the output files, 
    // The delta_cn, delta_lcn, repeat_ion_match, and tailor score calculation happens in PrintResults
    PrintResults(sc, spectrum_file_name, input_file_source, &psm_scores);

    delete spectrum;
    delete sc;
  }
}

void TideSearchApplication::XCorrScoring(int charge, ObservedPeakSet& observed, ActivePeptideQueue* active_peptide_queue, TideMatchSet& psm_scores){

  // Score the inactive peptides in the peptide queue if the number of nCadPeptides 
  // is less than the minimum. This is needed for Tailor scoring to get enough PSMS scores for statistics
  bool score_inactive_peptides = true;
  if (active_peptide_queue->min_candidates_ < active_peptide_queue->nCandPeptides_)
    score_inactive_peptides = false;
  
  //Actual Xcorr Scoring        
  int cnt = 0;
  for (deque<Peptide*>::const_iterator iter = active_peptide_queue->begin_; 
    iter != active_peptide_queue->end_;
    ++iter, ++cnt) {
    if ((*iter)->active_ == false && score_inactive_peptides == false) 
      continue;
    int xcorr = 0;
    int match_cnt = 0;
    int temp = 0;
    // int repeat_ion_match = 0;
  
    // Score with single charged b-y ion theoretical peaks
    xcorr += PeakMatching(observed, (*iter)->peaks_0, match_cnt, temp);

    if (charge > 2){
      // Score with double charged b-y ion theoretical peaks
      xcorr += PeakMatching(observed, (*iter)->peaks_1, match_cnt, temp);      
    }
    psm_scores.psm_scores_[cnt].peptide_itr_ = iter;
    psm_scores.psm_scores_[cnt].ordinal_ = cnt;    
    psm_scores.psm_scores_[cnt].xcorr_score_ = (double)xcorr/XCORR_SCALING;
    psm_scores.psm_scores_[cnt].by_ion_matched_ = match_cnt;
    psm_scores.psm_scores_[cnt].active_ = (*iter)->active_;
    psm_scores.psm_scores_[cnt].by_ion_total_ = (*iter)->peaks_0.size();
    if (charge > 2){
      psm_scores.psm_scores_[cnt].by_ion_total_ += (*iter)->peaks_1.size();
    }
  } 
}

int TideSearchApplication::PeakMatching(ObservedPeakSet& observed, vector<unsigned int>& peak_list, int& matching_peaks, int& repeat_matching_peaks) {
  bool previous_ion_matched = false;
  int score = 0;
  const int* cache = observed.GetCache(); 
  int end = observed.getCacheEnd();  

  for (vector<unsigned int>::const_iterator iter_uint = peak_list.begin(); iter_uint != peak_list.end(); iter_uint++) {
    if (*iter_uint >= end)
      continue;
    score += cache[*iter_uint];    // sum of the intensity of matching peaks, the xcorr score
    if (previous_ion_matched && cache[*iter_uint] > 0) {
      ++repeat_matching_peaks;
    }
    if (cache[*iter_uint] > 0) {
      previous_ion_matched = true;
      ++matching_peaks;
    } else {
      previous_ion_matched = false;
    }
  }
  return score;
}

void TideSearchApplication::PValueScoring(const SpectrumCollection::SpecCharge* sc, ActivePeptideQueue* active_peptide_queue, TideMatchSet& psm_scores){
  // 1. Calculate the REFACTORED XCORR SCORE and its EXACT P-VALUE

  // preprocess spectrum for refactored XCorr score calculation
  long num_range_skipped = 0;
  long num_precursors_skipped = 0;
  long num_isotopes_skipped = 0;
  long num_retained = 0;

  //For each candidate peptide, determine which discretized mass bin it is in
  //pepMassInt contains the corresponding mass bin for each candidate peptide
  //pepMassIntUnique contains the unique set of mass bins that candidate peptides fall in 
  vector<int> pepMassInt;
  pepMassInt.reserve(active_peptide_queue->nPeptides_);
  vector<int> pepMassIntUnique;
  pepMassIntUnique.reserve(active_peptide_queue->nPeptides_);
  getMassBin(pepMassInt, pepMassIntUnique, active_peptide_queue);
  int nPepMassIntUniq = (int)pepMassIntUnique.size();
  int maxPrecurMassBin = MassConstants::mass2bin(sc->neutral_mass + 250);

  vector< vector<int> > evidenceObs(nPepMassIntUniq, vector<int>(maxPrecurMassBin, 0));
  for (int pe = 0; pe < nPepMassIntUniq; pe++) {
    int pepMaInt = pepMassIntUnique[pe]; // TODO should be accessed with an iterator

    //preprocess to create one integerized evidence vector for each cluster of masses among selected peptides
    double pepMassMonoMean = (pepMaInt - 0.5 + bin_offset_) * bin_width_;
    evidenceObs[pe] = sc->spectrum->CreateEvidenceVectorDiscretized(
      bin_width_, bin_offset_, sc->charge, pepMassMonoMean, maxPrecurMassBin,
      &num_range_skipped, &num_precursors_skipped, &num_isotopes_skipped, &num_retained);
  }
  // The num_range_skipped etc. counts are counted only in the XCorr scoring. 

  // Calculate the null distribution OF PSM scores with dynamic programming method
  vector<double> nullDistribution;  // The score null distribution comes to this vector.
  //Score offset indicates the position in the vector corresponding to the score value 0.
  int score_offset = calcScoreCount(pepMassIntUnique, evidenceObs, nullDistribution);   

  // Calculate refactored XCorr scores
  // between a spectrum and all possible peptide candidates
  int scoreRefactInt;
  int cnt = 0;
  for (deque<Peptide*>::const_iterator iter = active_peptide_queue->begin_; 
      iter != active_peptide_queue->end_; 
      ++iter, ++cnt) {

      // if (!active_peptide_queue->candidatePeptideStatus_[cnt])
      if (! (*iter)->active_)
        continue;

    int pepMassIntIdx = 0;
    int curPepMassInt;
    for (int ma = 0; ma < nPepMassIntUniq; ++ma ) { //TODO should probably use iterator instead
      if (pepMassIntUnique[ma] == pepMassInt[cnt]) { //TODO pepMassIntUnique should be accessed with an interator
        pepMassIntIdx = ma;
        curPepMassInt = pepMassIntUnique[ma];
        break;
      }
    }
    // The actual scoring. Refactored XCorr Score calculation
    scoreRefactInt = 0;
    for (vector<unsigned int>::const_iterator iter_uint = (*iter)->peaks_1b.begin(); iter_uint != (*iter)->peaks_1b.end(); iter_uint++) {
      if (*iter_uint < maxPrecurMassBin)
        scoreRefactInt += evidenceObs[pepMassIntIdx][*iter_uint];
    }

    // Get the p-value of the refactored xcorr score 
    double pValue_xcorr = 1.0;
    if ((int)(scoreRefactInt) + score_offset >= 0) {
      pValue_xcorr = nullDistribution[(int)(scoreRefactInt) + score_offset];
    }

    // Store the scores. 
    psm_scores.psm_scores_[cnt].ordinal_ = cnt;
    psm_scores.psm_scores_[cnt].refactored_xcorr_ = scoreRefactInt / RESCALE_FACTOR;
    psm_scores.psm_scores_[cnt].exact_pval_ = pValue_xcorr;
    psm_scores.psm_scores_[cnt].active_ = (*iter)->active_;
  }

  // // 2. Calculate the RES-EV SCORE and its RES-EV P-VALUE, developed by Andy Lin
  //TODO so this includes ALL amino acids seen (including modified, NTerm mod, CTerm Mod)
  //as a result -- we will look for NTerm mod amino acids throughout spectrum instead of
  //just amino acids without NTerm mod

  //Creates a 3D vector representing 3D matrix
  //nPepMassIntUniq: number of mass bins candidate are in
  //nAARes: number of amino acids
  //maxPrecurMassBin: max number of mass bins
  int nAARes = iAAMass_.size();
  vector<vector<vector<double> > > residueEvidenceMatrix(nPepMassIntUniq,
      vector<vector<double> >(nAARes, vector<double>(maxPrecurMassBin, 0)));

  //Stores the score offset needed calculating res-ev p-values
  vector<int> scoreResidueOffsetObs(nPepMassIntUniq, -1);

  //For each mass bin, a vector hold the p-values for each corresponding res-ev score
  vector<vector<double> > pValuesResidueObs(maxPrecurMassBin);

//  map<int, bool> calcDPMatrix; //for each precursor mass bin, bool determines whether to calc DP matrix
  vector<bool> calcDPMatrix(nPepMassIntUniq, false); //for each precursor mass bin, bool determines whether to calc DP matrix

  //Create a residue evidence matrix and evidence vector
  //for each mass bin candidate peptides are in
  ObservedPeakSet observed(use_neutral_loss_peaks_, use_flanking_peaks_);

  for (int pe = 0; pe < nPepMassIntUniq; pe++) {
    // note: dAAMass_ contains amino acids masses in double form
    // precursorMass is the neutral mass
    observed.CreateResidueEvidenceMatrix(*(sc->spectrum), sc->charge, maxPrecurMassBin, sc->neutral_mass,
                                          nAARes, dAAMass_, fragTol_, granularityScale_,
                                          nTermMass_, cTermMass_, &num_range_skipped, 
                                          &num_precursors_skipped, &num_isotopes_skipped, &num_retained,
                                          residueEvidenceMatrix[pe]);
    vector<vector<double> > curResidueEvidenceMatrix = residueEvidenceMatrix[pe];

    //Get rid of values larger than curPepMassInt
    int curPepMassInt = pepMassIntUnique[pe];
    for (int i = 0; i < curResidueEvidenceMatrix.size(); i++) {
      curResidueEvidenceMatrix[i].resize(curPepMassInt);
    }
    residueEvidenceMatrix[pe] = curResidueEvidenceMatrix;
  }

  //Calculates a residue evidence score 
  //between a spectrum and all possible peptide candidates
  //based upon the residue evidence matrix and the theoretical spectrum
  int scoreResidueEvidence;
  vector<int> resEvScores;
  cnt = 0;
  for (deque<Peptide*>::const_iterator iter = active_peptide_queue->begin_; 
      iter != active_peptide_queue->end_; 
      ++iter, ++cnt) {

    if (!(*iter)->active_) {
      resEvScores.push_back(-1);
      continue;
    }

    int pepMassIntIdx = 0;
    int curPepMassInt = 0;
    for (int ma = 0; ma < nPepMassIntUniq; ma++ ) { //TODO should probably use iterator instead
      if (pepMassIntUnique[ma] == pepMassInt[cnt]) { //TODO pepMassIntUnique should be accessed with an interator
        pepMassIntIdx = ma;
        curPepMassInt = pepMassIntUnique[ma];
        break;
      }
    }     

    vector<vector<double> > curResidueEvidenceMatrix = residueEvidenceMatrix[pepMassIntIdx];
    scoreResidueEvidence = calcResEvScore(curResidueEvidenceMatrix, (*iter));
    resEvScores.push_back(scoreResidueEvidence);

    if (scoreResidueEvidence > 0) { // if > 0, set bool to true to create DP matrix
      calcDPMatrix[pepMassIntIdx] = true;
    }
  }
  //Create dyanamic programming matrix if there is a res-ev score greater than 0
  //and if user specified as a score function either 'residue-evidence matrix' or 'both'
  for (int pe=0 ; pe < nPepMassIntUniq; pe++) {
    int curPepMassInt = pepMassIntUnique[pe];
    if (calcDPMatrix[pe] == false) {
      continue;
    }

    vector<vector<double> > curResidueEvidenceMatrix = residueEvidenceMatrix[pe];
    vector<int> maxColEvidence(curPepMassInt, 0);

    //maxColEvidence is edited by reference
    int maxEvidence = getMaxColEvidence(curResidueEvidenceMatrix, maxColEvidence, curPepMassInt);
    int maxNResidue = floor((double)curPepMassInt / dAAMass_[0]); // dAAMass_[0] is the mass of the lightest amino acid Glutamine, (i.e. Glutamine 57)

    std::sort(maxColEvidence.begin(), maxColEvidence.end(), greater<int>());
    int maxScore = 0;
    for(int i = 0; i < maxNResidue; i++) { //maxColEvidence has been sorted
      maxScore += maxColEvidence[i];
    }

    int scoreOffset;
    vector<double> scoreResidueCount;

    calcResidueScoreCount(curPepMassInt, curResidueEvidenceMatrix, maxEvidence, maxScore,
                          scoreResidueCount, scoreOffset);
    scoreResidueOffsetObs[pe] = scoreOffset;

    double totalCount = 0;
    for (int i=scoreOffset ; i < scoreResidueCount.size(); i++) {
      totalCount += scoreResidueCount[i];
    }
    for (int i=scoreResidueCount.size()-2 ; i > -1; i--) {
      scoreResidueCount[i] = scoreResidueCount[i] + scoreResidueCount[i+1];
    }
    for (int i = 0; i < scoreResidueCount.size(); i++) {
      //Avoid potential underflow
      scoreResidueCount[i] = exp(log(scoreResidueCount[i]) - log(totalCount));
    }
    pValuesResidueObs[curPepMassInt] = scoreResidueCount;
  }

  /************ calculate p-values for PSMs using residue evidence matrix ****************/
  int curPepMassInt;
  double pValue_xcorr;
  double pValue_resEv;
  double pValue_combined = 0.3;

  cnt = 0;
  for (deque<Peptide*>::const_iterator iter = active_peptide_queue->begin_; 
      iter != active_peptide_queue->end_; 
      ++iter, ++cnt) {

    if (! (*iter)->active_ )
      continue;

    int pepMassIntIdx = 0;

    for (int ma = 0; ma < nPepMassIntUniq; ma++ ) { //TODO should probably use iterator instead
      if (pepMassIntUnique[ma] == pepMassInt[cnt]) { //TODO pepMassIntUnique should be accessed with an interator
        pepMassIntIdx = ma;
        curPepMassInt = pepMassIntUnique[ma];
        break;
      }
    }

    int scoreCountIdx;

    scoreResidueEvidence = resEvScores[cnt];
    if (calcDPMatrix[pepMassIntIdx]) {
      scoreCountIdx = scoreResidueEvidence + scoreResidueOffsetObs[pepMassIntIdx];
      pValue_resEv = pValuesResidueObs[curPepMassInt][scoreCountIdx];
    } else {
      pValue_resEv = 1.0;
    }
    // if (pValue_resEv == 0.0)
    //   carp(CARP_FATAL, "PSM Res-EV p-value should not be equal to 0.0");

    /************ calculate comnbined p-values  ****************/    
    pValue_xcorr = psm_scores.psm_scores_[cnt].exact_pval_;  // get the p-value od the XCorr from previous calculations
    if (pValue_resEv > 0.0) {
      double cPval = pValue_xcorr * pValue_resEv;    
      pValue_combined = calcCombinedPval(1.2, cPval, 2); //2 is the # of p-values that are combined // 1.2 value has been empircally determined for combining two independent p-values
    } else {
      pValue_combined = pValue_xcorr; //in case resEV pval == 0
    }

    // Store the scores. 
    psm_scores.psm_scores_[cnt].resEv_score_   = scoreResidueEvidence;
    psm_scores.psm_scores_[cnt].resEv_pval_    = pValue_resEv;
    psm_scores.psm_scores_[cnt].combined_pval_ = pValue_combined;
    psm_scores.psm_scores_[cnt].ordinal_       = cnt;
    psm_scores.psm_scores_[cnt].active_        = (*iter)->active_;  
  }
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
  Peptide* curPeptide
) {
  //Make sure the number of theoretical peaks match pepLen
  int pepLen = curPeptide->Len();

  int scoreResidueEvidence = 0;
  vector<double> residueMasses = curPeptide->getAAMasses(); //retrieves the amino acid masses, modifications included
  for (int res = 0; res < pepLen - 1; res++) {
    // Perform some rounding, because the AA masses are rounded in dAAMass
    double tmpAAMass = MassConstants::ToDouble(MassConstants::ToFixPt(residueMasses[res]));  

    vector<double>::const_iterator mass_itr = find(dAAMass_.begin(), dAAMass_.end(), tmpAAMass);

    if (mass_itr == dAAMass_.end()){
      // for (vector<double>::const_iterator mass_itr = dAAMass_.begin(); mass_itr != dAAMass_.end(); ++mass_itr) {
      //   carp(CARP_INFO, "AA mass: %lf", *mass_itr);
      // }
      carp(CARP_FATAL, "'%lf' does not exist. residue mass: %lf", tmpAAMass, residueMasses[res]);
    }
    
    int tmpAA = mass_itr - dAAMass_.begin();
    if ( curPeptide->peaks_1b.size() > res && curResidueEvidenceMatrix[tmpAA].size() >  curPeptide->peaks_1b[res]-1) {
      scoreResidueEvidence += curResidueEvidenceMatrix[tmpAA][curPeptide->peaks_1b[res]-1];
    }
  }
  return scoreResidueEvidence;
}



// Calculate the null distribution with dynamic programming
int TideSearchApplication::calcScoreCount(vector<int>& pepMassIntUnique, vector<vector<int>>& evidenceObs, vector<double>& nullDistribution) {
  const int nDeltaMass = iAAMass_.size();
  int minDeltaMass = iAAMass_[0];
  int maxDeltaMass = iAAMass_[nDeltaMass - 1];

  int nPepMassIntUniq = evidenceObs.size();
  int maxPrecurMassBin = evidenceObs[0].size();
  // local variables
  int row;
  int col;
  int ma;
  int evidence;
  int de;
  int evidenceRow;
  double sumScore;

  // For each Unique Mass Int we calculate different 
  int* nRows = new int[nPepMassIntUniq];   // Stores the length (rows) of the dynamic programming table
  int* scoreOffsets = new int[nPepMassIntUniq];
  double** pValueScoreObs = new double*[nPepMassIntUniq];

  for (int pe = 0; pe < nPepMassIntUniq; ++pe) { 
    int pepMassInt = pepMassIntUnique[pe]; // TODO should be accessed with an iterator

    // NOTE: will have to go back to separate dynamic programming for
    //       target and decoy if they have different probNI and probC
    int maxEvidence = *std::max_element(evidenceObs[pe].begin(), evidenceObs[pe].end());
    int minEvidence = *std::min_element(evidenceObs[pe].begin(), evidenceObs[pe].end());

    // estimate maxScore and minScore
    int maxNResidue = (int)floor((double)pepMassInt / (double)minDeltaMass);
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

    // Initialize variables for the dynamic programming table
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

    // Allocate memory for the DP table
    double** dynProgArray = new double*[nCol];
    for (col = 0; col < colFirst+maxDeltaMass; col++) {
      dynProgArray[col] = new double[nRow];
      memset(dynProgArray[col], 0.0, nRow*sizeof(double));
    }
    dynProgArray[colStart][initCountRow] = 1.0; // initial count of peptides with mass = 1 
    vector<int> deltaMassCol(nDeltaMass);
    // populate matrix with scores for first (i.e. N-terminal) amino acid in sequence
    for (de = 0; de < nDeltaMass; de++) {
      ma = iAAMass_[de];
      col = colStart + ma;
      row = initCountRow + evidenceObs[pe][col];
      if (col <= maxDeltaMass + colLast) {
        dynProgArray[col][row] += dynProgArray[colStart][initCountRow] * dAAFreqN_[de];
      }
    }
    for (ma = 0; ma < colFirst; ++ma )
        delete dynProgArray[ma];
    for (ma = colFirst; ma < colLast; ma++) {
      dynProgArray[ma+maxDeltaMass ] = new double[nRow];
      memset(dynProgArray[ma+maxDeltaMass], 0.0, nRow*sizeof(double));
      
      for (de = 0; de < nDeltaMass; de++) {
        deltaMassCol[de] = ma + iAAMass_[de];
      }
      for (row = rowFirst; row <= rowLast; row++) {
        if (dynProgArray[ma][row] == 0.0) {
          continue;
        }
        for (de = 0; de < nDeltaMass; de++) {
          col = deltaMassCol[de];
          if (col < colLast) {
            evidenceRow = row + evidenceObs[pe][col];
            dynProgArray[col][evidenceRow] += dynProgArray[ma][row]*dAAFreqI_[de]; 
          } else if (col == colLast) { 
            evidenceRow = row;
            dynProgArray[col][evidenceRow] += dynProgArray[ma][row]*dAAFreqC_[de]; 
          }
        }
      }
      delete dynProgArray[ma];
    }
    // The final null distribution is stored in pValueScoreObs[pe]
    nRows[pe] = nRow;
    pValueScoreObs[pe] = new double[nRow];
    scoreOffsets[pe] = scoreOffsetObs;
    for (row = 0; row < nRow; ++row) {
      pValueScoreObs[pe][row] = dynProgArray[colLast][row];
    }
    // clean up
    for (col = colLast; col < colLast+maxDeltaMass; ++col) {
      delete [] dynProgArray[col];    
    }
    delete [] dynProgArray;
  }

  // Merge separate score distirbutions.
  int max_offset = 0;
  int max_row = 0;
  for (int pe = 0; pe < nPepMassIntUniq; ++pe) {
    if (max_offset < scoreOffsets[pe] )
      max_offset = scoreOffsets[pe];
    if (max_row < nRows[pe])
      max_row = nRows[pe];
  }
  
  int score_idx;        
  max_row += 1;
  nullDistribution.reserve(max_row*2);
  for (int i = 0; i < max_row*2; ++i)
    nullDistribution[i] = 0.0;
  std::fill(nullDistribution.begin(), nullDistribution.begin()+max_row*2, 0);
  double *scoreCountBinAdjust = new double[max_row*2];
  memset(scoreCountBinAdjust, 0.0, sizeof(double)*max_row*2);
  
  // Merges the separated partial score histograms.
  double totalCount = 0.0;
  for (int pe = 0 ; pe < nPepMassIntUniq; ++pe) {
    int offset_diff = max_offset - scoreOffsets[pe];
    for (score_idx = 0; score_idx < nRows[pe]; ++score_idx) {
      nullDistribution[score_idx + offset_diff] += pValueScoreObs[pe][score_idx];
      totalCount += pValueScoreObs[pe][score_idx];
    }
  }

  if (totalCount == 0.0) {
    for (score_idx = 0; score_idx < max_row*2; score_idx++) {
      nullDistribution[score_idx] = 1.0;
    }
  } else {      
    // Normalize the raw score histogram to get a probabilisitc distribution
    double logTotalCount = log(totalCount);
    for (score_idx = max_row*2-2; score_idx >= 0; --score_idx) {
      scoreCountBinAdjust[score_idx] = nullDistribution[score_idx]/2.0;
      nullDistribution[score_idx] += nullDistribution[score_idx + 1];
    }
    for (score_idx = max_row*2-2; score_idx >= 0; --score_idx) {
      nullDistribution[score_idx] -= scoreCountBinAdjust[score_idx];
    }
    // normalize distribution; use exp( log ) to avoid potential underflow
    for (score_idx = max_row*2-2; score_idx >= 0; --score_idx) {
      if (nullDistribution[score_idx] > 0.0) {
        nullDistribution[score_idx] = exp(log(nullDistribution[score_idx]) - logTotalCount);
      }
    }
  }
  delete [] scoreCountBinAdjust;  
  delete [] nRows;
  delete [] scoreOffsets;
  for (int pe = 0; pe < nPepMassIntUniq; ++pe)
    delete [] pValueScoreObs[pe];
  delete [] pValueScoreObs;
  return max_offset;
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
  int pepMassInt,
  vector<vector<double> >& residueEvidenceMatrix,
  int maxEvidence,
  int maxScore,
  vector<double>& scoreCount, //this is returned for later use
  int& scoreOffset //this is returned for later use
) {
  const int nAa = iAAMass_.size();
  int minAaMass = iAAMass_[0];
  int maxAaMass = iAAMass_[nAa - 1];

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
  int colStart = nTermMassBin_;
  int nRow = bottomRowBuffer - minScore + 1 + maxScore + topRowBuffer;
  int nCol = colBuffer + pepMassInt + maxAaMass;
  int rowFirst = bottomRowBuffer + 1;
  int rowLast = rowFirst - minScore + maxScore;
  int colFirst = colStart + 1;
  int colLast = pepMassInt - cTermMassBin_;
  int initCountRow = bottomRowBuffer - minScore + 1;
  int initCountCol = colStart;

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
    ma = iAAMass_[de];

    //&& -1 is to account for zero-based indexing in evidence vector
    //row = initCountRow + residueEvidueMatrix[ de ][ ma + nTermMass - 1 ]; //original
    row = initCountRow + residueEvidenceMatrix[de][ma + nTermMassBin_ - 1]; //+nTermMassBin for N-Term mod and -1 for 0 indexing

    //TODO need to change this to based off bool
    // if (nTermMassBin_ == 1) { //N-Term not modified
    //   col = initCountCol + ma;
    // } else { //N-Term is modified
    //   col = initCountCol + ma - nTermMassBin_ + 1;
    // }
    col = initCountCol + ma - nTermMassBin_ + 1;

    if (col <= maxAaMass + colLast && col >= initCountCol) { //TODO not sure if below or above is correct
      dynProgArray[row][col] += dynProgArray[initCountRow][initCountCol] * dAAFreqN_[de];
    }
  }

  // Set to zero now that score counts for first amino acid are in matrix
  dynProgArray[initCountRow][initCountCol] = 0.0;

  //  The following code was reorganized by AKF to make the DP calculation ~3 times faster
  int newCol;
  for (ma = initCountCol+1; ma < colLast; ++ma) {
    col = ma;  
    for (de = 0; de < nAa; de++) {
      aaMassCol[de] = col + iAAMass_[de];
    }
    for (row = rowFirst; row <= rowLast; ++row) {
      if (dynProgArray[row][col] == 0.0) {
        continue;
      }
      for (de = 0; de < nAa; de++) {
        newCol = aaMassCol[de];
        if (newCol < colLast) {
          evidRow = row + residueEvidenceMatrix[de][newCol];
          dynProgArray[evidRow][newCol] += dynProgArray[row][col] * dAAFreqI_[de];
        } else if (newCol == colLast) {
          evidRow = row;
          dynProgArray[evidRow][newCol] += dynProgArray[row][col] * dAAFreqC_[de];            
        } 
      }        
    }      
  } 
  int colScoreCount = colLast;

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


void TideSearchApplication::computeWindow(
  const SpectrumCollection::SpecCharge& sc,
  vector<double>* out_min,
  vector<double>* out_max,
  double* min_range,
  double* max_range
) {
  double unit_dalton = BIN_WIDTH;
  switch (window_type_) {
  case WINDOW_MASS:
    for (vector<int>::const_iterator ie = negative_isotope_errors_.begin(); ie != negative_isotope_errors_.end(); ++ie) {
      out_min->push_back(sc.neutral_mass + (*ie * unit_dalton) - precursor_window_);
      out_max->push_back(sc.neutral_mass + (*ie * unit_dalton) + precursor_window_);
    }
    *min_range = (sc.neutral_mass + (negative_isotope_errors_.front() * unit_dalton)) - precursor_window_;
    *max_range = (sc.neutral_mass + (negative_isotope_errors_.back() * unit_dalton)) + precursor_window_;
    break;
  case WINDOW_MZ: {
    double mz_minus_proton = sc.spectrum->PrecursorMZ() - MASS_PROTON;
    for (vector<int>::const_iterator ie = negative_isotope_errors_.begin(); ie != negative_isotope_errors_.end(); ++ie) {
      out_min->push_back((mz_minus_proton - precursor_window_) * sc.charge + (*ie * unit_dalton));
      out_max->push_back((mz_minus_proton + precursor_window_) * sc.charge + (*ie * unit_dalton));
    }
    // N.B. Next two lines used to have max_charge instead of sc.charge. Not sure why. --WSN 9 Jun 2021
    *min_range = (mz_minus_proton*sc.charge + (negative_isotope_errors_.front() * unit_dalton)) - precursor_window_*sc.charge;
    *max_range = (mz_minus_proton*sc.charge + (negative_isotope_errors_.back() * unit_dalton)) + precursor_window_*sc.charge;
    break;
  }
  case WINDOW_PPM: {
    double tiny_precursor = precursor_window_ * 1e-6;
    for (vector<int>::const_iterator ie = negative_isotope_errors_.begin(); ie != negative_isotope_errors_.end(); ++ie) {
      out_min->push_back((sc.neutral_mass + (*ie * unit_dalton)) * (1.0 - tiny_precursor));
      out_max->push_back((sc.neutral_mass + (*ie * unit_dalton)) * (1.0 + tiny_precursor));
    }
    *min_range = (sc.neutral_mass + (negative_isotope_errors_.front() * unit_dalton)) * (1.0 - tiny_precursor);
    *max_range = (sc.neutral_mass + (negative_isotope_errors_.back() * unit_dalton)) * (1.0 + tiny_precursor);
    break;
  }
  default:
    carp(CARP_FATAL, "Invalid window type");
  }
  carp(CARP_DETAILED_DEBUG, "Scan=%d Charge=%d Mass window=[%f, %f]",
       sc.spectrum->SpectrumNumber(), sc.charge, (*out_min)[0], (*out_max)[0]);
}

vector<int> TideSearchApplication::getNegativeIsotopeErrors() {
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


void TideSearchApplication::getPeptideIndexData(const string input_index, ProteinVec& proteins, vector<const pb::AuxLocation*>& locations, pb::Header& peptides_header){

  string peptides_file = FileUtils::Join(input_index, "pepix");  
  string proteins_file = FileUtils::Join(input_index, "protix");
  string auxlocs_file = FileUtils::Join(input_index, "auxlocs");  
  string residue_stats_file = FileUtils::Join(input_index, "residue_stat");  

  // Read protein index file
  carp(CARP_INFO, "Reading index %s", input_index.c_str());

  pb::Header protein_header;

  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins, proteins_file, &protein_header)) {
    carp(CARP_FATAL, "Error reading index (%s)", proteins_file.c_str());
  }
  // There shouldn't be more than one header in the protein pb.
  pb::Header_Source headerSource = protein_header.source(0);  
  string decoy_prefix = "";
  if (headerSource.has_decoy_prefix()){
    decoy_prefix = headerSource.decoy_prefix();
  } else {
    carp(CARP_WARNING, "You are using an index generated by an old version of tide-index."
                       "This will not affect your results, but this index may need to be "
                       "re-created to work with future versions of tide-index. ");
  }
  TideMatchSet::decoy_prefix_ = decoy_prefix;  
  if (headerSource.has_filename()){
    TideMatchSet::fasta_file_name_ = headerSource.filename();
  }
  
  // Read auxlocs index file
  ReadRecordsToVector<pb::AuxLocation>(&locations, auxlocs_file);

  // Read peptides index file

  if ((peptides_header.file_type() != pb::Header::PEPTIDES) ||
      !peptides_header.has_peptides_header()) {
    carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
  }

  const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();
  decoy_num_ = pepHeader.has_decoys_per_target() ? pepHeader.decoys_per_target() : 0;

  // Initizalize the Mass Constants class
  MassConstants::Init(&pepHeader.mods(), &pepHeader.nterm_mods(), &pepHeader.cterm_mods(),
      &pepHeader.nprotterm_mods(), &pepHeader.cprotterm_mods(), bin_width_, bin_offset_);

  if (curScoreFunction_ == PVALUES) {

    // Get the terminal mass bins for RES-EV
    //TODO assumption is that there is one nterm static mod per peptide
    if (pepHeader.nterm_mods().static_mod_size() > 0) {
      nTermMassBin_ = MassConstants::mass2bin(
                        MassConstants::mono_h + pepHeader.nterm_mods().static_mod(0).delta());
      nTermMass_ = MassConstants::mono_h + pepHeader.nterm_mods().static_mod(0).delta();
    } else {
      nTermMassBin_ = MassConstants::mass2bin(MassConstants::mono_h);
      nTermMass_ = MassConstants::mono_h;
    }

    //TODO assumption is that there is one cterm mod per peptide
    if (pepHeader.cterm_mods().static_mod_size() > 0) {
      cTermMassBin_ = MassConstants::mass2bin(MassConstants::mono_oh + pepHeader.cterm_mods().static_mod(0).delta());
      cTermMass_ = MassConstants::mono_oh + pepHeader.cterm_mods().static_mod(0).delta();
    } else {
      cTermMassBin_ = MassConstants::mass2bin(MassConstants::mono_oh);
      cTermMass_ = MassConstants::mono_oh;
    }

    if (FileUtils::Exists(residue_stats_file)) {  // The amino acid frequency calculation is done in tide-index.
      RecordReader residue_stat_reader = RecordReader(residue_stats_file);
      bool done;
      while (!(done = residue_stat_reader.Done())) {
        pb::ResidueStats last_residue_stat;
        residue_stat_reader.Read(&last_residue_stat);    
        double aa_mass = last_residue_stat.aamass();
        dAAMass_.push_back(aa_mass);
        dAAFreqN_.push_back(last_residue_stat.aafreqn());
        dAAFreqI_.push_back(last_residue_stat.aafreqi());
        dAAFreqC_.push_back(last_residue_stat.aafreqc());
        iAAMass_.push_back(MassConstants::mass2bin(aa_mass));
        string aa_str = last_residue_stat.aa_str();
        mMass2AA_[aa_mass] = aa_str;
      }
    } else {
      carp(CARP_INFO, "You are using an old format of the peptide index data. Please recreate your peptide index data");
      // The following (this whole else branch) part should be removed later.
      // Calculate the Amino Acid Frequencies for the P-value calculation if it wasn't done with tide-index.
      unsigned int len;
      unsigned int i;
      unsigned int cntTerm = 0;
      unsigned int cntInside = 0;
      string peptide_seq;
      unsigned int residue_bin;
      string tempAA;
      int mod_precision  = Params::GetInt("mod-precision");
      const unsigned int MaxModifiedAAMassBin = MassConstants::ToFixPt(2000.0);;   //2000 is the maximum mass of a modified amino acid
      unsigned int* nvAAMassCounterN = new unsigned int[MaxModifiedAAMassBin];   //N-terminal amino acids
      unsigned int* nvAAMassCounterC = new unsigned int[MaxModifiedAAMassBin];   //C-terminal amino acids
      unsigned int* nvAAMassCounterI = new unsigned int[MaxModifiedAAMassBin];   //inner amino acids in the peptides
      memset(nvAAMassCounterN, 0, MaxModifiedAAMassBin * sizeof(unsigned int));
      memset(nvAAMassCounterC, 0, MaxModifiedAAMassBin * sizeof(unsigned int));
      memset(nvAAMassCounterI, 0, MaxModifiedAAMassBin * sizeof(unsigned int));

      pb::Peptide current_pb_peptide_;
      HeadedRecordReader peptide_reader = HeadedRecordReader(peptides_file, &peptides_header);
      while (!(peptide_reader.Done())) { //read all peptides form index
        peptide_reader.Read(&current_pb_peptide_);
        Peptide peptide(current_pb_peptide_, proteins);
        vector<double> residue_masses = peptide.getAAMasses(); //retrieves the amino acid masses, modifications included
        peptide_seq = peptide.Seq();
        len = current_pb_peptide_.length();

        vector<double> residue_mods(len, 0);  // Initialize a vecotr of peptide length  with zeros. 
        
        // Handle variable modifications
        if (current_pb_peptide_.has_nterm_mod()){ // Handle N-terminal modifications
          int index;
          double delta;
          MassConstants::DecodeMod(ModCoder::Mod(current_pb_peptide_.nterm_mod()), &index, &delta);
          residue_mods[index] = delta;
        }

        for (i = 0; i < current_pb_peptide_.modifications_size(); ++i) {
          int index;
          double delta;
          MassConstants::DecodeMod(current_pb_peptide_.modifications(i), &index, &delta);
          residue_mods[index] = delta;
        }
        if (current_pb_peptide_.has_cterm_mod()){  // Handle C-terminal modifications
          int index;
          double delta;
          MassConstants::DecodeMod(ModCoder::Mod(current_pb_peptide_.cterm_mod()), &index, &delta);
          residue_mods[index] = delta;
        }
        // count AA masses
        residue_bin = MassConstants::ToFixPt(residue_masses[0]);
        ++nvAAMassCounterN[residue_bin];  // N-temrianl
        if (nvAAMassCounterN[residue_bin] == 1){
          tempAA = peptide_seq[0];
          if (residue_mods[0] != 0) {
            tempAA += "[" + StringUtils::ToString(residue_mods[0], mod_precision) + ']';
          }
          mMass2AA_[MassConstants::ToDouble(residue_bin)] = tempAA;
        }
        for (i = 1; i < len-1; ++i) {
          residue_bin = MassConstants::ToFixPt(residue_masses[i]);
          ++nvAAMassCounterI[residue_bin];  // non-terminal
          if (nvAAMassCounterI[residue_bin] == 1){
            tempAA = peptide_seq[i];
            if (residue_mods[i] != 0) {
              tempAA += "[" + StringUtils::ToString(residue_mods[i], mod_precision) + ']';
            }
            mMass2AA_[MassConstants::ToDouble(residue_bin)] = tempAA;
          }		
          ++cntInside;
        }
        residue_bin = MassConstants::ToFixPt(residue_masses[len - 1]);
        ++nvAAMassCounterC[residue_bin];  // C-temrinal
        if (nvAAMassCounterC[residue_bin] == 1){
          tempAA = peptide_seq[len - 1];
          if (residue_mods[len - 1] != 0) {
            tempAA += "[" + StringUtils::ToString(residue_mods[len - 1], mod_precision) + ']';
          }
          mMass2AA_[MassConstants::ToDouble(residue_bin)] = tempAA;
        }
        ++cntTerm;
      }
      unsigned int uiUniqueMasses = 0;
      double aa_mass;
      for (i = 0; i < MaxModifiedAAMassBin; ++i) {
        if (nvAAMassCounterN[i] || nvAAMassCounterI[i] || nvAAMassCounterC[i]) {
          ++uiUniqueMasses;
          aa_mass = MassConstants::ToDouble(i);
          dAAMass_.push_back(aa_mass);
          dAAFreqN_.push_back((double)nvAAMassCounterN[i] / cntTerm);
          dAAFreqI_.push_back((double)nvAAMassCounterI[i] / cntInside);
          dAAFreqC_.push_back((double)nvAAMassCounterC[i] / cntTerm);
          iAAMass_.push_back(MassConstants::mass2bin(aa_mass));          
        }
      }
      // Write the statistics to a Protocol Buffer, so that next time it is at hand
      RecordWriter residue_stat_writer = RecordWriter(residue_stats_file);
      CHECK(residue_stat_writer.OK());  
      for (int i = 0; i < dAAMass_.size(); ++i){
        pb::ResidueStats last_residue_stat;
        last_residue_stat.set_aamass(dAAMass_[i]);
        last_residue_stat.set_aafreqn(dAAFreqN_[i]);
        last_residue_stat.set_aafreqi(dAAFreqI_[i]);
        last_residue_stat.set_aafreqc(dAAFreqC_[i]);
        string aa_str = mMass2AA_[dAAMass_[i]];
        last_residue_stat.set_aa_str(aa_str);
        CHECK(residue_stat_writer.Write(&last_residue_stat));
      }       
      delete[] nvAAMassCounterN;
      delete[] nvAAMassCounterI;
      delete[] nvAAMassCounterC;
    } // Finished calculating AA frequencies
  }  
}

string TideSearchApplication::getOutputFileName() {
  return output_file_name_;
}

// In order to add more options, you need to add them to ./src/util/Params.cpp
vector<string> TideSearchApplication::getOptions() const {
  string arr[] = {
    "auto-mz-bin-width",
    "auto-precursor-window",
    "concat",
    "deisotope",
    "elution-window-size",
    "fileroot",
    "fragment-tolerance",
    "isotope-error",
    "mass-precision",
    "max-precursor-charge",
    "min-precursor-charge",
    "min-peaks",
    "mod-precision",
    "mz-bin-offset",
    "mz-bin-width",
    "mzid-output",
    "mztab-output",
    "num-threads",
    "output-dir",
    "override-charges",
    "overwrite",
    "parameter-file",
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
    "score-function",
    "skip-preprocessing",
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
/*
Parameters to be removed from ./Params.cpp"
brief-output
peptide-centric-search
use-tailor-calibration
exact-p-value
charge-state
evidence-granularity
*/

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
    ParamMedic::RunAttributeResult errorCalcResult;
    ParamMedicApplication::processFiles(Params::GetStrings("tide spectra file"),
      true, false, &errorCalcResult, NULL);

    if (autoPrecursor != "false") {
      string fail = errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_FAILURE);
      if (fail.empty()) {
        double sigma = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_SIGMA));
        double prediction = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_PREDICTION));
        carp(CARP_INFO, "precursor ppm standard deviation: %f", sigma);
        carp(CARP_INFO, "precursor error estimate (ppm): %.2f", prediction);
        Params::Set("precursor-window", prediction);
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
        Params::Set("mz-bin-width", prediction);
      } else {
        carp(autoFragment == "fail" ? CARP_FATAL : CARP_ERROR,
             "failed to calculate fragment error: %s", fail.c_str());
      }
    }
  }  
  
}
void TideSearchApplication::getInputFiles(int thread_id) {
  // Try to read all spectrum files as spectrumrecords, convert those that fail
  if (thread_id > inputFiles_.size())
    return;
  for (vector<TideSearchApplication::InputFile>::iterator original_file_name = inputFiles_.begin()+thread_id; 
       original_file_name < inputFiles_.begin() + (inputFiles_.size()); 
       original_file_name = original_file_name + num_threads_) 
    {
    carp(CARP_DEBUG, "Start processing input files");
    bool keepSpectrumrecords = true;
    string original_name = (*original_file_name).OriginalName;
    string spectrumrecords = original_name;
    // Check if the input file is spectrum records of google protocol buffer
    pb::Header header;
    HeadedRecordReader reader(original_name, &header);
    if (header.file_type() != pb::Header::SPECTRA) {
      // converting to spectrumrecords file 

      carp(CARP_INFO, "Converting %s to spectrumrecords format", original_name.c_str());
      carp(CARP_DEBUG, "Elapsed time starting conversion: %.3g s", wall_clock() / 1e6);
      
      spectrumrecords = Params::GetString("store-spectra");
      keepSpectrumrecords = !spectrumrecords.empty();
      if (!keepSpectrumrecords) {
        spectrumrecords = make_file_path(FileUtils::BaseName( original_name) + ".spectrumrecords.tmp");
      } else if (inputFiles_.size() > 1) {
        carp(CARP_FATAL, "Cannot use store-spectra option with multiple input "
                         "spectrum files");
      }
      carp(CARP_DEBUG, "New spectrumrecords filename: %s", spectrumrecords.c_str());
      int spectra_num = 0;
      if (!SpectrumRecordWriter::convert(original_name, spectrumrecords, spectra_num)) {
        carp(CARP_FATAL, "Error converting %s to spectrumrecords format", original_name.c_str());
      }
      locks_array_[LOCK_SPECTRUM_READING]->lock();
      total_spectra_num_ += spectra_num;
      locks_array_[LOCK_SPECTRUM_READING]->unlock();

    }
    (*original_file_name).SpectrumRecords  = spectrumrecords;
    (*original_file_name).Keep = keepSpectrumrecords;
    carp(CARP_DEBUG, "Finish converting");
  }
}

void TideSearchApplication::createOutputFiles() {
  
  // Create output files for the search results
   
  bool overwrite = Params::GetBool("overwrite");  
  bool concat = Params::GetBool("concat");

  string concat_file_name;
  string target_file_name;
  string decoy_file_name;

  // Get output files for tsv format
  concat_file_name = make_file_path("tide-search.txt");
  target_file_name = make_file_path("tide-search.target.txt");
  decoy_file_name  = make_file_path("tide-search.decoy.txt");

  if (overwrite) {
    remove(concat_file_name.c_str());  
    remove(target_file_name.c_str());  
    remove(decoy_file_name.c_str());  
  }

  if (Params::GetBool("txt-output") == true) {  // original tide-search output format in tab-delimited text files (txt)
    string header = TideMatchSet::getHeader(TIDE_SEARCH_TSV, tide_index_mzTab_file_path_);

    if (concat) {

      out_tsv_target_ = create_stream_in_path(concat_file_name.c_str(), NULL, overwrite);
      output_file_name_ = concat_file_name;
      *out_tsv_target_ << header; 

    } else {

      out_tsv_target_ = create_stream_in_path(target_file_name.c_str(), NULL, overwrite);
      output_file_name_ = target_file_name;
      *out_tsv_target_ << header; 
      if (decoy_num_ > 0) {
        out_tsv_decoy_ = create_stream_in_path(decoy_file_name.c_str(), NULL, overwrite);
        *out_tsv_decoy_ << header;
      }
    }  
  }

  // Get output files for mzTAB format
  concat_file_name = make_file_path("tide-search.mzTab");
  target_file_name = make_file_path("tide-search.target.mzTab");
  decoy_file_name  = make_file_path("tide-search.decoy.mzTab");

  if (overwrite) {
    remove(concat_file_name.c_str());  
    remove(target_file_name.c_str());  
    remove(decoy_file_name.c_str());  
  }

  if ( Params::GetBool("mztab-output") == true) {   // mzTAB output format

    // Get column headers
    string header = TideMatchSet::getHeader(TIDE_SEARCH_MZTAB_TSV, tide_index_mzTab_file_path_);  // Gets the column headers
    if (concat) {

      out_mztab_target_ = create_stream_in_path(concat_file_name.c_str(), NULL, overwrite);
      output_file_name_ = concat_file_name;
      *out_mztab_target_ << header; 

    } else {

      out_mztab_target_ = create_stream_in_path(target_file_name.c_str(), NULL, overwrite);
      output_file_name_ = target_file_name;
      *out_mztab_target_ << header; 

      if (decoy_num_ > 0) {
        out_mztab_decoy_ = create_stream_in_path(decoy_file_name.c_str(), NULL, overwrite);
        *out_mztab_decoy_ << header;
      }
    }  
  }

}

void TideSearchApplication::PrintResults(const SpectrumCollection::SpecCharge* sc, string spectrum_file_name, int spectrum_file_cnt, TideMatchSet* psm_scores) {
  string concat_or_target_report;
  string decoy_report;

  if (out_mztab_target_ != NULL) {
    psm_scores->getReport(TIDE_SEARCH_MZTAB_TSV, spectrum_file_name, sc, spectrum_file_cnt, concat_or_target_report, decoy_report); 
    locks_array_[LOCK_RESULTS]->lock();
    *out_mztab_target_ << concat_or_target_report;
    locks_array_[LOCK_RESULTS]->unlock();

    if (out_mztab_decoy_ != NULL) {
      locks_array_[LOCK_RESULTS]->lock();
      *out_mztab_decoy_ << decoy_report;
      locks_array_[LOCK_RESULTS]->unlock();
    }
  }

  if ( out_tsv_target_ != NULL) {
    psm_scores->getReport(TIDE_SEARCH_TSV, spectrum_file_name, sc, spectrum_file_cnt, concat_or_target_report, decoy_report); 
    locks_array_[LOCK_RESULTS]->lock();
    *out_tsv_target_ << concat_or_target_report;
    locks_array_[LOCK_RESULTS]->unlock();
    
    if (out_tsv_decoy_ != NULL) {
      locks_array_[LOCK_RESULTS]->lock();
      *out_tsv_decoy_ << decoy_report;
      locks_array_[LOCK_RESULTS]->unlock();
    }
  }
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
  ActivePeptideQueue* active_peptide_queue
) {
  int pe = 0;
  for (deque<Peptide*>::const_iterator iter = active_peptide_queue->begin_;
      iter != active_peptide_queue->end_; 
      ++iter) {
    double pepMass = (*iter)->Mass();
    int pepMaInt = MassConstants::mass2bin(pepMass);
    pepMassInt[pe] = pepMaInt;
    if ( (*iter)->active_) {
      pepMassIntUnique.push_back(pepMaInt);
    }
    pe++;
  }

  //For pepMassIntUnique vector
  //Sort vector, take unique of vector, get rid of extra space in vector
  std::sort(pepMassIntUnique.begin(), pepMassIntUnique.end());
  vector<int>::iterator last = std::unique(pepMassIntUnique.begin(),
                                           pepMassIntUnique.end());
  pepMassIntUnique.erase(last, pepMassIntUnique.end());
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
    firstTerm += pow(lnpy, i) / double(MathUtil::factorial(i));
  }
  firstTerm = pow(p, y) * firstTerm;

  //compute second term (p-values are completely independent)
  double secondTerm = pow(p, y) * realPartofM * pow(lnpy, intPartofM) / (double)MathUtil::factorial(intPartofM);

  return firstTerm + secondTerm;
}

void TideSearchApplication::setSpectrumFlag(map<pair<string, unsigned int>, bool>* spectrum_flag) {
  spectrum_flag_ = spectrum_flag;
}