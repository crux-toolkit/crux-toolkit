/*********************************************************************//**
 * \file MatchCollection.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: 11/27 2006
 * \brief A set of peptide spectrum matches for one spectrum.
 *
 * Methods for creating and manipulating match_collections.   
 * Creating a match collection generates all matches (searches a
 * spectrum against a database.
 ****************************************************************************/
#include "MatchCollection.h"
#include "MatchCollectionIterator.h"
#include "MatchIterator.h"
#include <string>
#include "MatchFileReader.h"
#include "SQTReader.h"
#include "WinCrux.h"

using namespace std;
using namespace Crux;

/**
 * \returns An (empty) match_collection object.
 */
void MatchCollection::init() {

  match_total_ = 0;
  experiment_size_ = 0;
  target_experiment_size_ = 0;
  zstate_ = SpectrumZState();
  null_peptide_collection_ = false;

  // loop over to set all score type to false
  for(int score_type_idx = 0 ; 
    score_type_idx < NUMBER_SCORER_TYPES ; 
    ++score_type_idx){

    scored_type_[score_type_idx] = false;
  }

  // set last score to -1, thus nothing has been done yet
  last_sorted_ = (SCORER_TYPE_T)-1;
  iterator_lock_ = false;
  delta_cn_ = 0;
  sp_scores_sum_ = 0;
  sp_scores_mean_ = 0;
  mu_ = 0;
  l_value_ = 0;
  top_fit_sp_ = 0;
  base_score_sp_ = 0;
  eta_ = 0;
  beta_ = 0;
  shift_ = 0;
  correlation_ = 0;
  num_samples_ = 0;
  num_xcorrs_ = 0;

  post_process_collection_ = false;
  post_scored_type_set_ = false;
  top_scoring_sp_ = NULL;
  exact_pval_search_ = false;
  has_distinct_matches_ = false;
}

/**
 * /brief Free the memory allocated for a match collection
 * Deep free; each match is freed which, in turn, frees each spectrum
 * and peptide. 
 */
MatchCollection::~MatchCollection() {

  // decrement the pointer count in each match object
  while(match_total_ > 0){
    --match_total_;
    Match::freeMatch(match_[match_total_]);
    match_[match_total_] = NULL;
  }

  // and free the sample matches
  while(num_samples_ > 0){
    --num_samples_;
    Match::freeMatch(sample_matches_[num_samples_]);
    sample_matches_[num_samples_] = NULL;
  }

  if(top_scoring_sp_){
    Match::freeMatch(top_scoring_sp_);
  }
}

/**
 * \brief Creates a new match collection with no matches in it.  Sets
 * the member variable indicating if the matches are to real peptides
 * or to decoy (shuffled) peptides. Other member variables are set to
 * default values.  The method add_matches() can be used to search a
 * spectrum and store the matches in this collection.  
 *
 * \returns A newly allocated match collection
 */
MatchCollection::MatchCollection(
  bool is_decoy
  ){

  init();
  null_peptide_collection_ = is_decoy;
}

void MatchCollection::preparePostProcess() {

  // prepare the match_collection
  init();
  // set this as a post_process match collection
  post_process_collection_ = true;
}


/**
 * \brief Creates a new match_collection from the match collection
 * iterator. 
 *
 * Used in the post_processing extension.  Also used by
 * setup_match_collection_iterator which is called by next to find,
 * open, and parse the next psm file(s) to process.  If there are
 * multiple target psm files, it reads in all of them when set_type is
 * 0 and puts them all into one match_collection.  If the fileroot
 * parameter is non-null, only reads files with that prefix.
 *\returns A heap allocated match_collection.
 */
MatchCollection::MatchCollection(
  MatchCollectionIterator* match_collection_iterator, 
    ///< the working match_collection_iterator -in
  SET_TYPE_T set_type  
    ///< what set of match collection are we creating? (TARGET, DECOY1~3) -in 
  )
{ 

  Database* database = match_collection_iterator->getDatabase();
  Database* decoy_database = match_collection_iterator->getDecoyDatabase();

  preparePostProcess();

  // get the list of files to open
  vector<string> file_names;
  get_target_decoy_filenames(file_names, 
                             match_collection_iterator->getWorkingDirectory(),
                             set_type);

  // open each file and add psms to match collection
  for(int file_idx = 0; file_idx < (int)file_names.size(); file_idx++){
    char* full_filename = 
      get_full_filename(match_collection_iterator->getDirectoryName(),
                        file_names[file_idx].c_str());
    MatchFileReader delimited_result_file(full_filename);
    carp(CARP_DEBUG, "Creating new match collection from '%s' file.",
         full_filename);
    free(full_filename);

    extendTabDelimited(database, delimited_result_file, decoy_database);

    // for the first target file, set headers based on input files
    if( set_type == SET_TARGET && file_idx == 0 ){
      delimited_result_file.getMatchColumnsPresent(
                                  match_collection_iterator->getColsInFile()); 
    }
  } // next file
}

/**
 * \brief The main search function.  All peptides in the peptide
 * iterator are compared to the spectrum and the resulting score(s)
 * are stored in a match.  All matches are stored in the
 * match_collection.  Can be called on an empty match_collection or
 * one already containing matches.  No checks to confirm that the same
 * spectrum is being searched in subsequent calls.  Assumes that the
 * charge has been set for the match_collection.
 *
 * First, a match with no scores is generated for each peptide.  If
 * is_decoy is true, the peptides will be shuffled.  Then if
 * do_prelim_scoring is true, the matches are all scored with SP, 
 * ranked, and the collection truncated.
 *
 * Next, all (remaining) matches are scored with XCORR, the matches
 * are ranked, and the collection is truncated.
 *
 * If store_scores is true, each XCORR will be kept by the match
 * collection in a separate array so that it is retained even after
 * matches are deleted with truncation.  The scores can later be used
 * to estimate weibull parameters.
 * 
 * \returns The number of matches added.
 */
int MatchCollection::addMatches(
  Spectrum* spectrum,  ///< compare peptides to this spectrum
  SpectrumZState& zstate,            ///< use this charge state for spectrum
  ModifiedPeptidesIterator* peptide_iterator, ///< use these peptides
  bool is_decoy,     ///< are peptides to be shuffled
  bool store_scores, ///< save scores for p-val estimation
  bool do_sp_scoring, ///< start with SP scoring
  bool filter_by_sp  ///< truncate matches based on Sp scores
){

  if( peptide_iterator == NULL || spectrum == NULL ){
    carp(CARP_FATAL, "Cannot add matches to a collection when match " 
         "collection, spectrum and/or peptide iterator are NULL.");
  }

  //assert(matches->zstate == zstate);

  // generate a match for each peptide in the iterator, storing them
  // in the match collection
  int num_matches_added = addUnscoredPeptides(spectrum, zstate,
                                              peptide_iterator, is_decoy);

  if( num_matches_added == 0 ){
    if (do_sp_scoring) {
      scored_type_[SP] = true;
    }
    scored_type_[XCORR] = true;
    return num_matches_added;
  }

  int xcorr_max_rank = get_int_parameter("psms-per-spectrum-reported");

  // optional Sp score on all candidate peptides
  if( do_sp_scoring ){
    scoreMatchesOneSpectrum(SP, spectrum, zstate.getCharge(),
                               false); // don't store scores
    populateMatchRank(SP);
    saveTopSpMatch();
    if( filter_by_sp ){ // keep only high-ranking sp psms
      int sp_max_rank = get_int_parameter("max-rank-preliminary");
      truncate(sp_max_rank + 1, // extra for deltacn of last
               SP);
      xcorr_max_rank = sp_max_rank;
    }
  }

  // main scoring
  scoreMatchesOneSpectrum(XCORR, spectrum, zstate.getCharge(), 
                             store_scores); 
  populateMatchRank(XCORR);
  truncate(xcorr_max_rank + 1,// extra for deltacn of last
           XCORR);

  return num_matches_added;
}

/**
 * \brief Put all the matches from the source match collection in the
 * destination. Only copies the pointers of the matches so use with
 * caution. 
 * \returns The number of matches added.
 */
int MatchCollection::merge(MatchCollection* source,
                           MatchCollection* destination){
  if( source == NULL || destination == NULL ){
    carp(CARP_ERROR, "Cannot merge null match collections.");
    return 0;
  }
  carp(CARP_DETAILED_DEBUG, "Merging match collections.");

  // what is the index of the next insert position in destination
  int dest_idx = destination->match_total_;

  // if these are the first being added to the destination, set the
  // scored_type
  if( dest_idx == 0 ){
    int type_idx = 0;
    for(type_idx = 0; type_idx < NUMBER_SCORER_TYPES; type_idx++){
      destination->scored_type_[type_idx] = source->scored_type_[type_idx];
    }
  }else{ // check that same types are scored
    int type_idx = 0;
    for(type_idx = 0; type_idx < NUMBER_SCORER_TYPES; type_idx++){
      if( destination->scored_type_[type_idx] != source->scored_type_[type_idx]){
        const char* dest_str = (destination->scored_type_[type_idx]) ? "" : " not";
        const char* src_str = (source->scored_type_[type_idx]) ? "" : " not";
        const char* type_str = scorer_type_to_string((SCORER_TYPE_T)type_idx);
        carp(CARP_FATAL, "Cannot merge match collections scored for "
             "different types.  Trying to add matches%s scored for %s "
             "to matches%s scored for %s", 
             src_str, type_str, dest_str, type_str);
      }
    }
  }
  

  // make sure destination has room for more matches
  int src_num_matches = source->match_total_;
  if( dest_idx + src_num_matches > _MAX_NUMBER_PEPTIDES ){
    carp(CARP_FATAL, "Cannot merge match collections, insufficient capacity "
         "in destnation collection.");
  }
  carp(CARP_DETAILED_DEBUG, "Merging %d matches into a collection of %d",
       src_num_matches, dest_idx );

  int src_idx = 0;
  // for each match in source
  for(src_idx = 0; src_idx < src_num_matches; src_idx++){
    Match* cur_match = source->match_[src_idx];

    // copy pointer and add to destination
    cur_match->incrementPointerCount();
    destination->match_[dest_idx] = cur_match;

    dest_idx++;
  }

  // update destination count
  // a target match collection may not have the target experiment size set
  if( destination->target_experiment_size_ == 0 ){ 
    destination->target_experiment_size_ = 
      destination->getTargetExperimentSize();
  }
  destination->match_total_ += src_num_matches;
  destination->experiment_size_ += source->experiment_size_;
  destination->target_experiment_size_ += source->getTargetExperimentSize();
  destination->last_sorted_ = (SCORER_TYPE_T)-1;  // unset any last-sorted flag

  return src_num_matches;
}


/**
 * \brief Store the xcorr for each psm that was added in this
 * iteration.  Assumes that the matches with scores needing storing
 * are between indexes start_index and match_collection->match_total.
 * The xcorrs will used for the weibull parameter estimations for
 * p-values.  If keep_matches == false, the matches between indexes
 * start_index and match_collection->match_total will be deleted and
 * match_total will be updated.
 * 
 */
void MatchCollection::storeNewXcorrs(
  int start_index, ///< get first score from match at this index
  bool keep_matches ///< false=delete the matches after storing score
){

  int score_idx = num_xcorrs_;
  int psm_idx = start_index;

  carp(CARP_DETAILED_DEBUG, 
       "Adding to xcors[%i] scores from psm index %i to %i", 
       score_idx, psm_idx, match_total_);

  if( score_idx+(match_total_-psm_idx) 
      > _MAX_NUMBER_PEPTIDES ){
    carp(CARP_FATAL, "Too many xcorrs to store.");
  }

  for(psm_idx=start_index; psm_idx < match_total_; psm_idx++){
    FLOAT_T score = match_[psm_idx]->getScore(XCORR);
    xcorrs_[score_idx] = score;
    score_idx++;

    if( keep_matches == false ){
      Match::freeMatch(match_[psm_idx]);
      match_[psm_idx] = NULL;
      experiment_size_ -= 1;  // these should be decoys and 
                                               // we are not counting them
                                               
    }
  }

  num_xcorrs_ = score_idx;
  if( keep_matches == false ){
    match_total_ = start_index; // where we started deleting
  }
  carp(CARP_DETAILED_DEBUG, "There are now %i xcorrs.", score_idx);
}


/**
 * \brief After psms have been added to a match collection but before
 * the collection has been truncated, go through the list of matches
 * and combine those that are for the same peptide sequence.
 *
 * Requires that the match_collection was sorted by Sp so that
 * matches with identical peptides will be listed together.
 */
void MatchCollection::collapseRedundantMatches(){

  // must not be empty
  int match_total = match_total_;
  if( match_total == 0 ){
    return;
  }  

  carp(CARP_DETAILED_DEBUG, "Collapsing %i redundant matches.", match_total_);

  // must be sorted by Sp or xcorr
  assert( (last_sorted_ == SP) || 
          (last_sorted_ == XCORR) );

  Match** matches = match_;
  int match_idx = 0;
  FLOAT_T cur_score = matches[match_idx]->getScore(SP);

  // for entire list of matches
  while(match_idx < match_total-1){
    FLOAT_T next_score = matches[match_idx+1]->getScore(SP);

    // find the index of the last match with the same score
    int cur_score_last_index = match_idx;
    
    while(next_score == cur_score && cur_score_last_index < match_total-2){
      cur_score_last_index++;
      next_score = matches[cur_score_last_index+1]->getScore(SP);
    }
    // if the last two were equal, the last index was not incremented
    if( next_score == cur_score ){ cur_score_last_index++; }

    if( cur_score_last_index > match_idx ){
      consolidateMatches(matches, match_idx, cur_score_last_index);
    }

    match_idx = cur_score_last_index+1;
    cur_score = next_score;
  }// next match

  // shift contents of the match array to fill in deleted matches
  int opening_idx = 0;
  while( matches[opening_idx] != NULL && opening_idx < match_total){
    opening_idx++;
  }

  for(match_idx=opening_idx; match_idx<match_total; match_idx++){
    if( matches[match_idx] != NULL ){ // then move to opening
      matches[opening_idx] = matches[match_idx];
      opening_idx++;
    }
  }

  carp(CARP_DETAILED_DEBUG, "Removing duplicates changed count from %i to %i",
       match_total_, opening_idx);
  // reset total number of matches in the collection
  match_total_ = opening_idx;
  // remove duplicate peptides from the overall count
  int diff = match_total - opening_idx;
  carp(CARP_DETAILED_DEBUG, "Removing %i from total count %i",
       diff, experiment_size_);

  experiment_size_ -= diff;
}

/**
 * \brief For a list of matches with the same scores, combine those
 * that are the same peptide and delete redundant matches.
 *
 * Since there may be different peptide sequences with the same score,
 * compare each match to the remaining matches.
 */
void MatchCollection::consolidateMatches(
  Match** matches, 
  int start_idx, 
  int end_idx
  ){

  carp(CARP_DETAILED_DEBUG, "Consolidating index %i to %i.", start_idx, end_idx);
  int cur_match_idx = 0;
  for(cur_match_idx=start_idx; cur_match_idx < end_idx; cur_match_idx++){
    carp(CARP_DETAILED_DEBUG, "Try consolidating with match[%i].", 
         cur_match_idx);

    if(matches[cur_match_idx] == NULL){
      carp(CARP_DETAILED_DEBUG, "Can't consolodate with %i, it's null.", 
           cur_match_idx);
      continue;
    }    

    char* cur_seq = 
      matches[cur_match_idx]->getModSequenceStrWithSymbols();
    carp(CARP_DETAILED_DEBUG, "cur seq is %s.", cur_seq);
    int next_match_idx = cur_match_idx+1;
    for(next_match_idx=cur_match_idx+1; next_match_idx<end_idx+1; 
        next_match_idx++){
      carp(CARP_DETAILED_DEBUG, "Can match[%i] be added to cur.", 
           next_match_idx);

      if(matches[next_match_idx] == NULL){
        continue;
      }    

      char* next_seq = 
        matches[next_match_idx]->getModSequenceStrWithSymbols();
      carp(CARP_DETAILED_DEBUG, "next seq is %s.", next_seq);

      if( strcmp(cur_seq, next_seq) == 0){
        carp(CARP_DETAILED_DEBUG, 
             "Seqs %s and %s match.  Consolidate match[%i] into match[%i].", 
             cur_seq, next_seq, next_match_idx, cur_match_idx);

        // add peptide src of next to cur
        Peptide::mergePeptidesCopySrc( matches[cur_match_idx]->getPeptide(),
                        matches[next_match_idx]->getPeptide());
        // this frees the second peptide, so set what pointed to it to NULL
        //set_match_peptide(matches[next_match_idx], NULL);

        // delete match
        Match::freeMatch(matches[next_match_idx]);
        matches[next_match_idx] = NULL;
      }

      free(next_seq);
    }// next match to delete

    free(cur_seq);
  }// next match to consolidate to
}


/**
 * Sort the match collection by score type.
 */
void MatchCollection::sort(
  SCORER_TYPE_T score_type ///< the score type to sort by -in
  )
{
  carp(CARP_DEBUG, "Sorting match collection.");

  // check if we are allowed to alter match_collection
  if(iterator_lock_){
    carp(CARP_FATAL,
         "Cannot sort a match collection when a match iterator is already"
         " instantiated");
  }

  // Switch to the equivalent sort key.
  SCORER_TYPE_T sort_by = NUMBER_SCORER_TYPES; // Initialize to nonsense.
  int (*compare_match_function)(const void*, const void*) 
    = (QSORT_COMPARE_METHOD)compareSp;
  switch (score_type) {
  case SP: 
    carp(CARP_DEBUG, "Sorting match collection by Sp.");
    sort_by = SP;
    compare_match_function = (QSORT_COMPARE_METHOD)compareSp;
    break;
  case EVALUE:
    sort_by = EVALUE;
    compare_match_function = (QSORT_COMPARE_METHOD)compareEValue;
    break;
  case XCORR:
  case DECOY_XCORR_QVALUE:
  case LOGP_WEIBULL_XCORR: 
  case DECOY_XCORR_PEPTIDE_QVALUE:
  case DECOY_XCORR_PEP:
    carp(CARP_DEBUG, "Sorting match collection by XCorr.");
    sort_by = XCORR;
    compare_match_function = (QSORT_COMPARE_METHOD)compareXcorr;
    break;

  case LOGP_BONF_WEIBULL_XCORR: 
  case LOGP_QVALUE_WEIBULL_XCORR:
  case LOGP_PEPTIDE_QVALUE_WEIBULL:
  case LOGP_WEIBULL_PEP:
    carp(CARP_DEBUG, "Sorting match collection by p-value.");
    sort_by = LOGP_BONF_WEIBULL_XCORR;
    compare_match_function = (QSORT_COMPARE_METHOD)comparePValue;
    break;
/*
  case PERCOLATOR_SCORE:
    carp(CARP_INFO, "Sorting match collection by Percolator score.");
    sort_by = PERCOLATOR_SCORE;
    compare_match_function = (QSORT_COMPARE_METHOD)comparePercolatorScore;
    break;
*/
  case PERCOLATOR_QVALUE:
  case PERCOLATOR_PEPTIDE_QVALUE:
  case PERCOLATOR_PEP:
    carp(CARP_DEBUG, "Sorting match collection by Percolator q-value.");
    sort_by = PERCOLATOR_QVALUE;
    compare_match_function = (QSORT_COMPARE_METHOD)comparePercolatorQValue;
    break;
/*
  case PERCOLATOR_SCAN:
    carp(CARP_INFO, "Sorting match collection by Percolator scan.");
    sort_by = PERCOLATOR_SCAN;
    compare_match_function = (QSORT_COMPARE_METHOD)compareSpectrum;
    break;
*/
  case QRANKER_SCORE:
    carp(CARP_DEBUG, "Sorting match collection by Q-ranker score.");
    sort_by = QRANKER_SCORE;
    compare_match_function = (QSORT_COMPARE_METHOD)compareQRankerScore;
    break;

  case QRANKER_QVALUE:
  case QRANKER_PEPTIDE_QVALUE:
  case QRANKER_PEP:
    carp(CARP_DEBUG, "Sorting match collection by Q-ranker q-value.");
    compare_match_function = (QSORT_COMPARE_METHOD)compareQRankerQValue;
    sort_by = QRANKER_QVALUE;
    break;

  case BARISTA_SCORE:
    carp(CARP_DEBUG, "Sorting match collection by barista score.");
    sort_by = BARISTA_SCORE;
    compare_match_function = (QSORT_COMPARE_METHOD)compareBaristaScore;
    break;

  case BARISTA_QVALUE:
  case BARISTA_PEPTIDE_QVALUE:
  case BARISTA_PEP:
    carp(CARP_DEBUG, "Sorting match collection by barista q-value.");
    compare_match_function = (QSORT_COMPARE_METHOD)compareBaristaQValue;
    sort_by = BARISTA_QVALUE;
    break;

  // Should never reach this point.
  case NUMBER_SCORER_TYPES:
  case INVALID_SCORER_TYPE:
    carp(CARP_FATAL, "Something is terribly wrong in the sorting code!");
  }

  // Don't sort if it's already sorted.
  if (last_sorted_ == sort_by) {
    return;
  }

  // Do the sort.
  qsortMatch(match_,
             match_total_,
             compare_match_function);
  last_sorted_ = sort_by;
}

/**
 * \brief Sort a match_collection by the given score type, grouping
 * matches by spectrum (if multiple spectra are present).
 */
void MatchCollection::spectrumSort(
  SCORER_TYPE_T score_type ///< the score type to sort by -in
  ){

  // check if we are allowed to alter match_collection
  if(iterator_lock_){
    carp(CARP_FATAL,
         "Cannot alter match_collection when a match iterator is already"
         " instantiated");
  }

  switch(score_type){
  case SP: 
    qsortMatch(match_, match_total_,
                (QSORT_COMPARE_METHOD)compareSpectrumSp);
    last_sorted_ = SP;
    break;

  case XCORR:
  case LOGP_WEIBULL_XCORR: 
  case LOGP_BONF_WEIBULL_XCORR: 
  case LOGP_QVALUE_WEIBULL_XCORR: 
  case LOGP_PEPTIDE_QVALUE_WEIBULL:
  case LOGP_WEIBULL_PEP:
  case DECOY_XCORR_QVALUE:
  case DECOY_XCORR_PEPTIDE_QVALUE:
  case DECOY_XCORR_PEP:
    /* If we are sorting on a per-spectrum basis, then the xcorr is
       good enough, even in the presence of p-values. */
    qsortMatch(match_, match_total_,
                (QSORT_COMPARE_METHOD)compareSpectrumXcorr);
    last_sorted_ = XCORR;
    break;

  case PERCOLATOR_SCORE:
    qsortMatch(match_, match_total_,
                (QSORT_COMPARE_METHOD)compareSpectrumPercolatorScore);
    last_sorted_ = PERCOLATOR_SCORE;
    break;

  case PERCOLATOR_QVALUE:
  case PERCOLATOR_PEPTIDE_QVALUE:
  case PERCOLATOR_PEP:
    qsortMatch(match_, match_total_,
                (QSORT_COMPARE_METHOD)compareSpectrumPercolatorQValue);
    last_sorted_ = PERCOLATOR_QVALUE;
    break;
/*
  case PERCOLATOR_SCAN: 
    qsortMatch(
      match_,
      match_total_,
      (QSORT_COMPARE_METHOD)compareSpectrumScan
    );
    last_sorted_=PERCOLATOR_SCAN;
    break; 
*/
  case QRANKER_SCORE:
    qsortMatch(match_, match_total_,
                (QSORT_COMPARE_METHOD)compareSpectrumQRankerScore);
    last_sorted_ = QRANKER_SCORE;
    break;

  case QRANKER_QVALUE:
  case QRANKER_PEPTIDE_QVALUE:
  case QRANKER_PEP:
    qsortMatch(match_, match_total_,
                (QSORT_COMPARE_METHOD)compareSpectrumQRankerQValue);
    last_sorted_ = QRANKER_QVALUE;
    break;

  case BARISTA_SCORE:
    qsortMatch(match_, match_total_,
                (QSORT_COMPARE_METHOD)compareSpectrumBaristaScore);
    last_sorted_ = BARISTA_SCORE;
    break;

  case BARISTA_QVALUE:
  case BARISTA_PEPTIDE_QVALUE:
  case BARISTA_PEP:
    qsortMatch(match_, match_total_,
                (QSORT_COMPARE_METHOD)compareSpectrumBaristaQValue);
    last_sorted_ = BARISTA_QVALUE;
    break;


  // Should never reach this point.
  case NUMBER_SCORER_TYPES:
  case INVALID_SCORER_TYPE:
    carp(CARP_FATAL, "Something is terribly wrong in the sorting code!");
 }
}

/**
 * \brief Reduces the number of matches in the match_collection so
 * that only the [max_rank] highest scoring (by score_type) remain.
 *
 * Matches ranking up to max_rank are retained and those ranking
 * higher are freed.  The value of match_collection->total_matches is
 * adjusted to reflect the remaining number of matches.  The max rank
 * and total_matches may not be the same value if there are multiple
 * matches with the same rank.  Sorts match collection by score_type,
 * if necessary.
 */
void MatchCollection::truncate(
  int max_rank,     ///< max rank of matches to keep from SP -in
  SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
  )
{
  carp(CARP_DETAILED_DEBUG, "Truncating match collection to rank %d.", max_rank);
  if ( match_total_ == 0){
    carp(CARP_DETAILED_DEBUG, "No matches in collection, so not truncating");
    return;
  }

  // sort match collection by score type
  sort(score_type);

  // Free high ranking matches
  int highest_index = match_total_ - 1;
  Match** matches = match_;
  int cur_rank = matches[highest_index]->getRank(score_type);

  while( cur_rank > max_rank ){
    Match::freeMatch(matches[highest_index]);
    highest_index--;
    cur_rank = matches[highest_index]->getRank(score_type);
  }
  match_total_ = highest_index + 1;

  carp(CARP_DETAILED_DEBUG, "Truncated collection now has %d matches.", 
       match_total_);

}

/**
 * Assigns a rank for the given score type to each match.  First sorts
 * by the score type (if not already sorted).  Overwrites any existing
 * rank values, so it can be performed on a collection with matches
 * newly added to previously ranked matches.  Rank 1 is highest
 * score.  Matches with the same score will be given the same rank.
 *
 * \returns true, if populates the match rank in the match collection
 */
bool MatchCollection::populateMatchRank(
 SCORER_TYPE_T score_type ///< score type (SP, XCORR) by which to rank -in
 )
{
  carp(CARP_DETAILED_DEBUG, "Ranking matches by %i.", score_type);
  carp(CARP_DETAILED_DEBUG, "Collection currently ranked by %d", last_sorted_);
  // check if the match collection is in the correct sorted order
  if(last_sorted_ != score_type){
    // sort match collection by score type
    carp(CARP_DETAILED_DEBUG, "Sorting by score_type %i", score_type);
    sort(score_type);
  }

  // set match rank for all match objects that have been scored for
  // this type
  int match_index;
  int cur_rank = 0;
  FLOAT_T cur_score = NOT_SCORED;
  for(match_index=0; match_index<match_total_; ++match_index){
    Match* cur_match = match_[match_index];
    FLOAT_T this_score = cur_match->getScore(score_type);
    
    if( NOT_SCORED == cur_match->getScore(score_type) ){
      char* seq = cur_match->getModSequenceStrWithMasses(MOD_MASS_ONLY);
      carp(CARP_WARNING, 
           "PSM spectrum %i charge %i sequence %s was NOT scored for type %i",
           cur_match->getSpectrum()->getFirstScan(),
           cur_match->getCharge(), seq,
           (int)score_type);
      free(seq);
    }

    // does this match have a higher score?
    if( this_score != cur_score ){
      cur_score = this_score;
      cur_rank++;
    }

    //    set_match_rank( cur_match, score_type, match_index+1);
    cur_match->setRank(score_type, cur_rank);

    carp(CARP_DETAILED_DEBUG, "Match rank %i, score %f", cur_rank, cur_score);
  }
  
  return true;
}

/**
 * Keep track of the top-scoring Sp match.  It should be printed to
 * the sqt file even if its XCORR rank is not high enough to be
 * printed.  Requires that ranks have been set for Sp.
 *
 */
void MatchCollection::saveTopSpMatch(){

  assert(match_total_ > 0);
  Match* cur_rank_one_match = match_[0];

  // confirm that matches are sorted and ranks are set
  if( last_sorted_ != SP || 
      cur_rank_one_match->getRank(SP) != 1 ){
    populateMatchRank(SP);
  }

  // if no top sp yet, set it
  if( top_scoring_sp_ == NULL ){
    top_scoring_sp_ = cur_rank_one_match;
    cur_rank_one_match->incrementPointerCount();
    return;
  }

  // otherwise, see if the current top-ranked match has a higher score
  // the rank of top_scoring_sp should have a new rank
  if( top_scoring_sp_->getRank(SP) > 1 ){
    Match::freeMatch(top_scoring_sp_);
    top_scoring_sp_ = cur_rank_one_match;
    cur_rank_one_match->incrementPointerCount();
  }
}

/**
 * Create a new match_collection by randomly sampling matches 
 * from match_collection upto count number of matches
 * Must not free the matches
 * \returns a new match_collection of randomly sampled matches 
 */
MatchCollection* MatchCollection::randomSample(
  int count_max ///< the number of matches to randomly select -in
  )
{
  int count_idx = 0;
  int match_idx = 0;
  int score_type_idx = 0;
  
  MatchCollection* sample_collection = new MatchCollection();
  mysrandom(time(NULL));

  // make sure we don't sample more than the matches in the match collection
  if (count_max >= match_total_){
    delete sample_collection;
    return this;
  }

  // ranomly select matches upto count_max
  for(; count_idx < count_max; ++count_idx){
    match_idx =
      (int)(
        (double)myrandom()/((double)UNIFORM_INT_DISTRIBUTION_MAX + (double)1)
      ) * match_total_;
    
    // match_idx = random() % match_collection->match_total;
    sample_collection->match_[count_idx] = match_[match_idx];
    // increment pointer count of the match object 
    sample_collection->match_[count_idx]->incrementPointerCount();
  }
  
  // set total number of matches sampled
  sample_collection->match_total_ = count_idx;

  sample_collection->experiment_size_ = experiment_size_;
  sample_collection->target_experiment_size_ = target_experiment_size_;

  // set scored types in the sampled matches
  for(; score_type_idx < NUMBER_SCORER_TYPES ;  ++score_type_idx){
    sample_collection->scored_type_[score_type_idx] 
      = scored_type_[score_type_idx];
  }
  
  return sample_collection;
}

/**
 * This function is a transformation of the partial derivatives of
 * the log likelihood of the data given an extreme value distribution
 * with location parameter mu and scale parameter 1/L. The transformation 
 * has eliminated the explicit dependence on the location parameter, mu, 
 * leaving only the scale parameter, 1/L.
 *
 * The zero crossing of this function will correspond to the maximum of the 
 * log likelihood for the data.
 *
 * See equations 10 and 11 of "Maximum Likelihood fitting of extreme value 
 * distributions".
 *
 * The parameter values contains a list of the data values.
 * The parameter L is the reciprocal of the scale parameters.
 *
 *\returns the final exponential values of the score and sets the value of the function and its derivative.
 */
void MatchCollection::constraintFunction(
  SCORER_TYPE_T score_type, ///< score_type to estimate EVD distribution -in
  FLOAT_T l_value,  ///< L value -in
  FLOAT_T* function,  ///< the output function value -out
  FLOAT_T* derivative,  ///< the output derivative value -out
  FLOAT_T* exponential_sum ///< the final exponential array sum -out
  )
{
  int idx = 0;
  FLOAT_T* exponential = (FLOAT_T*)mycalloc(match_total_, sizeof(FLOAT_T));
  FLOAT_T numerator = 0;
  FLOAT_T second_numerator = 0;
  FLOAT_T score = 0;
  FLOAT_T denominator = 0;
  FLOAT_T score_sum = 0;
  Match** matches = match_;

  // iterate over the matches to calculate numerator, exponential value, denominator
  for(; idx < match_total_; ++idx){
    score = matches[idx]->getScore(score_type);
    exponential[idx] = exp(-l_value * score);
    numerator += (exponential[idx] * score);
    denominator += exponential[idx];
    score_sum += score;
    second_numerator += (score * score * exponential[idx]);
  }

  // assign function value
  *function = (1.0 / l_value) - (score_sum / match_total_) 
    + (numerator / denominator);

  // assign derivative value
  *derivative =  ((numerator * numerator) / (denominator * denominator)) 
    - ((second_numerator / denominator)) - (1.0 / (l_value * l_value));

  // assign the total sum of the exponential values
  *exponential_sum = denominator;

  // free exponential array
  free(exponential);
}


/**
 * For the #top_count ranked peptides, calculate the Weibull parameters
 *\returns true, if successfully calculates the Weibull parameters
 */
static const FLOAT_T MIN_WEIBULL_MATCHES = 40;
static const FLOAT_T MIN_XCORR_SHIFT = -5.0;
static const FLOAT_T MAX_XCORR_SHIFT  = 5.0;
//#define CORR_THRESHOLD 0.995   // Must achieve this correlation, else punt.
static const FLOAT_T CORR_THRESHOLD = 0.0;       // For now, turn off the threshold.
static const FLOAT_T XCORR_SHIFT = 0.05;
static const FLOAT_T MIN_SP_SHIFT = -100.0;
static const FLOAT_T MAX_SP_SHIFT = 300.0;
static const FLOAT_T SP_SHIFT = 5.0;

/**
 * \brief Check that a match collection has a sufficient number of
 * matches for estimating Weibull parameters.
 * \returns true if their are enough xcorrs for estimating Weibull
 * parameters or false if not.
 */
bool MatchCollection::hasEnoughWeibullPoints()
{
  return (num_xcorrs_ >= MIN_WEIBULL_MATCHES );
}

/**
 * \brief Use the xcorrs saved in the match_collection to estimate the
 * weibull parameters to be used for computing p-values. 
 *
 * Requires that main score be XCORR, but with relatively few changes
 * other scores could be accomodated.
 * Implementation of Weibull distribution parameter estimation from 
 * http:// www.chinarel.com/onlincebook/LifeDataWeb/rank_regression_on_y.htm
 */
bool MatchCollection::estimateWeibullParametersFromXcorrs(
  Spectrum* spectrum,
  int charge
  ){

  if( spectrum == NULL ){
    carp(CARP_ERROR, "Cannot estimate parameters from null inputs.");
    return false;
  }

  // check that we have the minimum number of matches
  FLOAT_T* scores = xcorrs_;
  int num_scores = num_xcorrs_;
  if( num_scores < MIN_WEIBULL_MATCHES ){
    carp(CARP_DETAILED_DEBUG, "Too few psms (%i) to estimate "
         "p-value parameters for spectrum %i, charge %i",
         num_scores, spectrum->getFirstScan(), charge);
    // set eta, beta, and shift to something???
    return false;
  }

  // reverse sort the scores
  std::sort(scores, scores + num_scores, greater<FLOAT_T>());

  // use only a fraction of the samples, the high-scoring tail
  // this parameter is hidden from the user
  double fraction_to_fit = get_double_parameter("fraction-top-scores-to-fit");
  assert( fraction_to_fit >= 0 && fraction_to_fit <= 1 );
  int num_tail_samples = (int)(num_scores * fraction_to_fit);
  carp(CARP_DEBUG, "Estimating Weibull params with %d psms (%.2f of %i)", 
       num_tail_samples, fraction_to_fit, num_scores);

  // do the estimation
  fit_three_parameter_weibull(scores, num_tail_samples, num_scores,
      MIN_XCORR_SHIFT, MAX_XCORR_SHIFT, XCORR_SHIFT, CORR_THRESHOLD,
      &(eta_), &(beta_),
      &(shift_), &(correlation_));
  carp(CARP_DEBUG, 
      "Corr: %.6f  Eta: %.6f  Beta: %.6f  Shift: %.6f", 
       correlation_, eta_, beta_, shift_);
  
  return true;
}

// TODO (BF 16-mar-09): use this instead of score_peptides
/**
 * \brief Add all peptides from iterator to match collection.
 * Additional matches will not be scored for any type.
 * \returns The number of peptides added.
 */
int MatchCollection::addUnscoredPeptides(
  Spectrum* spectrum, 
  SpectrumZState& zstate, 
  ModifiedPeptidesIterator* peptide_iterator,
  bool is_decoy
){

  if( spectrum == NULL || peptide_iterator == NULL ){
    carp(CARP_FATAL, "Cannot score peptides with NULL inputs.");
  }
  carp(CARP_DETAILED_DEBUG, "Adding decoy peptides to match collection? %i", 
       is_decoy);

  int starting_number_of_psms = match_total_;

  while( peptide_iterator->hasNext() ){
    // get peptide
    Peptide* peptide = peptide_iterator->next();

    // create a match
    Match* match = new Match(peptide, spectrum, zstate, is_decoy);

    // add to match collection
    if(match_total_ >= _MAX_NUMBER_PEPTIDES){
      carp(CARP_ERROR, "peptide count of %i exceeds max match limit: %d", 
          match_total_, _MAX_NUMBER_PEPTIDES);

      return false;
    }

    match_[match_total_] = match;
    match_total_++;

  }// next peptide

  int matches_added = match_total_ - starting_number_of_psms;
  experiment_size_ += matches_added;

  // Set ln experiment size for matches
  FLOAT_T ln_experiment_size = logf((FLOAT_T)getTargetExperimentSize());
  for (int i = 0; i < match_total_; i++) {
    match_[i]->setLnExperimentSize(ln_experiment_size);  
  }

  // matches are no longer correctly sorted
  last_sorted_ = (SCORER_TYPE_T)-1; // unsorted
  return matches_added;
}
/**
 * \brief Use the score type to compare the spectrum and peptide in
 * the matches in match collection.  
 *
 * If the match has already been scored for this type, it is not
 * scored at again.  Requires that the given spectrum  and charge
 * state are the same as the spectrum and charge state in each of the
 * matches.  
 *
 * \returns true, if matches are successfully scored.
 */
bool MatchCollection::scoreMatchesOneSpectrum(
  SCORER_TYPE_T score_type, 
  Spectrum* spectrum,
  int charge,
  bool store_scores
  ){

  if( spectrum == NULL ){
    carp(CARP_ERROR, "Cannot score matches in a NULL match collection.");
    return false;
  }
  
  Match** matches = match_;
  int num_matches = match_total_;

  carp(CARP_DETAILED_DEBUG, "Scoring matches for %s", 
       scorer_type_to_string(score_type));

  // create ion constraint
  IonConstraint* ion_constraint = 
    IonConstraint::newIonConstraintSmart(score_type, charge);

  // create scorer
  Scorer* scorer = new Scorer(score_type);

  // create a generic ion_series that will be reused for each peptide sequence
  IonSeries* ion_series = new IonSeries(ion_constraint, charge);  
  
  // score all matches
  int match_idx;

  for(match_idx = 0; match_idx < num_matches; match_idx++){

    Match* match = matches[match_idx];
    assert( match != NULL );

    // skip it if it's already been scored
    if( NOT_SCORED != match->getScore(score_type)){
      continue;
    }

    // make sure it's the same spec and charge
    assert( spectrum == match->getSpectrum());
    assert( charge == match->getCharge());
    char* sequence = match->getSequence();
    MODIFIED_AA_T* modified_sequence = match->getModSequence();

    // create ion series for this peptide
    ion_series->update(sequence, modified_sequence);
    ion_series->predictIons();

    // get the score
    FLOAT_T score = scorer->scoreSpectrumVIonSeries(spectrum, ion_series);

    // set score in match
    match->setScore(score_type, score);
    if( score_type == SP ){
      match->setBYIonInfo(scorer);
    }

    // save score in collection
    if( store_scores ){
      xcorrs_[num_xcorrs_] = score;
      num_xcorrs_++; 
    }

    IF_CARP_DETAILED_DEBUG(
      char* mod_seq = 
      modified_aa_string_to_string_with_masses(modified_sequence,
                                               strlen(sequence),
                                               MOD_MASS_ONLY);
      carp(CARP_DETAILED_DEBUG, "Second score %f for %s (null:%i)",
           score, mod_seq, match->getNullPeptide());
      free(mod_seq);
    )
    free(sequence);
    free(modified_sequence);
  }// next match

  // set the match_collection as having been scored
  scored_type_[score_type] = true;

  // clean up
  IonConstraint::free(ion_constraint);
  delete ion_series;
  delete scorer;
  return true;
}

/**
 * \brief  Uses the Weibull parameters estimated by
 * estimate_weibull_parameters() to compute a p-value for each psm in
 * the collection.
 *
 * Computes the p-value for score-type set in parameter.c (which should
 * have been used for estimating the parameters).  Stores scores at
 * match->match_scores[LOGP_BONF_WEIBULL_XCORR].  This function
 * previous performed in score_collection_logp_bonf_weibull_[xcorr,sp]. 
 * 
 * \returns true if p-values successfully computed for all matches,
 * else false.
 */
// FIXME (BF 8-Dec-2008): create new score-type P_VALUE to replace LOG...XCORR
bool MatchCollection::computePValues(
  FILE* output_pvalue_file ///< If non-NULL, file for storing p-values -in
  ){

  int scan_number = 
    match_[0]->getSpectrum()->getFirstScan();
  carp(CARP_DEBUG, "Computing p-values for %s spec %d charge %d "
       "with eta %f beta %f shift %f",
       (null_peptide_collection_) ? "decoy" : "target",
       scan_number,
       zstate_.getCharge(),
       eta_, beta_, shift_);

  SCORER_TYPE_T main_score = get_scorer_type_parameter("score-type");

  // check that the matches have been scored
  if(!scored_type_[main_score]){
    const char* type_str = scorer_type_to_string(main_score);
    carp(CARP_FATAL, 
         "Match collection was not scored by %s prior to computing p-values.",
         type_str);
  }

  // Print separator in the decoy p-value file.
  if (output_pvalue_file) {
    fprintf(output_pvalue_file, "# scan: %d charge: %d candidates: %d\n", 
            scan_number, zstate_.getCharge(),
            experiment_size_);
    fprintf(output_pvalue_file, 
            "# eta: %g beta: %g shift: %g correlation: %g\n",
            eta_, 
            beta_,
            shift_,
            correlation_);
  }

  // iterate over all matches 
  int match_idx =0;
  for(match_idx=0; match_idx < match_total_; match_idx++){
    Match* cur_match = match_[match_idx];

    // Get the Weibull p-value.
    double pvalue = compute_weibull_pvalue(cur_match->getScore(main_score),
                                           eta_, beta_, shift_);

    // Print the pvalue, if requested
    if (output_pvalue_file) {
      fprintf(output_pvalue_file, "%g\n", pvalue);
    }

    // Apply the Bonferroni correction.
    pvalue = bonferroni_correction(pvalue, experiment_size_);

    // set pvalue in match
    cur_match->setScore(LOGP_BONF_WEIBULL_XCORR, -log(pvalue));
    //#endif

  }// next match

  carp(CARP_DETAILED_DEBUG, "Computed p-values for %d PSMs.", match_idx);
  populateMatchRank(XCORR);

  // mark p-values as having been scored
  scored_type_[LOGP_BONF_WEIBULL_XCORR] = true;
  return true;
}

/**
 * \brief Use the matches collected from all spectra to compute FDR
 * and q_values from the ranked list of target and decoy scores.
 * Assumes the match_collection has an appropriate number of
 * target/decoy matches per spectrum (e.g. one target and one decoy
 * per spec).
 * \returns true if q-values successfully computed, else false.
 */
bool MatchCollection::computeDecoyQValues(){

  carp(CARP_DEBUG, "Computing decoy q-values.");

  // sort by score
  sort(XCORR);

  // compute FDR from a running total of number targets/decoys
  // FDR = #decoys / #targets
  FLOAT_T num_targets = 0;
  FLOAT_T num_decoys = 0;
  int match_idx = 0;
  for(match_idx = 0; match_idx < match_total_; match_idx++){
    Match* cur_match = match_[match_idx];

    if ( cur_match->getNullPeptide() == true ){
      num_decoys += 1;
    }else{
      num_targets += 1;
    }
    FLOAT_T score = num_decoys/num_targets;
    if( num_targets == 0 ){ score = 1.0; }

    /*
    if (peptide_level) {
      // Skip PSMs that are not top-scoring for their peptide.
      if (!is_peptide_level(cur_match)) {
        score = NOT_SCORED;
        set_match_score(cur_match, DECOY_XCORR_PEPTIDE_QVALUE, NOT_SCORED);
      } else {
        set_match_score(cur_match, DECOY_XCORR_PEPTIDE_QVALUE, SCORE);
      }
    } else {
    */
    cur_match->setScore(DECOY_XCORR_QVALUE, score);
    carp(CARP_DETAILED_DEBUG, 
         "match %i xcorr or pval %f num targets %i, num decoys %i, score %f",
         match_idx, cur_match->getScore(XCORR), 
         (int)num_targets, (int)num_decoys, score);
  }

  // compute q-value: go through list in reverse and use min FDR seen
  FLOAT_T min_fdr = 1.0;
  for(match_idx = match_total_-1; match_idx >= 0; match_idx--){
    Match* cur_match = match_[match_idx];
    FLOAT_T cur_fdr = cur_match->getScore(DECOY_XCORR_QVALUE);
    if( cur_fdr == P_VALUE_NA ){ continue; }

    if( cur_fdr < min_fdr ){
      min_fdr = cur_fdr;
    }

    cur_match->setScore(DECOY_XCORR_QVALUE, min_fdr);
    carp(CARP_DETAILED_DEBUG, 
         "match %i cur fdr %f min fdr %f is decoy %i",
         match_idx, cur_fdr, min_fdr, cur_match->getNullPeptide() );
  }

  scored_type_[DECOY_XCORR_QVALUE] = true;
  return true;
}


/**
 * match_collection get, set method
 */

/**
 *\returns true, if the match collection has been scored by score_type
 */
bool MatchCollection::getScoredType(
  SCORER_TYPE_T score_type ///< the score_type (MATCH_SP, MATCH_XCORR) -in
  )
{
  return scored_type_[score_type];
}

/**
 * sets the score_type to value
 */
void MatchCollection::setScoredType(
  SCORER_TYPE_T score_type, ///< the score_type (MATCH_SP, MATCH_XCORR) -in
  bool value
  )
{
  scored_type_[score_type] = value;
}

/**
 *
 */
void MatchCollection::getCustomScoreNames(
  vector<string>& custom_score_names
  ) {
  custom_score_names.clear();

  if (match_total_ > 0) {

    match_[0]->getCustomScoreNames(custom_score_names);

  }

}

/**                                                                                                    
 * Set the filepath for all matches in the collection                                                  
 *\returns the associated file idx                                                                    
 */
int MatchCollection::setFilePath(
  const string& file_path  ///< File path to set                                                  
  ) {

  if (match_total_ > 0) {
    int file_idx = match_[0]->setFilePath(file_path);
    for (int match_idx = 1;match_idx < match_total_;match_idx++) {
      match_[match_idx]->setFileIndex(file_idx);
    }
    return file_idx;
  } else {
    carp(CARP_WARNING, "MatchCollection::setFilePath(): No matches in %s",file_path.c_str());
    return -1;
  }
}

/**
 *\returns true, if there is a  match_iterators instantiated by match collection 
 */
bool MatchCollection::getIteratorLock()
{
  return iterator_lock_;
}

/**
 *\returns the total match objects avaliable in current match_collection
 */
int MatchCollection::getMatchTotal()
{
  return match_total_;
}

bool MatchCollection::getHasDistinctMatches() {
  return has_distinct_matches_;
}

void MatchCollection::setHasDistinctMatches(bool distinct) {
  has_distinct_matches_ = distinct;
}


void MatchCollection::setExperimentSize(int size)
{
  experiment_size_ = size;
}

/**
 * \returns The total number of peptides searched for this spectrum,
 * target peptides for a target collection or decoy peptides for a
 * decoy collection.
 */
int MatchCollection::getExperimentSize()
{
  return experiment_size_;
}

/**
 * Sets the total number of target peptides searched for this
 * spectrum.  Only needs to be used by decoy match collections.
 */
void MatchCollection::setTargetExperimentSize(int numMatches){
  target_experiment_size_ = numMatches;
}

/**
 * \returns The number of target matches that this spectrum had.
 * Different than getExperimentSize() for decoy match collections.
 */
int MatchCollection::getTargetExperimentSize(){
  if( null_peptide_collection_ ){
    return target_experiment_size_;
  }
  return experiment_size_;
}

/**
 *\returns the top peptide count used in the logp_exp_sp in match_collection
 */
int MatchCollection::getTopFitSp()
{
  return top_fit_sp_;
}

/**
 *\returns the charge of the spectrum that the match collection was created
 */
int MatchCollection::getCharge()
{
  return zstate_.getCharge();
}

/**
 * Must have been scored by Xcorr, returns error if not scored by Xcorr
 *\returns the delta cn value(difference in top and second ranked Xcorr values)
 */
FLOAT_T MatchCollection::getDeltaCn()
{
  // Check if xcorr value has been scored, thus delta cn value is valid
  if(scored_type_[XCORR]){
    return delta_cn_;
  }
  else{
    carp(CARP_ERROR, "must score match_collection with XCORR to get delta cn value");
    return 0.0;
  }
}

/**
 * \brief Transfer the Weibull distribution parameters, including the
 * correlation from one match_collection to another.  No check to see
 * that the parameters have been estimated.
 */
void MatchCollection::transferWeibull(
  MatchCollection* from_collection,
  MatchCollection* to_collection
  ){
  to_collection->eta_ = from_collection->eta_;
  to_collection->beta_ = from_collection->beta_;
  to_collection->shift_ = from_collection->shift_;
  to_collection->correlation_ = from_collection->correlation_;
}

/**
 * \brief Prints out the pepxml header to the output stream
 * passed in as a parameter.
 */

void MatchCollection::printXmlHeader(
  FILE* output
  ){
  if (output == NULL ){
    return;
  }
  time_t hold_time;
  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  char* enz_str = enzyme_type_to_string(enzyme);
  char* database = get_string_parameter("protein-database");
  char* msms_file = get_string_parameter("ms2 file");
  char* absolute_msms_path;
  if (msms_file == NULL){
    absolute_msms_path = (char*) malloc(sizeof(char)*3);
    strcpy(absolute_msms_path, "NA");
  } else {
    #if DARWIN
    char path_buffer[PATH_MAX];
    absolute_msms_path =  realpath(msms_file, path_buffer);
    #else
    absolute_msms_path =  realpath(msms_file, NULL);
    #endif
    free(msms_file);
  }
  // Removes the extension from ms2 file path
  char* extension = strstr(absolute_msms_path, ".ms2");
  if (extension != NULL) (*extension) = '\0';
  
  
  MASS_TYPE_T isotopic_mass_type = get_mass_type_parameter("isotopic-mass");
  MASS_TYPE_T fragment_mass_type = get_mass_type_parameter("fragment-mass");

  const char* isotopic_mass;
  const char* fragment_mass;
  DIGEST_T digest = get_digest_type_parameter("digestion");
  int max_num_internal_cleavages;
  int min_number_termini;
  int missed_cleavage = get_int_parameter("missed-cleavages");
  if (missed_cleavage){
    max_num_internal_cleavages = get_int_parameter("max-length");
  } else {
    max_num_internal_cleavages = 0;
  }

  if (digest == FULL_DIGEST){
    min_number_termini = 2;
  } else if (digest == PARTIAL_DIGEST){
    min_number_termini = 1;
  } else {
    min_number_termini = 0;
  }
  
  if (isotopic_mass_type == AVERAGE){
    isotopic_mass = "average";
  } else {
    isotopic_mass = "monoisotopic";
  }

  if (fragment_mass_type == AVERAGE){
    fragment_mass =  "average";
  } else {
    fragment_mass =  "monoisotopic";
  }

  char* absolute_database_path = NULL;
  if( database != NULL ){
    bool use_index = is_directory(database);
    if( use_index == true ){
      char* fasta_name  = Index::getBinaryFastaName(database);
      free(database);
      database = fasta_name;
    }
#if DARWIN
    char path_buffer[PATH_MAX];
    absolute_database_path =  realpath(database, path_buffer);
#else
    absolute_database_path =  realpath(database, NULL);
#endif
    free(database);
  } 
  
  hold_time = time(0);

  const char* xsi = "http://www.w3.org/2001/XMLSchema-instance";
  const char* xmlns = "http://regis-web.systemsbiology.net/pepXML";
  const char* schema_location = "/usr/local/tpp/schema/pepXML_v110.xsd";
  fprintf(output, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(output, "<?xml-stylesheet type=\"text/xsl\" href=\"\"?>\n");
  fprintf(output, "<msms_pipeline_analysis date=\"%s\" xmlns=\"%s\""
          " xmlns:xsi=\"%s\" xsi:schemaLocation=\"%s %s\""
          " summary_xml=\"\">\n",
          ctime(&hold_time),
          xmlns, xsi, xmlns, schema_location);

  fprintf(output, "<msms_run_summary base_name=\"%s\" msManufacturer=\"%s\" "
          "msModel=\"%s\" msIonization=\"%s\" msAnalyzer=\"%s\" "
          "msDectector=\"%s\" raw_data_type=\"%s\" raw_data=\"%s\" >\n",
          absolute_msms_path,
          "NA", // TODO, dummy value
          "NA", // TODO, dummy value
          "NA", // TODO, dummy value
          "NA", // TODO, dummy value
          "NA", // TODO, dummy value
          "NA", // TODO, dummy value
          "NA" // TODO, dummy value
          );
  

  fprintf(output, "<sample_enzyme name=\"%s\">\n</sample_enzyme>\n", enz_str);

  fprintf(output, "<search_summary base_name=\"%s\" search_engine=\"%s\" "
          "precursor_mass_type=\"%s\" fragment_mass_type=\"%s\" "
          "out_data_type=\"%s\" out_data=\"%s\" search_id=\"%i\" >\n",
          absolute_msms_path,
          "Crux",
          isotopic_mass, // isotopic mass type is precursor mass type?
          fragment_mass,
          "NA", // TODO, dummy value
          "NA",
          1 // TODO, dummy value
          );

  
  fprintf(output, "<search_database local_path=\"%s\" type=\"%s\" />\n", 
          absolute_database_path, 
          "AA"
          );
  fprintf(output, "<enzymatic_search_constraint enzyme=\"%s\" "
          "max_num_internal_cleavages=\"%i\" min_number_termini=\"%i\"/>\n",
          enz_str,
          max_num_internal_cleavages,
          min_number_termini
          );

#ifndef DARWIN
  free(absolute_msms_path);
  free(absolute_database_path);
#endif
  free(enz_str);


  char aa_str[2];
  aa_str[1] = '\0';
  int alphabet_size = (int)'A'+ ((int)'Z'-(int)'A');
  MASS_TYPE_T isotopic_type = get_mass_type_parameter("isotopic-mass");
  int aa = 0;

  // static amino acid modifications
  for (aa = (int)'A'; aa < alphabet_size-1; aa++){
    aa_str[0] = (char)aa;
    double mod = get_double_parameter(aa_str);
    double mass = get_mass_amino_acid(aa, isotopic_type);
    
    if (mod != 0 ){
      fprintf(output, "<aminoacid_modification aminoacid=\"%s\" mass=\"%f\" "
              "massdiff=\"%f\" variable=\"%s\" />\n",
              aa_str,
              mass + mod,
              mod,
              "N" // N if static modification
              );      
    }
  }
  
  // variable amino acid modifications
  AA_MOD_T** mod_list = NULL;
  int num_mods = get_aa_mod_list(&mod_list);
  for (int mod_idx = 0; mod_idx < num_mods; mod_idx++){
    FLOAT_T mod_mass = aa_mod_get_mass_change(mod_list[mod_idx]);
    
    bool* aas_modified = aa_mod_get_aa_list(mod_list[mod_idx]);
    for (int aa_idx = 0; aa_idx < AA_LIST_LENGTH; aa_idx++){
      if (aas_modified[aa_idx] == true ){
        int aa = (aa_idx+'A');
        FLOAT_T aa_mass = get_mass_amino_acid(aa , isotopic_type);
        fprintf(output, "<aminoacid_modification aminoacid=\"%c\" mass=\"%f\" "
                "massdiff=\"%f\" variable=\"%s\" />\n",
                aa,
                aa_mass + mod_mass,
                mod_mass,
                "Y" // Y if variable modification
                );    
      }
    }

  }

  // terminal modifciations
  // variable
  num_mods = get_c_mod_list(&mod_list); // variable c mods
  for(int mod_idx = 0; mod_idx < num_mods; mod_idx++){
    FLOAT_T mod_mass = aa_mod_get_mass_change(mod_list[mod_idx]);
    fprintf(output, "<terminal_modification terminus=\"c\" "
            "mass=\"%f\" massdiff=\"%f\" variable=\"Y\" />\n",
            MASS_OH + mod_mass,
            mod_mass
            );
  }
  num_mods = get_n_mod_list(&mod_list); // variable n mods
  for(int mod_idx = 0; mod_idx < num_mods; mod_idx++){
    FLOAT_T mod_mass = aa_mod_get_mass_change(mod_list[mod_idx]);
    fprintf(output, "<terminal_modification terminus=\"n\" "
            "mass=\"%f\" massdiff=\"%f\" variable=\"Y\" />\n",
            MASS_H_MONO + mod_mass,
            mod_mass
            );
  }
  // fixed
  if( get_num_fixed_mods() != 0 ){
    get_all_aa_mod_list(&mod_list);
    int fixed_mod_idx = get_fixed_mod_index(N_TERM); // fixed n mods
    if( fixed_mod_idx > -1 ){
      fprintf(output, "<terminal_modification terminus=\"n\" "
              "mass=\"?\" massdiff=\"%f\" variable=\"N\" />\n",
              aa_mod_get_mass_change(mod_list[fixed_mod_idx])
              );
    }

    fixed_mod_idx = get_fixed_mod_index(C_TERM); // fixed c mods
    if( fixed_mod_idx > -1 ){
      fprintf(output, "<terminal_modification terminus=\"c\" "
              "mass=\"?\" massdiff=\"%f\" variable=\"N\" />\n",
              aa_mod_get_mass_change(mod_list[fixed_mod_idx])
              );
    }

  }


  print_parameters_xml(output);
  
  fprintf(output, "</search_summary>\n");

}

/**
 * Write header for .sqt file.  Assumes only sequest-search is writing
 * this file type.
 */
void MatchCollection::printSqtHeader(
 FILE* output, 
 const char* type, 
 string database,
 int num_proteins,
 bool exact_pval_search){  
  if( output == NULL ){
    return;
  }

  time_t hold_time;
  hold_time = time(0);

  bool decoy = false;
  if( strcmp(type, "decoy") == 0 ){
    decoy = true;
  }

  fprintf(output, "H\tSQTGenerator Crux\n");
  fprintf(output, "H\tSQTGeneratorVersion 1.0\n");
  fprintf(output, "H\tComment Crux was written by...\n");
  fprintf(output, "H\tComment ref...\n");
  fprintf(output, "H\tStartTime\t%s", ctime(&hold_time));
  fprintf(output, "H\tEndTime                               \n");

  fprintf(output, "H\tDatabase\t%s\n", database.c_str());

  if(decoy){
  fprintf(output, "H\tComment\tDatabase shuffled; these are decoy matches\n");
  }
  fprintf(output, "H\tDBSeqLength\t?\n");
  fprintf(output, "H\tDBLocusCount\t%d\n", num_proteins);

  MASS_TYPE_T mass_type = get_mass_type_parameter("isotopic-mass");
  char temp_str[64];
  mass_type_to_string(mass_type, temp_str);
  fprintf(output, "H\tPrecursorMasses\t%s\n", temp_str);
  
  mass_type = get_mass_type_parameter("fragment-mass");
  mass_type_to_string(mass_type, temp_str);
  fprintf(output, "H\tFragmentMasses\t%s\n", temp_str); //?????????

  double tol = get_double_parameter("precursor-window");
  fprintf(output, "H\tAlg-PreMasTol\t%.1f\n",tol);
  fprintf(output, "H\tAlg-FragMassTol\t%.2f\n", 
          get_double_parameter("mz-bin-width") / 2.0);
  fprintf(output, "H\tAlg-XCorrMode\t0\n");

  SCORER_TYPE_T score = get_scorer_type_parameter("prelim-score-type");
  fprintf(output, "H\tComment\tpreliminary algorithm %s\n", 
          scorer_type_to_string(score));

  score = get_scorer_type_parameter("score-type");
  fprintf(output, "H\tComment\tfinal algorithm %s\n", 
          scorer_type_to_string(score));

  int aa = 0;
  char aa_str[2];
  aa_str[1] = '\0';
  int alphabet_size = (int)'A' + ((int)'Z'-(int)'A');
  MASS_TYPE_T isotopic_type = get_mass_type_parameter("isotopic-mass");

  for(aa = (int)'A'; aa < alphabet_size -1; aa++){
    aa_str[0] = (char)aa;
    double mod = get_double_parameter(aa_str);
    if( mod != 0 ){
      //      double mass = mod + get_mass_amino_acid(aa, isotopic_type);
      double mass = get_mass_amino_acid(aa, isotopic_type);
      fprintf(output, "H\tStaticMod\t%s=%.3f\n", aa_str, mass);
    }
  }
  // print dynamic mods, if any
  // format DiffMod <AAs><symbol>=<mass change>
  AA_MOD_T** aa_mod_list = NULL;
  int num_mods = get_all_aa_mod_list(&aa_mod_list);
  int mod_idx = 0;
  for(mod_idx = 0; mod_idx < num_mods; mod_idx++){
    
    AA_MOD_T* aamod = aa_mod_list[mod_idx];
    char* aa_list_str = aa_mod_get_aa_list_string(aamod);
    char aa_symbol = aa_mod_get_symbol(aamod);
    double mass_dif = aa_mod_get_mass_change(aamod);

    fprintf(output, "H\tDiffMod\t%s%c=%+.2f\n", aa_list_str, 
            aa_symbol, mass_dif);
    free(aa_list_str);
  }
  num_mods = get_c_mod_list(&aa_mod_list);
  for(mod_idx = 0; mod_idx < num_mods; mod_idx++){
    AA_MOD_T* aamod = aa_mod_list[mod_idx];
    char aa_symbol = aa_mod_get_symbol(aamod);

    fprintf(output, "H\tComment\tMod %c is a C-terminal modification\n",
            aa_symbol);
  }

  num_mods = get_n_mod_list(&aa_mod_list);
  for(mod_idx = 0; mod_idx < num_mods; mod_idx++){
    AA_MOD_T* aamod = aa_mod_list[mod_idx];
    char aa_symbol = aa_mod_get_symbol(aamod);

    fprintf(output, "H\tComment\tMod %c is a N-terminal modification\n",
            aa_symbol);
  }



  //for letters in alphabet
  //  double mod = get_double_parameter(letter);
  //  if mod != 0
  //     double mass = mod + getmass(letter);
  //     fprintf(output, "H\tStaticMod\t%s=%.3f\n", letter, mass);
  //  fprintf(output, "H\tStaticMod\tC=160.139\n");
  fprintf(output, "H\tAlg-DisplayTop\t%d\n", 
          //          get_int_parameter("max-sqt-result")); 
          get_int_parameter("top-match")); 
  // this is not correct for an sqt from analzyed matches

  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  DIGEST_T digestion = get_digest_type_parameter("digestion");
  char* enz_str = enzyme_type_to_string(enzyme);
  char* dig_str = digest_type_to_string(digestion);
  char custom_str[SMALL_BUFFER];
  if( enzyme == CUSTOM_ENZYME){
    char* rule = get_string_parameter("custom-enzyme");
    sprintf(custom_str, ", custom pattern: %s", rule);
  }else{
    custom_str[0] = 0;
  }
  fprintf(output, "H\tEnzymeSpec\t%s-%s%s\n", enz_str, dig_str, custom_str);
  free(enz_str);
  free(dig_str);
  // write a comment that says what the scores are
  fprintf(output, "H\tLine fields: S, scan number, scan number, "
          "charge, 0, server, experimental mass, total ion intensity, "
          "lowest Sp, number of matches\n");
  if (exact_pval_search) {
    fprintf(output, "H\tLine fields: M, rank by xcorr score, rank by sp score, "
            "peptide mass, deltaCn, exact P-value, recalibrated xcorr, sp score, number ions matched, "
            "total ions compared, sequence, validation status\n");
  } else {
    fprintf(output, "H\tLine fields: M, rank by xcorr score, rank by sp score, "
            "peptide mass, deltaCn, xcorr score, sp score, number ions matched, "
            "total ions compared, sequence, validation status\n");
  }
}

/**
 * Print the header line for a tab-delimited file.
 */
void MatchCollection::printTabHeader(FILE* output){

  if( output == NULL ){
    return;
  }

  int column_idx;
  for (column_idx = 0; column_idx < NUMBER_MATCH_COLUMNS; column_idx++) {
    fprintf(output, "%s", get_column_header(column_idx));
    if (column_idx < NUMBER_MATCH_COLUMNS - 1) {
      fprintf(output, "\t");
    } else {
      fprintf(output, "\n");
    }
  }
}


/**
 * Print Footer lines for xml files
 */
void MatchCollection::printXmlFooter(FILE* output){
  if (output == NULL ){
    return;
  }
  
  fprintf(output, "</msms_run_summary>\n");
  fprintf(output, "</msms_pipeline_analysis>\n");
}


/**
 * \brief Print the given match collection for several spectra to
 * xml files only. Takes the spectrum information from the
 * matches in the collection. At least for now, prints all matches in
 * the collection rather than limiting by top-match parameter. 
 */
void MatchCollection::printMultiSpectraXml(
  PepXMLWriter* output
){
  carp(CARP_DETAILED_DEBUG, "Writing matches to xml file");
  // reuse these for each match
  vector<string> protein_ids;
  vector<string> protein_descriptions;
  double* scores = new double[NUMBER_SCORER_TYPES];

  int match_idx = 0;
  int num_matches = match_total_;
  for (match_idx = 0; match_idx < num_matches; match_idx++){
    Match* cur_match = match_[match_idx];
    bool is_decoy = cur_match->getNullPeptide();
    Spectrum* spectrum = cur_match->getSpectrum();
    double cur_ln_experiment_size=0;
    if (! is_decoy){
      int* ranks =new int[NUMBER_SCORER_TYPES];
      ranks[XCORR]=-1; 
      if( scored_type_[XCORR] ){
        ranks[XCORR] = cur_match->getRank(XCORR);
      }else if(scored_type_[SP]){
        ranks[SP]=cur_match->getRank(SP);
        scores[SP]=cur_match->getScore(SP);
      }
      char* peptide_sequence = cur_match->getSequence();
      char* mod_peptide_sequence = cur_match->getModSequenceStrWithMasses(
                             get_mass_format_type_parameter("mod-mass-format"));
      Peptide* peptide = cur_match->getPeptide();
      char* flanking_aas = peptide->getFlankingAAs();
      int num_proteins = peptide->getProteinInfo(protein_ids, 
                                                 protein_descriptions);
      for(int score_idx = 0; score_idx < NUMBER_SCORER_TYPES; score_idx++){
        if( scored_type_[score_idx] == true ){
          scores[score_idx] = cur_match->getScore((SCORER_TYPE_T)score_idx);
          ranks[score_idx]=cur_match->getRank((SCORER_TYPE_T)score_idx);
         
        }
  
      }
      unsigned num_matches = getTargetExperimentSize(); 
      if(isDecoy())
        num_matches=getExperimentSize();
      else 
        num_matches=getTargetExperimentSize(); 
      output->writePSM(spectrum->getFirstScan(),
        spectrum->getFilename(),
        cur_match->getNeutralMass(),
        cur_match->getCharge(),
        ranks,
        peptide_sequence,
        mod_peptide_sequence,
        peptide->getPeptideMass(),
        num_proteins,
        flanking_aas,
        protein_ids,
        protein_descriptions,
        cur_match->getDeltaCn(),
        scored_type_,
        scores,
        cur_match->getBYIonMatched(),
        cur_match->getBYIonPossible(),
        num_matches
      );
    }
  }
  
}

/**
 * \brief Print the psm features to file in xml format
 *
 * Prints a spectrum_query tag which encompasses the search_hit tag
 * which represents peptide to spectra match.
 *
 * returns true, if succesfully printed xml format of PSMs, else false
 *
 */

bool MatchCollection::printXml(
  PepXMLWriter* output,
  int top_match,
  Spectrum* spectrum,
  SCORER_TYPE_T main_score
  )
{
  if ( output == NULL || spectrum == NULL || match_total_ == 0){
    return false;
  }

  // calculate delta_cn and populate fields in the matches
  calculateDeltaCn();

  // for deciding when to quit
  int count = 0;
  int last_rank = 0;

  // reuse these for each match
  vector<string> protein_ids;
  vector<string> protein_descriptions;
  bool* scores_computed = new bool[NUMBER_SCORER_TYPES];
  
  for(int score_idx = 0; score_idx < NUMBER_SCORER_TYPES; score_idx++){
    scores_computed[score_idx] = false;
  }
  scores_computed[main_score] = true;
  if( scored_type_[SP]) {
    scores_computed[SP] = true;
  }
  scores_computed[TIDE_SEARCH_EXACT_PVAL] = exact_pval_search_;
  scores_computed[TIDE_SEARCH_REFACTORED_XCORR] = exact_pval_search_;
  scored_type_[TIDE_SEARCH_EXACT_PVAL] = exact_pval_search_;
  scored_type_[TIDE_SEARCH_REFACTORED_XCORR] = exact_pval_search_;
  scores_computed[main_score] = !exact_pval_search_;

  double* scores = new double[NUMBER_SCORER_TYPES];
  int* ranks=new int[NUMBER_SCORER_TYPES];

  Match* match = NULL;
  // create match iterator
  // true: return match in sorted order of main_score type
  MatchIterator* match_iterator = 
    new MatchIterator(this, main_score, true);
  // iterate over matches
  while(match_iterator->hasNext()){
    match = match_iterator->next();
    int cur_rank = match->getRank(main_score);   
    if(scored_type_[XCORR])
      ranks[XCORR]=match->getRank(XCORR);
    if(scored_type_[SP]){
      ranks[SP]= match->getRank(SP);
      scores[SP]= match->getScore(SP);
    }
    // print if we haven't reached the limit
    // or if we are at the limit but this match is a tie with the last
    if( count < top_match || last_rank == cur_rank ){
      
      char* peptide_sequence = match->getSequence();
      char* mod_peptide_sequence = match->getModSequenceStrWithMasses(
                            get_mass_format_type_parameter("mod-mass-format"));
      Peptide* peptide = match->getPeptide();
      char* flanking_aas = peptide->getFlankingAAs();
      int num_proteins = peptide->getProteinInfo(protein_ids, 
                                                 protein_descriptions);
     for(int score_idx=0; score_idx < NUMBER_SCORER_TYPES; score_idx++){
      if(scored_type_[score_idx])
        scores[score_idx] = match->getScore((SCORER_TYPE_T)score_idx);
      
     }   
     unsigned num_matches= getTargetExperimentSize(); 
     if(isDecoy())
       num_matches= getExperimentSize(); 
     else
       num_matches= getTargetExperimentSize(); 
      output->writePSM(spectrum->getFirstScan(),
        spectrum->getFilename(),
        zstate_.getNeutralMass(),
        zstate_.getCharge(),
        ranks,
        peptide_sequence,
        mod_peptide_sequence,
        peptide->getPeptideMass(),
        num_proteins,
        flanking_aas,
        protein_ids,
        protein_descriptions,
        match->getDeltaCn(),
        scores_computed,
        scores,
        match->getBYIonMatched(),
        match->getBYIonPossible(), 
        num_matches
      );
      count++;
      last_rank = cur_rank;
      free(peptide_sequence);
      free(mod_peptide_sequence);
      free(flanking_aas);
    } else if( count >= top_match && last_rank != cur_rank ) {
      break;
    } // else see if there is one more tie to print

  }// next match
  
  carp(CARP_DETAILED_DEBUG, "printed %d out of %d xml matches", 
       count, match_total_);

  delete match_iterator;
  delete scores_computed;
  delete scores;
  delete ranks;

  return true;
}


/**
 * \brief Print the psm features to file in sqt format.
 *
 * Prints one S line, 'top_match' M lines, and one locus line for each
 * peptide source of each M line.
 * Assumes one spectrum per match collection.  Only crux
 * sequset-search produces sqt files so the two scores are always Sp
 * and xcorr.
 * Possible side effects: Calculates delta_cn and stores in each match
 *\returns true, if sucessfully print sqt format of the PSMs, else false 
 */
bool MatchCollection::printSqt(
  FILE* output,                  ///< the output file -out
  int top_match,                 ///< the top matches to output -in
  Spectrum* spectrum           ///< the spectrum to print sqt -in
  )
{

  if( output == NULL || spectrum == NULL || match_total_ == 0 ){
    return false;
  }

  SpectrumZState& zstate = zstate_; 
  int num_matches = experiment_size_;

  // calculate delta_cn and populate fields in the matches
  calculateDeltaCn();

  // First, print spectrum info
  spectrum->printSqt(output, num_matches, zstate);
  
  Match* match = NULL;
  
  // create match iterator
  // true: return match in sorted order of main_score type
  MatchIterator* match_iterator =
    new MatchIterator(this, XCORR, true);
  
  // Second, iterate over matches, prints M and L lines
  while(match_iterator->hasNext()){
    match = match_iterator->next();    

    // print only up to max_rank_result of the matches
    if( match->getRank(XCORR) > top_match ){
      break;
    }// else

    match->printSqt(output);

  }// next match
  
  // make sure top_scoring_sp_ has been set
  if( top_scoring_sp_ == NULL){
    carp(CARP_DEBUG, "Top scoring SP was not set.");
  } else if( top_scoring_sp_->getRank(XCORR) > top_match ){
    // print the match with Sp rank==1 if its xcorr rank > top_match rank.  
    top_scoring_sp_->printSqt(output);
  }
  
  delete match_iterator;
  
  return true;
}

/**
 * \brief Print the psm features to file in tab delimited format.
 *
 * Matches will be sorted by main_score and the ranks of those scores
 * will be used to determine how many matches are printed for each
 * spectrum.
 * \returns true, if sucessfully print tab-delimited format of the
 * PSMs, else false  
 */
bool MatchCollection::printTabDelimited(
  MatchFileWriter* output,                  ///< the output file -out
  int top_match,                 ///< the top matches to output -in
  Spectrum* spectrum,          ///< the spectrum to print sqt -in
  SCORER_TYPE_T main_score       ///< the main score to report -in
  )
{

  if( output == NULL || spectrum == NULL || match_total_ == 0 ){
    return false;
  }
  int num_target_matches = getTargetExperimentSize();
  int num_decoy_matches = getExperimentSize();
  int scan_num = spectrum->getFirstScan();
  FLOAT_T spectrum_precursor_mz = spectrum->getPrecursorMz();

  Match* match = NULL;
  
  // create match iterator
  // true: return match in sorted order of main_score type
  MatchIterator* match_iterator = 
    new MatchIterator(this, main_score, true);
  int count = 0;
  int last_rank = 0;

  // iterate over matches
  while(match_iterator->hasNext()){
    match = match_iterator->next();
    int cur_rank = match->getRank(main_score);

    // print if we haven't reached the limit
    // or if we are at the limit but this match is a tie with the last
    if( count < top_match || last_rank == cur_rank ){

      match->printTab(this, output, scan_num, 
                      spectrum_precursor_mz, 
                      num_target_matches, num_decoy_matches);
      count++;
      last_rank = cur_rank;
    } else if( count >= top_match && last_rank != cur_rank ) {
      break;
    } // else see if there is one more tie to print

  }// next match
  
  carp(CARP_DETAILED_DEBUG, "printed %d out of %d tab matches", 
       count, num_target_matches);

  delete match_iterator;
  
  return true;
}

/**
 * Retrieve the calibration parameter eta.
 */
FLOAT_T MatchCollection::getCalibrationEta()
{
  return eta_;
}

/**
 * Retrieve the calibration parameter beta.
 */
FLOAT_T MatchCollection::getCalibrationBeta()
{
  return beta_;
}

/**
 * Retrieve the calibration parameter shift.
 */
FLOAT_T MatchCollection::getCalibrationShift()
{
  return shift_;
}

/**
 * Retrieve the calibration correlation.
 */
FLOAT_T MatchCollection::getCalibrationCorr()
{
  return correlation_;
}

/**
 * \brief Print the given match collection for several spectra to
 * tab-delimited files only.  Takes the spectrum information from the
 * matches in the collection.  At least for now, prints all matches in
 * the collection rather than limiting by top-match parameter.  Uses
 * SP as preliminary score and XCORR as main score.
 */
void MatchCollection::printMultiSpectra(
 MatchFileWriter* tab_file, 
 MatchFileWriter* decoy_tab_file
  ){

  carp(CARP_DETAILED_DEBUG, "Writing matches to file");

  // if file location is target (i.e. tdc=T), print all to target
  MatchFileWriter* decoy_file = decoy_tab_file;
  if( get_boolean_parameter("tdc") == true ){
    decoy_file = tab_file;
  }

  // for each match, get spectrum info, determine if decoy, print
  int match_idx = 0;
  int num_matches = match_total_;
  for(match_idx = 0; match_idx < num_matches; match_idx++){
    Match* cur_match = match_[match_idx];
    bool is_decoy = cur_match->getNullPeptide();
    Spectrum* spectrum = cur_match->getSpectrum();
    int scan_num = spectrum->getFirstScan();
    FLOAT_T mz = spectrum->getPrecursorMz();
    FLOAT_T num_psm_per_spec = cur_match->getLnExperimentSize();
    num_psm_per_spec = expf(num_psm_per_spec) + 0.5; // round to nearest int
    int num_target_psm_per_spec = cur_match->getTargetExperimentSize();

    if( is_decoy ){
      cur_match->printTab(this, decoy_file, scan_num, mz, 
                          (int)num_psm_per_spec, num_target_psm_per_spec);
    }
    else{
      cur_match->printTab(this, tab_file, scan_num, mz,
                          (int)num_psm_per_spec, num_target_psm_per_spec);
    }

  }

}

/*******************************************
 * match_collection post_process extension
 ******************************************/

/**
 * Create a list of file names that we will look for in the psm directory.
 */
void set_possible_names(vector<string>& possible_names, SET_TYPE_T type){

  // if a specific file has been requested, return just that file name
  const char* psm_filename = get_string_parameter_pointer("input PSMs");
  if( psm_filename == NULL || strcmp(psm_filename, "__NULL_STR") != 0 ){
    possible_names.push_back(psm_filename);
    return;
  }

  // else, decide on the search string for targets or decoys
  switch(type){
  case SET_TARGET:
    possible_names.push_back("search.target.txt");
    possible_names.push_back("sequest.target.txt");
    break;
  case SET_DECOY1:
    possible_names.push_back("search.decoy.txt");
    possible_names.push_back("search.decoy-1.txt");
    possible_names.push_back("sequest.decoy.txt");
    possible_names.push_back("sequest.decoy-1.txt");
    break;
  case SET_DECOY2:
    possible_names.push_back("search.decoy-2.txt");
    possible_names.push_back("sequest.decoy-2.txt");
    break;
  case SET_DECOY3:
    possible_names.push_back("search.decoy-3.txt");
    possible_names.push_back("sequest.decoy-3.txt");
    break;
  }
}

/**
 * Read files in the directory and return the names of target or
 * decoy files to use for post-search commands.
 * \returns Vector parameter filled with names of target or decoy
 * files.
 */
void get_target_decoy_filenames(vector<string>& target_decoy_names,
                                DIR* directory,
                                SET_TYPE_T type){
  if( directory == NULL ){
    carp(CARP_FATAL, "Cannot read files from NULL directory.");
  }

  // first see if there is a specific file to open
  const char* psm_file = get_string_parameter_pointer("input PSMs");
  if( psm_file != NULL && strcmp(psm_file, "__NULL_STR") != 0 ){
    // strip off the path
    char** name_path = parse_filename_path(psm_file);
    target_decoy_names.push_back(name_path[0]);
    free(name_path[0]);
    free(name_path[1]);
    free(name_path);
    return;
  }

  // look for both files from search-for-matches and sequest-search
  vector<string> possible_names;

  set_possible_names(possible_names, type);

  // open the directory
  struct dirent* directory_entry = NULL;
  rewinddir(directory);
  // for each file, compare to each name, if it matches, add to the list
  while((directory_entry = readdir(directory))){
    for(int name_idx = 0; name_idx < (int)possible_names.size(); name_idx++){
      string filename = directory_entry->d_name;
      if( filename.find(possible_names[name_idx]) != string::npos ){
        target_decoy_names.push_back(filename);
      } 
    }
  }

  // check that files are only sequest or search, not both
  bool found_search = false;
  bool found_sequest = false;
  if( target_decoy_names.size() > 1 ){
    for(int name_idx = 0; name_idx < (int)target_decoy_names.size();name_idx++){
      // don't look for just "search" and "sequest" in case they are in
      // the fileroot
      if( target_decoy_names[name_idx].find(possible_names.front()) 
          != string::npos ){
        found_search = true;
      }
      if( target_decoy_names[name_idx].find(possible_names.back()) 
          != string::npos ){
        found_sequest = true;
      }
    }
  }

  if( found_search && found_sequest ){
    carp(CARP_FATAL, "Cannot analyze results from both crux search-for-matches "
         " and sequest-search.  Please remove one from the directory.");
  }
  // check that headers are all the same??
}


/**
 * parse all the match objects and add to match collection
 *\returns true, if successfully parse all PSMs in result_file, else false
 */
bool MatchCollection::extendTabDelimited(
  Database* database, ///< the database holding the peptides -in
  MatchFileReader& result_file,   ///< the result file to parse PSMs -in
  Database* decoy_database ///< the database holding the decoy peptides -in
  )
{
  Match* match = NULL;

  FLOAT_T delta_cn = 0;
  FLOAT_T ln_delta_cn = 0;
  FLOAT_T ln_experiment_size = 0;

  // only for post_process_collections
  if(!post_process_collection_){
    carp(CARP_ERROR, "Must be a post process match collection to extend.");
    return false;
  }

  while (result_file.hasNext()) {

    /*** get spectrum specific features ***/
    zstate_.setNeutralMass(
      result_file.getFloat(SPECTRUM_NEUTRAL_MASS_COL),
      result_file.getInteger(CHARGE_COL));
    scored_type_[DELTA_CN] = scored_type_[DELTA_CN] || !result_file.empty(DELTA_CN_COL);
    delta_cn = result_file.getFloat(DELTA_CN_COL);
    if (delta_cn <= 0.0) {
      ln_delta_cn = 0;
    } else {
      ln_delta_cn = logf(delta_cn);
    }
    if (!result_file.empty(DISTINCT_MATCHES_SPECTRUM_COL)) {
      has_distinct_matches_ = true;
      ln_experiment_size = log(result_file.getFloat(DISTINCT_MATCHES_SPECTRUM_COL));
    } else if (!result_file.empty(MATCHES_SPECTRUM_COL)) {
      ln_experiment_size = log(result_file.getFloat(MATCHES_SPECTRUM_COL));
    } else {
      ln_experiment_size = 0;
    }

    //TODO: Parse all boolean indicators for scores
    scored_type_[SP] = !result_file.empty(SP_SCORE_COL);

    scored_type_[XCORR] = !result_file.empty(XCORR_SCORE_COL);
    
    scored_type_[EVALUE] = !result_file.empty(EVALUE_COL);
    
    scored_type_[DECOY_XCORR_QVALUE] = !result_file.empty(DECOY_XCORR_QVALUE_COL);

/* TODO
    match_collection -> 
      scored_type[LOGP_WEIBULL_XCORR] = 
      result_file.getString("logp weibull xcorr") != "";
*/

    scored_type_[LOGP_BONF_WEIBULL_XCORR] = 
      !result_file.empty(PVALUE_COL);

    scored_type_[PERCOLATOR_QVALUE] = 
      !result_file.empty(PERCOLATOR_QVALUE_COL);

    scored_type_[PERCOLATOR_SCORE] = 
      !result_file.empty(PERCOLATOR_SCORE_COL);

    scored_type_[LOGP_QVALUE_WEIBULL_XCORR] = 
      !result_file.empty(WEIBULL_QVALUE_COL);
  
    scored_type_[QRANKER_SCORE] = 
      !result_file.empty(QRANKER_SCORE_COL);
    
    scored_type_[QRANKER_QVALUE] = 
      !result_file.empty(QRANKER_QVALUE_COL);

    scored_type_[BARISTA_SCORE] =
      !result_file.empty(BARISTA_SCORE_COL);

    scored_type_[BARISTA_QVALUE] =
      !result_file.empty(BARISTA_QVALUE_COL);

    post_scored_type_set_ = true;

    // parse match object
    match = Match::parseTabDelimited(result_file, database, decoy_database);
    if (match == NULL) {
      carp(CARP_ERROR, "Failed to parse tab-delimited PSM match");
      return false;
    }

    //set all spectrum specific features to parsed match
    match->setZState(zstate_);
    match->setDeltaCn(delta_cn);
    match->setDeltaLCn(ln_delta_cn);
    match->setLnExperimentSize(ln_experiment_size);    
    //add match to match collection.
    addMatchToPostMatchCollection(match);
    //increment pointer.
    result_file.next();
  }

  return true;
}



/**
 * \brief Adds the match to match_collection by copying the pointer.
 * 
 * No new match is allocated.  Match_collection total_matches must not
 * exceed the _MAX_NUMBER_PEPTIDES. 
 * \returns true if successfully adds the match to the
 * match_collection, else false 
 */
bool MatchCollection::addMatch(
  Match* match ///< the match to add -in
  )
{
  if( match == NULL ){
    carp(CARP_FATAL, "Cannot add NULL match.");
  }

  // check if enough space for peptide match
  if(match_total_ >= _MAX_NUMBER_PEPTIDES){
    carp(CARP_FATAL, "Cannot add to match collection; count exceeds limit: %d", 
         _MAX_NUMBER_PEPTIDES);
  }

  // add a new match to array
  match_[match_total_] = match;
  match->incrementPointerCount();
  
  // increment total rich match count
  ++match_total_;

  
  return true;
}

/**
 * Adds the match object to match_collection
 * Must not exceed the _MAX_NUMBER_PEPTIDES to be match added
 * Only for post_process (i.e. post search) match_collections.  Keeps
 * track of all peptides in a hash table.
 * \returns true if successfully adds the match to the
 * match_collection, else false 
 */
// this method renamed so that a more general add_match_to_match_collection could be implemented
bool MatchCollection::addMatchToPostMatchCollection(
  Match* match ///< the match to add -in
  )
{
  if( match == NULL ){
    carp(CARP_FATAL, "Cannot add NULL match to NULL collection.");
  }

  // only for post_process_collections
  if(!post_process_collection_){
    carp(CARP_ERROR, "Must be a post process match collection to add a match.");
    return false;
  }

  // check if enough space for peptide match
  if(match_total_ >= _MAX_NUMBER_PEPTIDES){
    carp(CARP_ERROR, "Rich match count exceeds max match limit: %d", 
         _MAX_NUMBER_PEPTIDES);
    return false;
  }

  // add a new match to array
  match_[match_total_] = match;
  match->incrementPointerCount();
  
  // increment total rich match count
  ++match_total_;
  
  // DEBUG, print total peptided scored so far
  if(match_total_ % 1000 == 0){
    carp(CARP_INFO, "parsed PSM: %d", match_total_);
  }
  
  return true;
}

/**
 * Fill the match objects score with the given the array, and populate
 * the corresponding ranks.  The match object order must not have been
 * altered since scoring.  The result array size must equal the number
 * of matches in the given match collection.  After the function
 * completes, the match collection is sorted by the specified score,
 * unless preserve_order is set to true.
 */
void MatchCollection::fillResult(
  double* results,           ///< array of scores -in
  SCORER_TYPE_T score_type,  ///< The score type of the results to fill -in
  bool preserve_order   ///< preserve match order?
  )
{
  Match** match_array = NULL;
  SCORER_TYPE_T score_type_old = last_sorted_;

  // iterate over match object in collection, set scores
  int match_idx = 0;
  for(; match_idx < match_total_; ++match_idx){
    Match* match = match_[match_idx];
    match->setScore(score_type, results[match_idx]);    
  }
  
  // if need to preserve order store a copy of array in original order 
  if(preserve_order){
    match_array = (Match**)mycalloc(match_total_, sizeof(Match*));
    for(match_idx=0; match_idx < match_total_; ++match_idx){
      match_array[match_idx] = match_[match_idx];
    }
  }

  // populate the rank of match_collection
  if(!populateMatchRank(score_type)){
    carp(CARP_FATAL, "failed to populate match rank in match_collection");
  }
  
  // restore match order.
  if(preserve_order){
    for(match_idx=0; match_idx < match_total_; ++match_idx){
      match_[match_idx] = match_array[match_idx];
    }
    last_sorted_ = score_type_old;
    free(match_array);
  }

  scored_type_[score_type] = true;
}

/**
 * Process run specific features from all the PSMs
 */
void MatchCollection::processRunSpecificFeatures() {
}


/**
 * \brief Calculate the delta_cn of each match and populate the field.
 * 
 * Delta_cn is the normalized difference between xcorrs of different ranks.
 * match[i] = (match[i] - match[i+1]) / match[i].
 * Sorts match_collection by xcorr, if necessary.
 * 
 */
bool MatchCollection::calculateDeltaCn(){

  if( scored_type_[XCORR] == false ){
    carp(CARP_WARNING, 
      "Delta_cn not calculated because match collection not scored for xcorr");
    return false;
  }

  // sort, if not already
  // N.B. Can't use sort_match_collection because iterator already exists!
  Match** matches = match_;
  int num_matches = match_total_;
  if( last_sorted_ != XCORR ){
    qsortMatch(matches, num_matches, (QSORT_COMPARE_METHOD)compareXcorr);
    last_sorted_ = XCORR;
  }

  FLOAT_T last_xcorr=0.0;
  FLOAT_T delta_cn = 0.0;
  FLOAT_T delta_lcn = 0.0;
  FLOAT_T next_xcorr=0.0;
  FLOAT_T current_xcorr = 0 ; 
  if(num_matches>1){
    last_xcorr = matches[num_matches-1]->getScore(XCORR);
    for (size_t idx = 0 ;idx < num_matches;idx++) { 
      current_xcorr = matches[idx]->getScore(XCORR);
      if (idx+1<=num_matches-1)
        next_xcorr=matches[idx+1]->getScore(XCORR);
      delta_cn = (current_xcorr - next_xcorr) / max(current_xcorr, (FLOAT_T)1.0);
      delta_lcn = (current_xcorr - last_xcorr) / max(current_xcorr, (FLOAT_T)1.0);
    
      if(fabs(delta_cn)== numeric_limits<FLOAT_T>::infinity()){
        carp(CARP_DEBUG, "delta_cn was %f and set to zero. XCorr score is %f", delta_cn, current_xcorr);
        delta_cn = 0.0;
      }   
      if(fabs(delta_lcn) == numeric_limits<FLOAT_T>::infinity()){
        carp(CARP_DEBUG, "delta_lcn was %f and set to zero. XCorr score is %f", delta_lcn, current_xcorr);
        delta_lcn = 0.0;
      }   
      matches[idx]->setDeltaCn(delta_cn);
      matches[idx]->setDeltaLCn(delta_lcn);
    
    }   
  }

  return true;
}


/**********************************
 * match_collection get, set methods
 **********************************/

/**
 * \returns true if the match_collection only contains decoy matches,
 * else (all target or mixed) returns false.
 */
bool MatchCollection::isDecoy()
{
  return null_peptide_collection_;
}

/**
 * Try setting the match collection's charge.  Successful if the
 * current charge is 0 (i.e. hasn't yet been set) or if the current
 * charge is the same as the given value.  Otherwise, returns false
 *
 * \returns true if the match_collection's charge state was changed.
 */

bool MatchCollection::setZState(
  SpectrumZState& zstate ///< new zstate
  ) {

  if (getCharge() == 0) {
    zstate_ = zstate;
    return true;
  } else {
    //error
    carp(CARP_WARNING, "Cannot change the zstate of a match collection "
        "once it has been set.");
    return false;
  }
}

/**
 * Search the given database or index using shuffled peptides and the
 * spectrum/charge in the target psm match collection.  Add those
 * scores to the target psm match collection for use in weibull
 * parameter estimation but do not save the matches.
 */
void MatchCollection::addDecoyScores(
  Spectrum* spectrum, ///< search this spectrum
  SpectrumZState& zstate, ///< search spectrum at this charge state
  ModifiedPeptidesIterator* peptides ///< use these peptides to search
){

  // reuse these for scoring all matches
  int charge = zstate.getCharge();
  IonConstraint* ion_constraint = 
    IonConstraint::newIonConstraintSmart(XCORR, charge);
 
  IonSeries* ion_series = new IonSeries(ion_constraint, charge);  
  Scorer* scorer = new Scorer(XCORR);
  
  // for each peptide in the iterator
  while( peptides->hasNext() ){

    // get peptide and sequence
    Peptide* peptide = peptides->next();
    char* decoy_sequence = peptide->getSequence();
    MODIFIED_AA_T* modified_seq = peptide->getModifiedAASequence();

    // create the ion series for this peptide
    ion_series->update(decoy_sequence, modified_seq);
    ion_series->predictIons();

    // get the score
    FLOAT_T score = scorer->scoreSpectrumVIonSeries(spectrum, ion_series);

    // add to collection's list of xcorrs
    xcorrs_[num_xcorrs_] = score;
    num_xcorrs_++;

    // clean up
    free(decoy_sequence);
    free(modified_seq);
    delete peptide;
  } // next peptide

  IonConstraint::free(ion_constraint);
  delete ion_series;
  delete scorer;

}

// cheater functions for testing

void MatchCollection::forceScoredBy(SCORER_TYPE_T type){
  scored_type_[type] = true;
}

/**
 * Extract a given type of score into an array.  The array is
 * allocated here and must be freed by the caller.
 */
FLOAT_T* MatchCollection::extractScores(
  SCORER_TYPE_T       score_type ///< Type of score to extract.
)
{
  FLOAT_T* return_value = (FLOAT_T*)mycalloc(match_total_,
                                             sizeof(FLOAT_T));

  MatchIterator* match_iterator =
    new MatchIterator(this, XCORR, false);
  int idx = 0;
  while(match_iterator->hasNext()){
    Match* match = match_iterator->next();
    return_value[idx] = match->getScore(score_type);
    idx++;
  }
  delete match_iterator;

  return(return_value);
}

/**
 * Given a hash table that maps from a score to its q-value, assign
 * q-values to all of the matches in a given collection.
 */
void MatchCollection::assignQValues(
  const map<FLOAT_T, FLOAT_T>* score_to_qvalue_hash,
  SCORER_TYPE_T score_type
){

  // Iterate over the matches filling in the q-values
  MatchIterator* match_iterator = 
    new MatchIterator(this, score_type, false);

  while(match_iterator->hasNext()){
    Match* match = match_iterator->next();
    FLOAT_T score = match->getScore(score_type);

    // Retrieve the corresponding q-value.
    map<FLOAT_T, FLOAT_T>::const_iterator map_position 
      = score_to_qvalue_hash->find(score);
    if (map_position == score_to_qvalue_hash->end()) {
      carp(CARP_FATAL,
           "Cannot find q-value corresponding to score of %g.",
           score);
    }
    FLOAT_T qvalue = map_position->second;

    /* If we're given a base score, then store the q-value.  If we're
       given a q-value, then store the peptide-level q-value. */
    SCORER_TYPE_T derived_score_type = INVALID_SCORER_TYPE;
    switch (score_type) {
    case XCORR:
      derived_score_type = DECOY_XCORR_QVALUE;
      break;
    case DECOY_XCORR_QVALUE:
      derived_score_type = DECOY_XCORR_PEPTIDE_QVALUE;
      break;
    case EVALUE:
      derived_score_type = DECOY_EVALUE_QVALUE;
      break;
    case DECOY_EVALUE_QVALUE:
      derived_score_type = DECOY_EVALUE_PEPTIDE_QVALUE;
      break;
    case LOGP_BONF_WEIBULL_XCORR: 
      derived_score_type = LOGP_QVALUE_WEIBULL_XCORR;
      break;
    case LOGP_QVALUE_WEIBULL_XCORR:
      derived_score_type = LOGP_PEPTIDE_QVALUE_WEIBULL;
      break;
    case PERCOLATOR_SCORE:
      derived_score_type = PERCOLATOR_QVALUE;
      break;
    case PERCOLATOR_QVALUE:
      derived_score_type = PERCOLATOR_PEPTIDE_QVALUE;
      break;
    case QRANKER_SCORE:
      derived_score_type = QRANKER_QVALUE;
      break;
    case QRANKER_QVALUE:
      derived_score_type = QRANKER_PEPTIDE_QVALUE;
      break;
    case BARISTA_SCORE:
      derived_score_type = BARISTA_QVALUE;
      break;
    case BARISTA_QVALUE:
      derived_score_type = BARISTA_PEPTIDE_QVALUE;
      break;
    // Should never reach this point.
    case SP: 
    case LOGP_WEIBULL_XCORR: 
    case DECOY_XCORR_PEPTIDE_QVALUE:
    case LOGP_PEPTIDE_QVALUE_WEIBULL:
    case PERCOLATOR_PEPTIDE_QVALUE:
    case QRANKER_PEPTIDE_QVALUE:
    case QRANKER_PEP:
    case BARISTA_PEPTIDE_QVALUE:
    case BARISTA_PEP:
    case DECOY_XCORR_PEP:
    case LOGP_WEIBULL_PEP:
    case PERCOLATOR_PEP:
    case NUMBER_SCORER_TYPES:
    case INVALID_SCORER_TYPE:
      carp(CARP_FATAL, "Something is terribly wrong!");
    }

    match->setScore(derived_score_type, qvalue);
    scored_type_[derived_score_type] = true;

  }
  delete match_iterator;
}

/**
 * Given a hash table that maps from a score to its PEP, assign
 * PEPs to all of the matches in a given collection.
 */
void MatchCollection::assignPEPs(
    const map<FLOAT_T, FLOAT_T>* score_to_pep_hash,
    SCORER_TYPE_T score_type )
{
  // Iterate over the matches filling in the q-values
  MatchIterator* match_iterator = 
    new MatchIterator(this, score_type, false);

  while(match_iterator->hasNext()){
    Match* match = match_iterator->next();
    FLOAT_T score = match->getScore(score_type);

    // Retrieve the corresponding PEP.
    map<FLOAT_T, FLOAT_T>::const_iterator map_position 
      = score_to_pep_hash->find(score);
    if (map_position == score_to_pep_hash->end()) {
      carp(CARP_FATAL,
           "Cannot find q-value corresponding to score of %g.",
           score);
    }
    FLOAT_T qvalue = map_position->second;

    /* If we're given a base score, then store the q-value.  If we're
       given a q-value, then store the peptide-level q-value. */
    SCORER_TYPE_T derived_score_type = INVALID_SCORER_TYPE;
    switch (score_type) {
    case XCORR:
      derived_score_type = DECOY_XCORR_PEP;
      break;
    case EVALUE:
      derived_score_type = DECOY_EVALUE_PEP;
      break;
    case LOGP_BONF_WEIBULL_XCORR: 
      derived_score_type = LOGP_WEIBULL_PEP;
      break;
    case PERCOLATOR_SCORE:
      derived_score_type = PERCOLATOR_PEP;
      break;
    case QRANKER_SCORE:
      derived_score_type = QRANKER_PEP;
      break;
    case BARISTA_SCORE:
      derived_score_type = BARISTA_PEP;
      break;
    // Should never reach this point.
    case SP: 
    case LOGP_WEIBULL_XCORR: 
    case DECOY_XCORR_PEPTIDE_QVALUE:
    case DECOY_XCORR_QVALUE:
    case LOGP_QVALUE_WEIBULL_XCORR:
    case LOGP_PEPTIDE_QVALUE_WEIBULL:
    case PERCOLATOR_PEPTIDE_QVALUE:
    case QRANKER_PEPTIDE_QVALUE:
    case QRANKER_PEP:
    case QRANKER_QVALUE:
    case BARISTA_PEPTIDE_QVALUE:
    case BARISTA_PEP:
    case BARISTA_QVALUE:
    case DECOY_XCORR_PEP:
    case LOGP_WEIBULL_PEP:
    case PERCOLATOR_QVALUE:
    case PERCOLATOR_PEP:
    case NUMBER_SCORER_TYPES:
    case INVALID_SCORER_TYPE:
      carp(CARP_FATAL, "Something is terribly wrong!");
    }

    match->setScore(derived_score_type, qvalue);
    scored_type_[derived_score_type] = true;

  }
  delete match_iterator;

}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


