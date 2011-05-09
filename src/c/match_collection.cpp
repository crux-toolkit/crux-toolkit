/*********************************************************************//**
 * \file match_collection.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: 11/27 2006
 * \brief A set of peptide spectrum matches for one spectrum.
 *
 * Methods for creating and manipulating match_collections.   
 * Creating a match collection generates all matches (searches a
 * spectrum against a database.
 ****************************************************************************/
#include "match_collection.h"
#include <string>

#include "MatchFileReader.h"

using namespace std;

/* Private data types (structs) */

/**
 * \struct match_collection
 * \brief An object that contains a set of match objects.
 *
 * May contain matches for one spectrum or many spectra. 
 * 
 * 
 */
struct match_collection{
  MATCH_T* match[_MAX_NUMBER_PEPTIDES]; ///< array of match object
  int match_total;      ///< size of match array, may vary w/truncation
  int experiment_size;  ///< total matches before any truncation
  SpectrumZState zstate; ///< zstate of the associated spectrum
  BOOLEAN_T null_peptide_collection; ///< are the peptides shuffled
  BOOLEAN_T scored_type[NUMBER_SCORER_TYPES]; 
                        ///< TRUE if matches have been scored by the type
  SCORER_TYPE_T last_sorted; 
    ///< the last type by which it's been sorted ( -1 if unsorted)
  BOOLEAN_T iterator_lock; 
    ///< has an itterator been created? if TRUE can't manipulate matches

  // values used for various scoring functions.
  // TODO this should be moved to match
  FLOAT_T delta_cn; ///< the difference in top and second Xcorr scores
  FLOAT_T sp_scores_sum; ///< for getting mean, backward compatible
  FLOAT_T sp_scores_mean;  ///< the mean value of the scored peptides sp score
  FLOAT_T mu;// obsolete 
  ///< EVD parameter Xcorr(characteristic value of extreme value distribution)
  FLOAT_T l_value; // obsolete
  ///< EVD parameter Xcorr(decay constant of extreme value distribution)
  int top_fit_sp; // obsolete
  ///< The top ranked sp scored peptides to use as EXP_SP parameter estimation
  FLOAT_T base_score_sp; // obsolete
 ///< The lowest sp score within top_fit_sp, used as the base to rescale sp
  // Values for fitting the Weibull distribution
  FLOAT_T eta;  ///< The eta parameter for the Weibull distribution.
  FLOAT_T beta; ///< The beta parameter for the Weibull distribution.
  FLOAT_T shift; ///< The location parameter for the Weibull distribution.
  FLOAT_T correlation; ///< The correlation parameter for the Weibull distribution.
  // replace this ...
  MATCH_T* sample_matches[_PSM_SAMPLE_SIZE];
  int num_samples;  // the number of items in the above array
  // ...with this
  FLOAT_T xcorrs[_MAX_NUMBER_PEPTIDES]; ///< xcorrs to be used for weibull
  int num_xcorrs;

  // The following features (post_*) are only valid when
  // post_process_collection boolean is TRUE 
  BOOLEAN_T post_process_collection; 
  ///< Is this a post process match_collection?
  int post_protein_counter_size; 
  ///< the size of the protein counter array, usually the number of proteins in database
  int* post_protein_counter; 
  ///< the counter for how many each protein has matches other PSMs
  int* post_protein_peptide_counter; 
  ///< the counter for how many each unique peptides each protein has matches other PSMs
  HASH_T* post_hash; ///< hash table that keeps tracks of the peptides
  BOOLEAN_T post_scored_type_set; 
  ///< has the scored type been confirmed for the match collection,
  // set after the first match collection is extended
  MATCH_T* top_scoring_sp; ///< the match with Sp rank == 1 
};

/**
 *\struct match_iterator
 *\brief An object that iterates over the match objects in the
 * specified match_collection for the specified score type (SP, XCORR)
 */
struct match_iterator{
  MATCH_COLLECTION_T* match_collection; 
                            ///< the match collection to iterate -out
  SCORER_TYPE_T match_mode; ///< the current working score (SP, XCORR)
  int match_idx;            ///< current match to return
  int match_total;          ///< total_match_count
};

/**
 * \struct match_collection_iterator
 * \brief An object that iterates over the match_collection objects in
 * the specified directory of serialized match_collections 
 */
struct match_collection_iterator{
  DIR* working_directory; 
  ///< the working directory for the iterator to find match_collections
  char* directory_name; ///< the directory name in char
  DATABASE_T* database; ///< the database for which the match_collection
  int number_collections; 
  ///< the total number of match_collections in the directory (target+decoy)
  int collection_idx;  ///< the index of the current collection to return
  MATCH_COLLECTION_T* match_collection; ///< the match collection to return
  BOOLEAN_T is_another_collection; 
  ///< is there another match_collection to return?
  vector<bool>* cols_in_file; ///< which columns were in the target file
};

/******* Private function declarations, described in definintions below ***/
int add_unscored_peptides(
  MATCH_COLLECTION_T* match_collection, 
  Spectrum* spectrum, 
  SpectrumZState& charge, 
  MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator,
  BOOLEAN_T is_decoy
);

BOOLEAN_T score_matches_one_spectrum(
  SCORER_TYPE_T score_type, 
  MATCH_COLLECTION_T* match_collection,
  Spectrum* spectrum,
  int charge,
  BOOLEAN_T store_scores
  );

BOOLEAN_T populate_match_rank_match_collection(
 MATCH_COLLECTION_T* match_collection, 
 SCORER_TYPE_T score_type 
 );

void store_new_xcorrs(MATCH_COLLECTION_T* match_collection, 
                      int start_index,
                      BOOLEAN_T keep_matches);

void collapse_redundant_matches(MATCH_COLLECTION_T* matches);
void consolidate_matches(MATCH_T** matches, int start_idx, int end_idx);

BOOLEAN_T extend_match_collection_tab_delimited(
  MATCH_COLLECTION_T* match_collection, ///< match collection to extend -out
  DATABASE_T* database, ///< the database holding the peptides -in
  MatchFileReader& result_file   ///< the result file to parse PSMs -in
  );

BOOLEAN_T add_match_to_post_match_collection(
  MATCH_COLLECTION_T* match_collection, 
  MATCH_T* match 
  );

void update_protein_counters(
  MATCH_COLLECTION_T* match_collection, 
  PEPTIDE_T* peptide  
  );

BOOLEAN_T calculate_delta_cn(
  MATCH_COLLECTION_T* match_collection,
  COMMAND_T search_type = SEARCH_COMMAND
);

/********* end of function declarations *******************/


/**
 * \returns An (empty) match_collection object.
 */
MATCH_COLLECTION_T* allocate_match_collection()
{
  MATCH_COLLECTION_T* match_collection =
    (MATCH_COLLECTION_T*)mycalloc(1, sizeof(MATCH_COLLECTION_T));
    
  // loop over to set all score type to FALSE
  int score_type_idx = 0;
  for(; score_type_idx < NUMBER_SCORER_TYPES ; ++score_type_idx){
    match_collection->scored_type[score_type_idx] = FALSE;
  }
  
  // set last score to -1, thus nothing has been done yet
  match_collection->last_sorted = (SCORER_TYPE_T)-1;
  match_collection->iterator_lock = FALSE;
  match_collection->post_process_collection = FALSE;
  match_collection->null_peptide_collection = FALSE;
  
  return match_collection;
}

/**
 * /brief Free the memory allocated for a match collection
 * Deep free; each match is freed which, in turn, frees each spectrum
 * and peptide. 
 */
void free_match_collection(
  MATCH_COLLECTION_T* match_collection ///< the match collection to free -out
  )
{
  // decrement the pointer count in each match object
  while(match_collection->match_total > 0){
    --match_collection->match_total;
    free_match(match_collection->match[match_collection->match_total]);
    match_collection->match[match_collection->match_total] = NULL;
  }

  // and free the sample matches
  while(match_collection->num_samples > 0){
    --match_collection->num_samples;
   free_match(match_collection->sample_matches[match_collection->num_samples]);
    match_collection->sample_matches[match_collection->num_samples] = NULL;
  }
  
  // free post_process_collection specific memory
  if(match_collection->post_process_collection){
    // free protein counter
    free(match_collection->post_protein_counter);
    
    // free protein peptide counter
    free(match_collection->post_protein_peptide_counter);
  
    // free hash table
    free_hash(match_collection->post_hash);
  }

  if(match_collection->top_scoring_sp){
    free_match(match_collection->top_scoring_sp);
  }

  free(match_collection);
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
MATCH_COLLECTION_T* new_empty_match_collection(BOOLEAN_T is_decoy){

  MATCH_COLLECTION_T* match_collection = allocate_match_collection();
  int idx = 0;  

  // set member variables according to parameter.c 
  match_collection->match_total = 0;
  match_collection->experiment_size = 0; 
  match_collection->zstate = SpectrumZState();
  match_collection->null_peptide_collection = is_decoy;

  for(idx=0; idx < NUMBER_SCORER_TYPES ; idx++){
    match_collection->scored_type[idx] = FALSE;
  }
  match_collection->last_sorted = (SCORER_TYPE_T)-1;
  match_collection->iterator_lock = FALSE;
  match_collection->num_samples = 0;
  match_collection->num_xcorrs = 0;

  match_collection->post_hash = NULL;
  match_collection->top_scoring_sp = NULL;

  return match_collection;
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
int add_matches(
  MATCH_COLLECTION_T* matches, ///< add matches to this
  Spectrum* spectrum,  ///< compare peptides to this spectrum
  SpectrumZState& zstate,            ///< use this charge state for spectrum
  MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator, ///< use these peptides
  BOOLEAN_T is_decoy,     ///< are peptides to be shuffled
  BOOLEAN_T store_scores, ///< save scores for p-val estimation
  BOOLEAN_T do_sp_scoring, ///< start with SP scoring
  BOOLEAN_T filter_by_sp  ///< truncate matches based on Sp scores
){

  if( matches == NULL || peptide_iterator == NULL || spectrum == NULL ){
    carp(CARP_FATAL, "Cannot add matches to a collection when match " 
         "collection, spectrum and/or peptide iterator are NULL.");
  }

  //assert(matches->zstate == zstate);

  // generate a match for each peptide in the iterator, storing them
  // in the match collection
  int num_matches_added = add_unscored_peptides(matches, spectrum, zstate,
                                                peptide_iterator, is_decoy);

  if( num_matches_added == 0 ){
    return num_matches_added;
  }

  int xcorr_max_rank = get_int_parameter("psms-per-spectrum-reported");

  // optional Sp score on all candidate peptides
  if( do_sp_scoring ){
    score_matches_one_spectrum(SP, matches, spectrum, zstate.getCharge(),
                               FALSE); // don't store scores
    populate_match_rank_match_collection(matches, SP);
    if( filter_by_sp ){ // keep only high-ranking sp psms
      save_top_sp_match(matches);
      int sp_max_rank = get_int_parameter("max-rank-preliminary");
      truncate_match_collection(matches, 
                                sp_max_rank + 1, // extra for deltacn of last
                                SP);
      xcorr_max_rank = sp_max_rank;
    }
  }

  // main scoring
  score_matches_one_spectrum(XCORR, matches, spectrum, zstate.getCharge(), 
                             store_scores); 
  populate_match_rank_match_collection(matches, XCORR);
  truncate_match_collection(matches, 
                            xcorr_max_rank + 1,// extra for deltacn of last
                            XCORR);

  return num_matches_added;
}

/**
 * \brief Put all the matches from the source match collection in the
 * destination. Only copies the pointers of the matches so use with
 * caution. 
 * \returns The number of matches added.
 */
int merge_match_collections(MATCH_COLLECTION_T* source,
                            MATCH_COLLECTION_T* destination){
  if( source == NULL || destination == NULL ){
    carp(CARP_ERROR, "Cannot merge null match collections.");
  }
  carp(CARP_DETAILED_DEBUG, "Merging match collections.");

  // what is the index of the next insert position in destination
  int dest_idx = destination->match_total;

  // if these are the first being added to the destination, set the
  // scored_type
  if( dest_idx == 0 ){
    int type_idx = 0;
    for(type_idx = 0; type_idx < NUMBER_SCORER_TYPES; type_idx++){
      destination->scored_type[type_idx] = source->scored_type[type_idx];
    }
  }else{ // check that same types are scored
    int type_idx = 0;
    for(type_idx = 0; type_idx < NUMBER_SCORER_TYPES; type_idx++){
      if( destination->scored_type[type_idx] != source->scored_type[type_idx]){
        char type_str[SMALL_BUFFER];
        const char* dest_str = (destination->scored_type[type_idx]) ? "" : " not";
        const char* src_str = (source->scored_type[type_idx]) ? "" : " not";
        scorer_type_to_string((SCORER_TYPE_T)type_idx, type_str);
        carp(CARP_FATAL, "Cannot merge match collections scored for "
             "different types.  Trying to add matches%s scored for %s "
             "to matches%s scored for %s", 
             src_str, type_str, dest_str, type_str);
      }
    }
  }
  

  // make sure destination has room for more matches
  int src_num_matches = source->match_total;
  if( dest_idx + src_num_matches > _MAX_NUMBER_PEPTIDES ){
    carp(CARP_FATAL, "Cannot merge match collections, insufficient capacity "
         "in destnation collection.");
  }
  carp(CARP_DETAILED_DEBUG, "Merging %d matches into a collection of %d",
       src_num_matches, dest_idx );

  int src_idx = 0;
  // for each match in source
  for(src_idx = 0; src_idx < src_num_matches; src_idx++){
    MATCH_T* cur_match = source->match[src_idx];

    // copy pointer and add to destination
    increment_match_pointer_count(cur_match);
    destination->match[dest_idx] = cur_match;

    dest_idx++;
  }

  // update destination count
  destination->match_total += src_num_matches;
  destination->experiment_size += source->experiment_size;
  destination->last_sorted = (SCORER_TYPE_T)-1;  // unset any last-sorted flag

  return src_num_matches;
}


/**
 * \brief Store the xcorr for each psm that was added in this
 * iteration.  Assumes that the matches with scores needing storing
 * are between indexes start_index and match_collection->match_total.
 * The xcorrs will used for the weibull parameter estimations for
 * p-values.  If keep_matches == FALSE, the matches between indexes
 * start_index and match_collection->match_total will be deleted and
 * match_total will be updated.
 * 
 */
void store_new_xcorrs(
  MATCH_COLLECTION_T* match_collection, ///< source and destination of scores
  int start_index, ///< get first score from match at this index
  BOOLEAN_T keep_matches ///< FALSE=delete the matches after storing score
){

  if( match_collection == NULL ){
    carp(CARP_FATAL, "Cannot store scores of NULL match collection.");
  }

  int score_idx = match_collection->num_xcorrs;
  int psm_idx = start_index;

  carp(CARP_DETAILED_DEBUG, 
       "Adding to xcors[%i] scores from psm index %i to %i", 
       score_idx, psm_idx, match_collection->match_total);

  if( score_idx+(match_collection->match_total-psm_idx) 
      > _MAX_NUMBER_PEPTIDES ){
    carp(CARP_FATAL, "Too many xcorrs to store.");
  }

  for(psm_idx=start_index; psm_idx < match_collection->match_total; psm_idx++){
    FLOAT_T score = get_match_score( match_collection->match[psm_idx], XCORR);
    match_collection->xcorrs[score_idx] = score;
    score_idx++;

    if( keep_matches == FALSE ){
      free_match(match_collection->match[psm_idx]);
      match_collection->match[psm_idx] = NULL;
      match_collection->experiment_size -= 1;  // these should be decoys and 
                                               // we are not counting them
                                               
    }
  }

  match_collection->num_xcorrs = score_idx;
  if( keep_matches == FALSE ){
    match_collection->match_total = start_index; // where we started deleting
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
void collapse_redundant_matches(MATCH_COLLECTION_T* match_collection){
  if( match_collection == NULL ){
    carp(CARP_FATAL, "Cannot collapse matches from null collection.");
  }

  // must not be empty
  int match_total = match_collection->match_total;
  if( match_total == 0 ){
    return;
  }  

  carp(CARP_DETAILED_DEBUG, "Collapsing %i redundant matches.", match_total);

  // must be sorted by Sp or xcorr
  assert( (match_collection->last_sorted == SP) || 
          (match_collection->last_sorted == XCORR) );

  MATCH_T** matches = match_collection->match;
  int match_idx = 0;
  FLOAT_T cur_score = get_match_score(matches[match_idx], SP);

  // for entire list of matches
  while(match_idx < match_total-1){
    FLOAT_T next_score = get_match_score(matches[match_idx+1], SP);

    // find the index of the last match with the same score
    int cur_score_last_index = match_idx;
    
    while(next_score == cur_score && cur_score_last_index < match_total-2){
      cur_score_last_index++;
      next_score = get_match_score(matches[cur_score_last_index+1], SP);
    }
    // if the last two were equal, the last index was not incremented
    if( next_score == cur_score ){ cur_score_last_index++; }

    if( cur_score_last_index > match_idx ){
      consolidate_matches(matches, match_idx, cur_score_last_index);
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
       match_collection->match_total, opening_idx);
  // reset total number of matches in the collection
  match_collection->match_total = opening_idx;
  // remove duplicate peptides from the overall count
  int diff = match_total - opening_idx;
  carp(CARP_DETAILED_DEBUG, "Removing %i from total count %i",
       diff, match_collection->experiment_size);

  match_collection->experiment_size -= diff;
}

/**
 * \brief For a list of matches with the same scores, combine those
 * that are the same peptide and delete redundant matches.
 *
 * Since there may be different peptide sequences with the same score,
 * compare each match to the remaining matches.
 */
void consolidate_matches(MATCH_T** matches, int start_idx, int end_idx){

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
      get_match_mod_sequence_str_with_symbols(matches[cur_match_idx]);
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
        get_match_mod_sequence_str_with_symbols(matches[next_match_idx]);
      carp(CARP_DETAILED_DEBUG, "next seq is %s.", next_seq);

      if( strcmp(cur_seq, next_seq) == 0){
        carp(CARP_DETAILED_DEBUG, 
             "Seqs %s and %s match.  Consolidate match[%i] into match[%i].", 
             cur_seq, next_seq, next_match_idx, cur_match_idx);

        // add peptide src of next to cur
        merge_peptides_copy_src( get_match_peptide(matches[cur_match_idx]),
                        get_match_peptide(matches[next_match_idx]));
        // this frees the second peptide, so set what pointed to it to NULL
        //set_match_peptide(matches[next_match_idx], NULL);

        // delete match
        free_match(matches[next_match_idx]);
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
void sort_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to sort -out
  SCORER_TYPE_T score_type ///< the score type to sort by -in
  )
{
  carp(CARP_DEBUG, "Sorting match collection.");

  // check if we are allowed to alter match_collection
  if(match_collection->iterator_lock){
    carp(CARP_FATAL,
         "Cannot sort a match collection when a match iterator is already"
         " instantiated");
  }

  // Switch to the equivalent sort key.
  SCORER_TYPE_T sort_by = NUMBER_SCORER_TYPES; // Initialize to nonsense.
  int (*compare_match_function)(const void*, const void*) 
    = (QSORT_COMPARE_METHOD)compare_match_sp;
  switch (score_type) {
  case SP: 
    carp(CARP_DEBUG, "Sorting match collection by Sp.");
    sort_by = SP;
    compare_match_function = (QSORT_COMPARE_METHOD)compare_match_sp;
    break;

  case XCORR:
  case DECOY_XCORR_QVALUE:
  case LOGP_WEIBULL_XCORR: 
  case DECOY_XCORR_PEPTIDE_QVALUE:
    carp(CARP_DEBUG, "Sorting match collection by XCorr.");
    sort_by = XCORR;
    compare_match_function = (QSORT_COMPARE_METHOD)compare_match_xcorr;
    break;

  case LOGP_BONF_WEIBULL_XCORR: 
  case LOGP_QVALUE_WEIBULL_XCORR:
  case LOGP_PEPTIDE_QVALUE_WEIBULL:
    carp(CARP_DEBUG, "Sorting match collection by p-value.");
    sort_by = LOGP_BONF_WEIBULL_XCORR;
    compare_match_function = (QSORT_COMPARE_METHOD)compare_match_p_value;
    break;

  case PERCOLATOR_SCORE:
    carp(CARP_INFO, "Sorting match collection by Percolator score.");
    sort_by = PERCOLATOR_SCORE;
    compare_match_function = (QSORT_COMPARE_METHOD)compare_match_percolator_score;
    break;

  case PERCOLATOR_QVALUE:
  case PERCOLATOR_PEPTIDE_QVALUE:
    carp(CARP_DEBUG, "Sorting match collection by Percolator q-value.");
    sort_by = PERCOLATOR_QVALUE;
    compare_match_function = (QSORT_COMPARE_METHOD)compare_match_percolator_qvalue;
    break;

  case QRANKER_SCORE:
    carp(CARP_DEBUG, "Sorting match collection by Q-ranker score.");
    sort_by = QRANKER_SCORE;
    compare_match_function = (QSORT_COMPARE_METHOD)compare_match_qranker_score;
    break;

  case QRANKER_QVALUE:
  case QRANKER_PEPTIDE_QVALUE:
    carp(CARP_DEBUG, "Sorting match collection by Q-ranker q-value.");
    compare_match_function = (QSORT_COMPARE_METHOD)compare_match_qranker_qvalue;
    sort_by = QRANKER_QVALUE;
    break;

  // Should never reach this point.
  case NUMBER_SCORER_TYPES:
  case INVALID_SCORER_TYPE:
    carp(CARP_FATAL, "Something is terribly wrong in the sorting code!");
  }

  // Don't sort if it's already sorted.
  if (match_collection->last_sorted == sort_by) {
    return;
  }

  // Do the sort.
  qsort_match(match_collection->match,
              match_collection->match_total,
              compare_match_function);
  match_collection->last_sorted = sort_by;
}

/**
 * \brief Sort a match_collection by the given score type, grouping
 * matches by spectrum (if multiple spectra are present).
 */
void spectrum_sort_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< match collection to sort -out
  SCORER_TYPE_T score_type ///< the score type to sort by -in
  ){

  // check if we are allowed to alter match_collection
  if(match_collection->iterator_lock){
    carp(CARP_FATAL,
         "Cannot alter match_collection when a match iterator is already"
         " instantiated");
  }

  switch(score_type){
  case SP: 
    qsort_match(match_collection->match, match_collection->match_total,
                (QSORT_COMPARE_METHOD)compare_match_spectrum_sp);
    match_collection->last_sorted = SP;
    break;

  case XCORR:
  case LOGP_WEIBULL_XCORR: 
  case LOGP_BONF_WEIBULL_XCORR: 
  case LOGP_QVALUE_WEIBULL_XCORR: 
  case LOGP_PEPTIDE_QVALUE_WEIBULL:
  case DECOY_XCORR_QVALUE:
  case DECOY_XCORR_PEPTIDE_QVALUE:
    /* If we are sorting on a per-spectrum basis, then the xcorr is
       good enough, even in the presence of p-values. */
    qsort_match(match_collection->match, match_collection->match_total,
                (QSORT_COMPARE_METHOD)compare_match_spectrum_xcorr);
    match_collection->last_sorted = XCORR;
    break;

  case PERCOLATOR_SCORE:
    qsort_match(match_collection->match, match_collection->match_total,
                (QSORT_COMPARE_METHOD)compare_match_spectrum_percolator_score);
    match_collection->last_sorted = PERCOLATOR_SCORE;
    break;

  case PERCOLATOR_QVALUE:
  case PERCOLATOR_PEPTIDE_QVALUE:
    qsort_match(match_collection->match, match_collection->match_total,
                (QSORT_COMPARE_METHOD)compare_match_spectrum_percolator_qvalue);
    match_collection->last_sorted = PERCOLATOR_QVALUE;
    break;

  case QRANKER_SCORE:
    qsort_match(match_collection->match, match_collection->match_total,
                (QSORT_COMPARE_METHOD)compare_match_spectrum_qranker_score);
    match_collection->last_sorted = QRANKER_SCORE;
    break;

  case QRANKER_QVALUE:
  case QRANKER_PEPTIDE_QVALUE:
    qsort_match(match_collection->match, match_collection->match_total,
                (QSORT_COMPARE_METHOD)compare_match_spectrum_qranker_qvalue);
    match_collection->last_sorted = QRANKER_QVALUE;
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
void truncate_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< match collection to truncate -out
  int max_rank,     ///< max rank of matches to keep from SP -in
  SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
  )
{
  carp(CARP_DETAILED_DEBUG, "Truncating match collection to rank %d.", max_rank);
  if (match_collection == NULL || match_collection->match_total == 0){
    carp(CARP_DETAILED_DEBUG, "No matches in collection, so not truncating");
    return;
  }

  // sort match collection by score type
  sort_match_collection(match_collection, score_type);

  // Free high ranking matches
  int highest_index = match_collection->match_total - 1;
  MATCH_T** matches = match_collection->match;
  int cur_rank = get_match_rank(matches[highest_index], score_type);

  while( cur_rank > max_rank ){
    free_match(matches[highest_index]);
    highest_index--;
    cur_rank = get_match_rank(matches[highest_index], score_type);
  }
  match_collection->match_total = highest_index + 1;

  carp(CARP_DETAILED_DEBUG, "Truncated collection now has %d matches.", 
       match_collection->match_total);

}

/**
 * Assigns a rank for the given score type to each match.  First sorts
 * by the score type (if not already sorted).  Overwrites any existing
 * rank values, so it can be performed on a collection with matches
 * newly added to previously ranked matches.  Rank 1 is highest
 * score.  Matches with the same score will be given the same rank.
 *
 * \returns TRUE, if populates the match rank in the match collection
 */
BOOLEAN_T populate_match_rank_match_collection(
 MATCH_COLLECTION_T* match_collection, ///< match collection to rank -out
 SCORER_TYPE_T score_type ///< score type (SP, XCORR) by which to rank -in
 )
{
  carp(CARP_DETAILED_DEBUG, "Ranking matches by %i.", score_type);
  carp(CARP_DETAILED_DEBUG, "Collection currently ranked by %d", match_collection->last_sorted);
  // check if the match collection is in the correct sorted order
  if(match_collection->last_sorted != score_type){
    // sort match collection by score type
    carp(CARP_DETAILED_DEBUG, "Sorting by score_type %i", score_type);
    sort_match_collection(match_collection, score_type);
  }

  // set match rank for all match objects that have been scored for
  // this type
  int match_index;
  int cur_rank = 0;
  FLOAT_T cur_score = NOT_SCORED;
  for(match_index=0; match_index<match_collection->match_total; ++match_index){
    MATCH_T* cur_match = match_collection->match[match_index];
    FLOAT_T this_score = get_match_score(cur_match, score_type);
    
    if( NOT_SCORED == get_match_score(cur_match, score_type) ){
      char* seq = get_match_mod_sequence_str_with_masses(cur_match, FALSE);
      carp(CARP_WARNING, 
           "PSM spectrum %i charge %i sequence %s was NOT scored for type %i",
           (get_match_spectrum(cur_match))->getFirstScan(),
           get_match_charge, seq,
           (int)score_type);
      free(seq);
    }

    // does this match have a higher score?
    if( this_score != cur_score ){
      cur_score = this_score;
      cur_rank++;
    }

    //    set_match_rank( cur_match, score_type, match_index+1);
    set_match_rank( cur_match, score_type, cur_rank);

    carp(CARP_DETAILED_DEBUG, "Match rank %i, score %f", cur_rank, cur_score);
  }
  
  return TRUE;
}

/**
 * Keep track of the top-scoring Sp match.  It should be printed to
 * the sqt file even if its XCORR rank is not high enough to be
 * printed.  Requires that ranks have been set for Sp.
 *
 */
void save_top_sp_match(MATCH_COLLECTION_T* match_collection){

  if( match_collection == NULL ){
    carp(CARP_FATAL, "Cannot set top Sp match for NULL match collection.");
  }

  assert(match_collection->match_total > 0);
  MATCH_T* cur_rank_one_match = match_collection->match[0];

  // confirm that matches are sorted and ranks are set
  if( match_collection->last_sorted != SP || 
      get_match_rank(cur_rank_one_match, SP) != 1 ){
    populate_match_rank_match_collection(match_collection, SP);
  }

  // if no top sp yet, set it
  if( match_collection->top_scoring_sp == NULL ){
    match_collection->top_scoring_sp = cur_rank_one_match;
    increment_match_pointer_count(cur_rank_one_match);
    return;
  }

  // otherwise, see if the current top-ranked match has a higher score
  // the rank of top_scoring_sp should have a new rank
  if( get_match_rank(match_collection->top_scoring_sp, SP) > 1 ){
    free_match(match_collection->top_scoring_sp);
    match_collection->top_scoring_sp = cur_rank_one_match;
    increment_match_pointer_count(cur_rank_one_match);
  }
}

/**
 * Create a new match_collection by randomly sampling matches 
 * from match_collection upto count number of matches
 * Must not free the matches
 * \returns a new match_collection of randomly sampled matches 
 */
MATCH_COLLECTION_T* random_sample_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to sample -out
  int count_max ///< the number of matches to randomly select -in
  )
{
  int count_idx = 0;
  int match_idx = 0;
  int score_type_idx = 0;
  MATCH_COLLECTION_T* sample_collection = allocate_match_collection();
  srandom(time(NULL));

  // make sure we don't sample more than the matches in the match collection
  if (count_max >= match_collection->match_total){
    free_match_collection(sample_collection);
    return match_collection;
  }

  // ranomly select matches upto count_max
  for(; count_idx < count_max; ++count_idx){
    match_idx = (int)((double)random()/((double)RAND_MAX + (double)1)) 
      * match_collection->match_total;
    
    // match_idx = random() % match_collection->match_total;
    sample_collection->match[count_idx] = match_collection->match[match_idx];
    // increment pointer count of the match object 
    increment_match_pointer_count(sample_collection->match[count_idx]);
  }
  
  // set total number of matches sampled
  sample_collection->match_total = count_idx;

  sample_collection->experiment_size = match_collection->experiment_size;

  // set scored types in the sampled matches
  for(; score_type_idx < NUMBER_SCORER_TYPES ;  ++score_type_idx){
    sample_collection->scored_type[score_type_idx] 
      = match_collection->scored_type[score_type_idx];
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
void constraint_function(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to estimate evd parameters -in
  SCORER_TYPE_T score_type, ///< score_type to estimate EVD distribution -in
  FLOAT_T l_value,  ///< L value -in
  FLOAT_T* function,  ///< the output function value -out
  FLOAT_T* derivative,  ///< the output derivative value -out
  FLOAT_T* exponential_sum ///< the final exponential array sum -out
  )
{
  int idx = 0;
  FLOAT_T* exponential = (FLOAT_T*)mycalloc(match_collection->match_total, sizeof(FLOAT_T));
  FLOAT_T numerator = 0;
  FLOAT_T second_numerator = 0;
  FLOAT_T score = 0;
  FLOAT_T denominator = 0;
  FLOAT_T score_sum = 0;
  MATCH_T** matches = match_collection->match;

  // iterate over the matches to calculate numerator, exponential value, denominator
  for(; idx < match_collection->match_total; ++idx){
    score = get_match_score(matches[idx], score_type);
    exponential[idx] = exp(-l_value * score);
    numerator += (exponential[idx] * score);
    denominator += exponential[idx];
    score_sum += score;
    second_numerator += (score * score * exponential[idx]);
  }

  // assign function value
  *function = (1.0 / l_value) - (score_sum / match_collection->match_total) 
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
 *\returns TRUE, if successfully calculates the Weibull parameters
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
 * \returns TRUE if their are enough xcorrs for estimating Weibull
 * parameters or FALSE if not.
 */
BOOLEAN_T has_enough_weibull_points(
  MATCH_COLLECTION_T* match_collection
){
  return (match_collection->num_xcorrs >= MIN_WEIBULL_MATCHES );
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
BOOLEAN_T estimate_weibull_parameters_from_xcorrs(
  MATCH_COLLECTION_T* match_collection, 
  Spectrum* spectrum,
  int charge
  ){

  if( match_collection == NULL || spectrum == NULL ){
    carp(CARP_ERROR, "Cannot estimate parameters from null inputs.");
    return FALSE;
  }

  // check that we have the minimum number of matches
  FLOAT_T* scores = match_collection->xcorrs;
  int num_scores = match_collection->num_xcorrs;
  if( num_scores < MIN_WEIBULL_MATCHES ){
    carp(CARP_DETAILED_DEBUG, "Too few psms (%i) to estimate "
         "p-value parameters for spectrum %i, charge %i",
         num_scores, spectrum->getFirstScan(), charge);
    // set eta, beta, and shift to something???
    return FALSE;
  }

  // reverse sort the scores
  sort(scores, scores + num_scores, compareDescending());

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
      &(match_collection->eta), &(match_collection->beta),
      &(match_collection->shift), &(match_collection->correlation));
  carp(CARP_DEBUG, 
      "Corr: %.6f  Eta: %.6f  Beta: %.6f  Shift: %.6f", 
       match_collection->correlation, match_collection->eta, 
       match_collection->beta, match_collection->shift);
  
  return TRUE;
}

// TODO (BF 16-mar-09): use this instead of score_peptides
/**
 * \brief Add all peptides from iterator to match collection.
 * Additional matches will not be scored for any type.
 * \returns The number of peptides added.
 */
int add_unscored_peptides(
  MATCH_COLLECTION_T* match_collection, 
  Spectrum* spectrum, 
  SpectrumZState& zstate, 
  MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator,
  BOOLEAN_T is_decoy
){

  if( match_collection == NULL || spectrum == NULL 
      || peptide_iterator == NULL ){
    carp(CARP_FATAL, "Cannot score peptides with NULL inputs.");
  }
  carp(CARP_DETAILED_DEBUG, "Adding decoy peptides to match collection? %i", 
       is_decoy);

  int starting_number_of_psms = match_collection->match_total;

  while( modified_peptides_iterator_has_next(peptide_iterator)){
    // get peptide
    PEPTIDE_T* peptide = modified_peptides_iterator_next(peptide_iterator);

    // create a match
    MATCH_T* match = new_match();

    // set match fields
    set_match_peptide(match, peptide);
    set_match_spectrum(match, spectrum);
    set_match_zstate(match, zstate);
    set_match_null_peptide(match, is_decoy);

    // add to match collection
    if(match_collection->match_total >= _MAX_NUMBER_PEPTIDES){
      carp(CARP_ERROR, "peptide count of %i exceeds max match limit: %d", 
          match_collection->match_total, _MAX_NUMBER_PEPTIDES);

      return FALSE;
    }

    match_collection->match[match_collection->match_total] = match;
    match_collection->match_total++;

  }// next peptide

  int matches_added = match_collection->match_total - starting_number_of_psms;
  match_collection->experiment_size += matches_added;

  // matches are no longer correctly sorted
  match_collection->last_sorted = (SCORER_TYPE_T)-1; // unsorted
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
 * \returns TRUE, if matches are successfully scored.
 */
BOOLEAN_T score_matches_one_spectrum(
  SCORER_TYPE_T score_type, 
  MATCH_COLLECTION_T* match_collection,
  Spectrum* spectrum,
  int charge,
  BOOLEAN_T store_scores
  ){

  if( match_collection == NULL || spectrum == NULL ){
    carp(CARP_ERROR, "Cannot score matches in a NULL match collection.");
    return FALSE;
  }
  
  MATCH_T** matches = match_collection->match;
  int num_matches = match_collection->match_total;

  char type_str[64];
  scorer_type_to_string(score_type, type_str);
  carp(CARP_DETAILED_DEBUG, "Scoring matches for %s", type_str);

  // create ion constraint
  IonConstraint* ion_constraint = 
    IonConstraint::newIonConstraintSmart(score_type, charge);

  // create scorer
  SCORER_T* scorer = new_scorer(score_type);

  // create a generic ion_series that will be reused for each peptide sequence
  IonSeries* ion_series = new IonSeries(ion_constraint, charge);  
  
  // score all matches
  int match_idx;

  for(match_idx = 0; match_idx < num_matches; match_idx++){

    MATCH_T* match = matches[match_idx];
    assert( match != NULL );

    // skip it if it's already been scored
    if( NOT_SCORED != get_match_score(match, score_type)){
      continue;
    }

    // make sure it's the same spec and charge
    assert( spectrum == get_match_spectrum(match));
    assert( charge == get_match_charge(match));
    char* sequence = get_match_sequence(match);
    MODIFIED_AA_T* modified_sequence = get_match_mod_sequence(match);

    // create ion series for this peptide
    ion_series->update(sequence, modified_sequence);
    ion_series->predictIons();

    // get the score
    FLOAT_T score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);

    // set score in match
    set_match_score(match, score_type, score);
    if( score_type == SP ){
      set_match_b_y_ion_info(match, scorer);
    }

    // save score in collection
    if( store_scores ){
      match_collection->xcorrs[match_collection->num_xcorrs] = score;
      match_collection->num_xcorrs++; 
    }

    IF_CARP_DETAILED_DEBUG(
      char* mod_seq = 
      modified_aa_string_to_string_with_masses(modified_sequence,
                                               strlen(sequence),
                                               FALSE);
      carp(CARP_DETAILED_DEBUG, "Second score %f for %s (null:%i)",
           score, mod_seq,get_match_null_peptide(match));
      free(mod_seq);
    )
    free(sequence);
    free(modified_sequence);
  }// next match

  // set the match_collection as having been scored
  match_collection->scored_type[score_type] = TRUE;

  // clean up
  IonConstraint::free(ion_constraint);
  delete ion_series;
  free_scorer(scorer);
  return TRUE;
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
 * \returns TRUE if p-values successfully computed for all matches,
 * else FALSE.
 */
// FIXME (BF 8-Dec-2008): create new score-type P_VALUE to replace LOG...XCORR
BOOLEAN_T compute_p_values(
  MATCH_COLLECTION_T* match_collection,
  FILE* output_pvalue_file ///< If non-NULL, file for storing p-values -in
  ){

  if(match_collection == NULL){
    carp(CARP_ERROR, "Cannot compute p-values for NULL match collection.");
    return FALSE;
  }

  int scan_number = 
    (get_match_spectrum(match_collection->match[0]))->getFirstScan();
  carp(CARP_DEBUG, "Computing p-values for %s spec %d charge %d "
       "with eta %f beta %f shift %f",
       (match_collection->null_peptide_collection) ? "decoy" : "target",
       scan_number,
       match_collection->zstate.getCharge(),
       match_collection->eta, match_collection->beta, match_collection->shift);

  SCORER_TYPE_T main_score = get_scorer_type_parameter("score-type");

  // check that the matches have been scored
  if(!match_collection->scored_type[main_score]){
    char type_str[64];
    scorer_type_to_string(main_score, type_str);
    carp(CARP_FATAL, 
         "Match collection was not scored by %s prior to computing p-values.",
         type_str);
  }

  // Print separator in the decoy p-value file.
  if (output_pvalue_file) {
    fprintf(output_pvalue_file, "# scan: %d charge: %d candidates: %d\n", 
            scan_number, match_collection->zstate.getCharge(),
            match_collection->experiment_size);
    fprintf(output_pvalue_file, 
            "# eta: %g beta: %g shift: %g correlation: %g\n",
            match_collection->eta, 
            match_collection->beta,
            match_collection->shift,
            match_collection->correlation);
  }

  // iterate over all matches 
  int match_idx =0;
  for(match_idx=0; match_idx < match_collection->match_total; match_idx++){
    MATCH_T* cur_match = match_collection->match[match_idx];

    // Get the Weibull p-value.
    double pvalue = compute_weibull_pvalue(get_match_score(cur_match, 
                                                           main_score),
                                           match_collection->eta, 
                                           match_collection->beta,
                                           match_collection->shift);

    // Print the pvalue, if requested
    if (output_pvalue_file) {
      fprintf(output_pvalue_file, "%g\n", pvalue);
    }

    // Apply the Bonferroni correction.
    pvalue = bonferroni_correction(pvalue, match_collection->experiment_size);

    // set pvalue in match
    set_match_score(cur_match, LOGP_BONF_WEIBULL_XCORR, -log(pvalue));
    //#endif

  }// next match

  carp(CARP_DETAILED_DEBUG, "Computed p-values for %d PSMs.", match_idx);
  populate_match_rank_match_collection(match_collection, XCORR);

  // mark p-values as having been scored
  match_collection->scored_type[LOGP_BONF_WEIBULL_XCORR] = TRUE;
  return TRUE;
}

/**
 * \brief Use the matches collected from all spectra to compute FDR
 * and q_values from the ranked list of target and decoy scores.
 * Assumes the match_collection has an appropriate number of
 * target/decoy matches per spectrum (e.g. one target and one decoy
 * per spec).
 * \returns TRUE if q-values successfully computed, else FALSE.
 */
BOOLEAN_T compute_decoy_q_values(
 MATCH_COLLECTION_T* match_collection ///< Set of PSMs to be sorted.
){

  if( match_collection == NULL ){
    carp(CARP_ERROR, "Cannot compute q-values for null match collection.");
    return FALSE;
  }

  carp(CARP_DEBUG, "Computing decoy q-values.");

  // sort by score
  sort_match_collection(match_collection, XCORR);

  // compute FDR from a running total of number targets/decoys
  // FDR = #decoys / #targets
  FLOAT_T num_targets = 0;
  FLOAT_T num_decoys = 0;
  int match_idx = 0;
  for(match_idx = 0; match_idx < match_collection->match_total; match_idx++){
    MATCH_T* cur_match = match_collection->match[match_idx];

    if ( get_match_null_peptide(cur_match) == TRUE ){
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
    set_match_score(cur_match, DECOY_XCORR_QVALUE, score);
    carp(CARP_DETAILED_DEBUG, 
         "match %i xcorr or pval %f num targets %i, num decoys %i, score %f",
         match_idx, get_match_score(cur_match, XCORR), 
         (int)num_targets, (int)num_decoys, score);
  }

  // compute q-value: go through list in reverse and use min FDR seen
  FLOAT_T min_fdr = 1.0;
  for(match_idx = match_collection->match_total-1; match_idx >= 0; match_idx--){
    MATCH_T* cur_match = match_collection->match[match_idx];
    FLOAT_T cur_fdr = get_match_score(cur_match, DECOY_XCORR_QVALUE);
    if( cur_fdr == P_VALUE_NA ){ continue; }

    if( cur_fdr < min_fdr ){
      min_fdr = cur_fdr;
    }

    set_match_score(cur_match, DECOY_XCORR_QVALUE, min_fdr);
    carp(CARP_DETAILED_DEBUG, 
         "match %i cur fdr %f min fdr %f is decoy %i",
         match_idx, cur_fdr, min_fdr, get_match_null_peptide(cur_match) );
  }

  match_collection->scored_type[DECOY_XCORR_QVALUE] = TRUE;
  return TRUE;
}


/**
 * match_collection get, set method
 */

/**
 *\returns TRUE, if the match collection has been scored by score_type
 */
BOOLEAN_T get_match_collection_scored_type(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to iterate -in
  SCORER_TYPE_T score_type ///< the score_type (MATCH_SP, MATCH_XCORR) -in
  )
{
  return match_collection->scored_type[score_type];
}

/**
 * sets the score_type to value
 */
void set_match_collection_scored_type(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to iterate -in
  SCORER_TYPE_T score_type, ///< the score_type (MATCH_SP, MATCH_XCORR) -in
  BOOLEAN_T value
  )
{
  match_collection->scored_type[score_type] = value;
}

/**
 *\returns TRUE, if there is a  match_iterators instantiated by match collection 
 */
BOOLEAN_T get_match_collection_iterator_lock(
  MATCH_COLLECTION_T* match_collection ///< working match collection -in
  )
{
  return match_collection->iterator_lock;
}

/**
 *\returns the total match objects avaliable in current match_collection
 */
int get_match_collection_match_total(
  MATCH_COLLECTION_T* match_collection ///< working match collection -in
  )
{
  return match_collection->match_total;
}

/**
 *\returns the total peptides searched in the experiment in match_collection
 */
int get_match_collection_experimental_size(
  MATCH_COLLECTION_T* match_collection ///< working match collection -in
  )
{
  return match_collection->experiment_size;
}

/**
 *\returns the top peptide count used in the logp_exp_sp in match_collection
 */
int get_match_collection_top_fit_sp(
  MATCH_COLLECTION_T* match_collection ///< working match collection -in
  )
{
  return match_collection->top_fit_sp;
}

/**
 *\returns the charge of the spectrum that the match collection was created
 */
int get_match_collection_charge(
  MATCH_COLLECTION_T* match_collection ///< working match collection -in
  )
{
  return match_collection->zstate.getCharge();
}

/**
 * Must have been scored by Xcorr, returns error if not scored by Xcorr
 *\returns the delta cn value(difference in top and second ranked Xcorr values)
 */
FLOAT_T get_match_collection_delta_cn(
  MATCH_COLLECTION_T* match_collection ///< working match collection -in
  )
{
  // Check if xcorr value has been scored, thus delta cn value is valid
  if(match_collection->scored_type[XCORR]){
    return match_collection->delta_cn;
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
void transfer_match_collection_weibull(
  MATCH_COLLECTION_T* from_collection,
  MATCH_COLLECTION_T* to_collection
  ){
  to_collection->eta = from_collection->eta;
  to_collection->beta = from_collection->beta;
  to_collection->shift = from_collection->shift;
  to_collection->correlation = from_collection->correlation;
}

/**
 * \brief Prints out the pepxml header to the output stream
 * passed in as a parameter.
 */

void print_xml_header(
  FILE* output
  ){
  if (output == NULL ){
    return;
  }
  time_t hold_time;
  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  char* enz_str = enzyme_type_to_string(enzyme);
  char* database = get_string_parameter("protein database");
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
  BOOLEAN_T missed_cleavage = get_boolean_parameter("missed-cleavages");
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




  BOOLEAN_T use_index = is_directory(database);
  if( use_index == TRUE ){
    char* fasta_name  = get_index_binary_fasta_name(database);
    free(database);
    database = fasta_name;
  }
  #if DARWIN
  char path_buffer[PATH_MAX];
  char* absolute_database_path =  realpath(database, path_buffer);
  #else
  char* absolute_database_path =  realpath(database, NULL);
  #endif
  free(database);

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
              mass,
              mod,
              "N" // N if static modification
              );      
    }
  }
  
  // variable amino acid modifications
  AA_MOD_T** mod_list = NULL;
  int num_mods = get_all_aa_mod_list(&mod_list);
  for (int mod_idx = 0; mod_idx < num_mods; mod_idx++){
    FLOAT_T mass = aa_mod_get_mass_change(mod_list[mod_idx]);
    
    BOOLEAN_T* aas_modified = aa_mod_get_aa_list(mod_list[mod_idx]);
    for (int aa_idx = 0; aa_idx < AA_LIST_LENGTH; aa_idx++){
      if (aas_modified[aa_idx] == TRUE ){
        int aa = (aa_idx+'A');
        FLOAT_T original_mass = get_mass_amino_acid(aa , isotopic_type);
        FLOAT_T mass_dif = mass - original_mass;
        fprintf(output, "<aminoacid_modification aminoacid=\"%c\" mass=\"%f\" "
                "massdiff=\"%f\" variable=\"%s\" />\n",
                aa,
                mass,
                mass_dif,
                "Y" // Y if variable modification
                );    

      }
    }

  }
  print_parameters_xml(output);
  
  fprintf(output, "</search_summary>\n");
  
}


void print_sqt_header(
 FILE* output, 
 const char* type, 
 int num_proteins, 
 BOOLEAN_T is_analysis){  // for analyze-matches look at algorithm param
  if( output == NULL ){
    return;
  }

  time_t hold_time;
  hold_time = time(0);

  BOOLEAN_T decoy = FALSE;
  if( strcmp(type, "decoy") == 0 ){
    decoy = TRUE;
  }

  fprintf(output, "H\tSQTGenerator Crux\n");
  fprintf(output, "H\tSQTGeneratorVersion 1.0\n");
  fprintf(output, "H\tComment Crux was written by...\n");
  fprintf(output, "H\tComment ref...\n");
  fprintf(output, "H\tStartTime\t%s", ctime(&hold_time));
  fprintf(output, "H\tEndTime                               \n");

  char* database = get_string_parameter("protein database");
  BOOLEAN_T use_index = is_directory(database);

  if( use_index == TRUE ){
    char* fasta_name  = get_index_binary_fasta_name(database);
    free(database);
    database = fasta_name;
  }
  fprintf(output, "H\tDatabase\t%s\n", database);
  free(database);

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
  scorer_type_to_string(score, temp_str);
  fprintf(output, "H\tComment\tpreliminary algorithm %s\n", temp_str);

  score = get_scorer_type_parameter("score-type");
  scorer_type_to_string(score, temp_str);
  fprintf(output, "H\tComment\tfinal algorithm %s\n", temp_str);

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

  //PEPTIDE_TYPE_T cleavages = get_peptide_type_parameter("cleavages");
  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  DIGEST_T digestion = get_digest_type_parameter("digestion");
  //peptide_type_to_string(cleavages, temp_str);
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
  fprintf(output, "H\tLine fields: S, scan number, scan number,"
          "charge, 0, precursor mass, 0, 0, number of matches\n");

  // fancy logic for printing the scores. see match.c:print_match_sqt

  SCORER_TYPE_T main_score = get_scorer_type_parameter("score-type");
  SCORER_TYPE_T other_score = get_scorer_type_parameter("prelim-score-type");
  ALGORITHM_TYPE_T analysis_score = get_algorithm_type_parameter("algorithm");
  BOOLEAN_T pvalues = get_boolean_parameter("compute-p-values");
  if( is_analysis == TRUE && analysis_score == PERCOLATOR_ALGORITHM){
    main_score = PERCOLATOR_SCORE; 
    other_score = PERCOLATOR_QVALUE;
  }else if( is_analysis == TRUE && analysis_score == QRANKER_ALGORITHM ){
    main_score = QRANKER_SCORE; 
    other_score = QRANKER_QVALUE;
  }else if( is_analysis == TRUE && analysis_score == QVALUE_ALGORITHM ){
    main_score = LOGP_QVALUE_WEIBULL_XCORR;
  }else if( pvalues == TRUE ){
    main_score = LOGP_BONF_WEIBULL_XCORR;
  }

  char main_score_str[64];
  scorer_type_to_string(main_score, main_score_str);
  char other_score_str[64];
  scorer_type_to_string(other_score, other_score_str);

  // ranks are always xcorr and sp
  // main/other scores from search are...xcorr/sp (OK as is)
  // ...p-val/xcorr
  if( main_score == LOGP_BONF_WEIBULL_XCORR ){
    strcpy(main_score_str, "-log(p-value)");
    strcpy(other_score_str, "xcorr");
  }// main/other scores from analyze are perc/q-val (OK as is)
   // q-val/xcorr
  if( main_score == LOGP_QVALUE_WEIBULL_XCORR ){
    strcpy(main_score_str, "q-value");  // to be changed to curve-fit-q-value
    strcpy(other_score_str, "xcorr");
  }

  fprintf(output, "H\tLine fields: M, rank by xcorr score, rank by sp score, "
          "peptide mass, deltaCn, %s score, %s score, number ions matched, "
          "total ions compared, sequence\n", main_score_str, other_score_str);
}

/**
 * Print the header line for a tab-delimited file.
 */
void print_tab_header(FILE* output){

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
void print_xml_footer(FILE* output){
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
void print_matches_multi_spectra_xml(
  MATCH_COLLECTION_T* match_collection,
  FILE* output){
  carp(CARP_DETAILED_DEBUG, "Writing matches to xml file");
  static int index_count = 1;
  int match_idx = 0;
  int num_matches = match_collection->match_total;
  for (match_idx = 0; match_idx < num_matches; match_idx++){
    MATCH_T* cur_match = match_collection->match[match_idx];
    BOOLEAN_T is_decoy = get_match_null_peptide(cur_match);
    Spectrum* spectrum = get_match_spectrum(cur_match);
    SpectrumZState& zstate = match_collection->zstate;
    spectrum->printXml(output, zstate, index_count);
    fprintf(output, "    <search_result>\n");
    if (! is_decoy){
      print_match_xml(cur_match, output, 
                      match_collection->scored_type );
      fprintf(output, "    </search_result>\n");
      fprintf(output, "    </spectrum_query>\n");
    }
    index_count++;
  }
  
}

/**
 * \brief Print the psm features to file in xml format
 *
 * Prints a spectrum_query tag which encompasses the search_hit tag
 * which represents peptide to spectra match.
 *
 * returns TRUE, if succesfully printed xml format of PSMs, else FALSE
 *
 */

BOOLEAN_T print_match_collection_xml(
  FILE* output,
  int top_match,
  MATCH_COLLECTION_T* match_collection,
  Spectrum* spectrum,
  SCORER_TYPE_T main_score,
  int index
  )
{
  if ( output == NULL || match_collection == NULL || spectrum == NULL ){
    return FALSE;
  }
  SpectrumZState zstate = match_collection->zstate; 
  int num_matches = match_collection->experiment_size;

  // calculate delta_cn and populate fields in the matches
  calculate_delta_cn(match_collection, SEARCH_COMMAND);

  /* print spectrum query */
  spectrum->printXml(output, zstate, index);


  MATCH_T* match = NULL;
  // create match iterator
  // TRUE: return match in sorted order of main_score type
  MATCH_ITERATOR_T* match_iterator = 
    new_match_iterator(match_collection, main_score, TRUE);
  int count = 0;
  int last_rank = 0;
  
  fprintf(output, "    <search_result>\n");
   // iterate over matches
  while(match_iterator_has_next(match_iterator)){
    match = match_iterator_next(match_iterator);    
    int cur_rank = get_match_rank(match, main_score);
  

    
    // print if we haven't reached the limit
    // or if we are at the limit but this match is a tie with the last
    if( count < top_match || last_rank == cur_rank ){
      
      print_match_xml(match, 
                      output, 
                      match_collection->scored_type);
      count++;
      last_rank = cur_rank;
    } else if( count >= top_match && last_rank != cur_rank ) {
      break;
    } // else see if there is one more tie to print

  }// next match
  
  carp(CARP_DETAILED_DEBUG, "printed %d out of %d xml matches", 
       count, num_matches);

  free_match_iterator(match_iterator);
  fprintf(output, "    </search_result>\n");
  fprintf(output, "    </spectrum_query>\n");
  
  return TRUE;

    
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
 *\returns TRUE, if sucessfully print sqt format of the PSMs, else FALSE 
 */
BOOLEAN_T print_match_collection_sqt(
  FILE* output,                  ///< the output file -out
  int top_match,                 ///< the top matches to output -in
  MATCH_COLLECTION_T* match_collection,
  ///< the match_collection to print sqt -in
  Spectrum* spectrum           ///< the spectrum to print sqt -in
  )
{

  if( output == NULL || match_collection == NULL || spectrum == NULL ){
    return FALSE;
  }
  time_t hold_time;
  hold_time = time(0);
  SpectrumZState& zstate = match_collection->zstate; 
  int num_matches = match_collection->experiment_size;

  // calculate delta_cn and populate fields in the matches
  calculate_delta_cn(match_collection, SEQUEST_COMMAND);

  // First, print spectrum info
  spectrum->printSqt(output, num_matches, zstate);
  
  MATCH_T* match = NULL;
  
  // create match iterator
  // TRUE: return match in sorted order of main_score type
  MATCH_ITERATOR_T* match_iterator = 
    new_match_iterator(match_collection, XCORR, TRUE);
  
  // Second, iterate over matches, prints M and L lines
  while(match_iterator_has_next(match_iterator)){
    match = match_iterator_next(match_iterator);    

    // print only up to max_rank_result of the matches
    if( get_match_rank(match, XCORR) > top_match ){
      break;
    }// else

    print_match_sqt(match, output);

  }// next match
  
  // print the match with Sp rank==1 if its xcorr rank > top_match rank.  
  if( get_match_rank(match_collection->top_scoring_sp, XCORR) > top_match ){
    print_match_sqt(match_collection->top_scoring_sp, output);
  }

  free_match_iterator(match_iterator);
  
  return TRUE;
}

/**
 * \brief Print the psm features to file in tab delimited format.
 *
 * Matches will be sorted by main_score and the ranks of those scores
 * will be used to determine how many matches are printed for each
 * spectrum.
 * \returns TRUE, if sucessfully print tab-delimited format of the
 * PSMs, else FALSE  
 */
BOOLEAN_T print_match_collection_tab_delimited(
  MatchFileWriter* output,                  ///< the output file -out
  int top_match,                 ///< the top matches to output -in
  MATCH_COLLECTION_T* match_collection,
  ///< the match_collection to print sqt -in
  Spectrum* spectrum,          ///< the spectrum to print sqt -in
  SCORER_TYPE_T main_score       ///< the main score to report -in
  )
{

  if( output == NULL || match_collection == NULL || spectrum == NULL ){
    return FALSE;
  }
  //int charge = match_collection->charge; 
  int num_matches = match_collection->experiment_size;
  int scan_num = spectrum->getFirstScan();
  //FLOAT_T spectrum_neutral_mass = spectrum->getNeutralMass(charge);
  FLOAT_T spectrum_precursor_mz = spectrum->getPrecursorMz();

  // calculate delta_cn and populate fields in the matches
  calculate_delta_cn(match_collection, SEARCH_COMMAND);

  MATCH_T* match = NULL;
  
  // create match iterator
  // TRUE: return match in sorted order of main_score type
  MATCH_ITERATOR_T* match_iterator = 
    new_match_iterator(match_collection, main_score, TRUE);
  int count = 0;
  int last_rank = 0;

  // iterate over matches
  while(match_iterator_has_next(match_iterator)){
    match = match_iterator_next(match_iterator);    
    int cur_rank = get_match_rank(match, main_score);

    // print if we haven't reached the limit
    // or if we are at the limit but this match is a tie with the last
    if( count < top_match || last_rank == cur_rank ){

      print_match_tab(match_collection, match, output, scan_num, 
                      spectrum_precursor_mz, 
                      num_matches);
      count++;
      last_rank = cur_rank;
    } else if( count >= top_match && last_rank != cur_rank ) {
      break;
    } // else see if there is one more tie to print

  }// next match
  
  carp(CARP_DETAILED_DEBUG, "printed %d out of %d tab matches", 
       count, num_matches);

  free_match_iterator(match_iterator);
  
  return TRUE;
}

/**
 * Retrieve the calibration parameter eta.
 */
FLOAT_T get_calibration_eta (
  MATCH_COLLECTION_T* my_collection ///< The collection -in
  )
{
  return(my_collection->eta);
}

/**
 * Retrieve the calibration parameter beta.
 */
FLOAT_T get_calibration_beta (
  MATCH_COLLECTION_T* my_collection ///< The collection -in
  )
{
  return(my_collection->beta);
}

/**
 * Retrieve the calibration parameter shift.
 */
FLOAT_T get_calibration_shift (
  MATCH_COLLECTION_T* my_collection ///< The collection -in
  )
{
  return(my_collection->shift);
}

/**
 * Retrieve the calibration correlation.
 */
FLOAT_T get_calibration_corr (
  MATCH_COLLECTION_T* my_collection ///< The collection -in
  )
{
  return(my_collection->correlation);
}

/**
 * match_iterator routines!
 *
 */

/**
 * create a new memory allocated match iterator, which iterates over
 * match iterator only one iterator is allowed to be instantiated per
 * match collection at a time 
 *\returns a new memory allocated match iterator
 */
MATCH_ITERATOR_T* new_match_iterator(
  MATCH_COLLECTION_T* match_collection,
  ///< the match collection to iterate -out
  SCORER_TYPE_T score_type,
  ///< the score type to iterate (LOGP_EXP_SP, XCORR) -in
  BOOLEAN_T sort_match  ///< should I return the match in sorted order?
  )
{
  // TODO (BF 06-Feb-08): Could we pass back an iterator with has_next==False
  if (match_collection == NULL){
    carp(CARP_FATAL, "Null match collection passed to match iterator");
  }
  // is there an existing iterator?
  if(match_collection->iterator_lock){
    carp(CARP_FATAL, 
         "Can only have one match iterator instantiated at a time");
  }
  
  // has the score type been populated in match collection?
  if(!match_collection->scored_type[score_type]){
    char score_str[64];
    scorer_type_to_string(score_type, score_str);
    carp(CARP_ERROR, "New match iterator for score type %s.", score_str);
    carp(CARP_FATAL, 
         "The match collection has not been scored for request score type.");
  }
  
  // allocate a new match iterator
  MATCH_ITERATOR_T* match_iterator = 
    (MATCH_ITERATOR_T*)mycalloc(1, sizeof(MATCH_ITERATOR_T));
  
  // set items
  match_iterator->match_collection = match_collection;
  match_iterator->match_mode = score_type;
  match_iterator->match_idx = 0;
  match_iterator->match_total = match_collection->match_total;

  // only sort if requested and match collection is not already sorted
  if(sort_match){
    sort_match_collection(match_collection, score_type);
  }


  // ok lock up match collection
  match_collection->iterator_lock = TRUE;
  
  return match_iterator;
}

/**
 * \brief Create a match iterator to return matches from a collection
 * grouped by spectrum and sorted by given score type.
 *
 * \returns A heap-allocated match iterator.
 */
MATCH_ITERATOR_T* new_match_iterator_spectrum_sorted(
  MATCH_COLLECTION_T* match_collection,  ///< for iteration -in
  SCORER_TYPE_T scorer ///< the score type to sort by -in
){

  MATCH_ITERATOR_T* match_iterator = 
    (MATCH_ITERATOR_T*)mycalloc(1, sizeof(MATCH_ITERATOR_T));

  // set up fields
  match_iterator->match_collection = match_collection;
  match_iterator->match_mode = scorer;
  match_iterator->match_idx = 0;
  match_iterator->match_total = match_collection->match_total;

  spectrum_sort_match_collection(match_collection, scorer);

  match_collection->iterator_lock = TRUE;

  return match_iterator;
}

/**
 * Does the match_iterator have another match struct to return?
 *\returns TRUE, if match iter has a next match, else False
 */
BOOLEAN_T match_iterator_has_next(
  MATCH_ITERATOR_T* match_iterator ///< the working  match iterator -in
  )
{
  return (match_iterator->match_idx < match_iterator->match_total);
}

/**
 * return the next match struct!
 *\returns the match in decreasing score order for the match_mode(SCORER_TYPE_T)
 */
MATCH_T* match_iterator_next(
  MATCH_ITERATOR_T* match_iterator ///< the working match iterator -in
  )
{
  return match_iterator->match_collection->match[match_iterator->match_idx++];
}

/**
 * free the memory allocated iterator
 */
void free_match_iterator(
  MATCH_ITERATOR_T* match_iterator ///< the match iterator to free
  )
{
  // free iterator
  if (match_iterator != NULL){
    if (match_iterator->match_collection != NULL){
      match_iterator->match_collection->iterator_lock = FALSE;
    }
    free(match_iterator);
  }
}


/**
 * \brief Print the given match collection for several spectra to
 * tab-delimited files only.  Takes the spectrum information from the
 * matches in the collection.  At least for now, prints all matches in
 * the collection rather than limiting by top-match parameter.  Uses
 * SP as preliminary score and XCORR as main score.
 */
void print_matches_multi_spectra
(MATCH_COLLECTION_T* match_collection, 
 MatchFileWriter* tab_file, 
 MatchFileWriter* decoy_tab_file){

  carp(CARP_DETAILED_DEBUG, "Writing matches to file");

  // if file location is target (i.e. tdc=T), print all to target
  MatchFileWriter* decoy_file = decoy_tab_file;
  if( get_boolean_parameter("tdc") == TRUE ){
    decoy_file = tab_file;
  }

  // for each match, get spectrum info, determine if decoy, print
  int match_idx = 0;
  int num_matches = match_collection->match_total;
  for(match_idx = 0; match_idx < num_matches; match_idx++){
    MATCH_T* cur_match = match_collection->match[match_idx];
    BOOLEAN_T is_decoy = get_match_null_peptide(cur_match);
    Spectrum* spectrum = get_match_spectrum(cur_match);
    int scan_num = spectrum->getFirstScan();
    FLOAT_T mz = spectrum->getPrecursorMz();
    FLOAT_T num_psm_per_spec = get_match_ln_experiment_size(cur_match);
    num_psm_per_spec = expf(num_psm_per_spec) + 0.5; // round to nearest int

    if( is_decoy ){
      print_match_tab(match_collection, cur_match, decoy_file, scan_num, mz, 
                      (int)num_psm_per_spec);
    }
    else{
      print_match_tab(match_collection, cur_match, tab_file, scan_num, mz,
                      (int)num_psm_per_spec);
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
MATCH_COLLECTION_T* new_match_collection_psm_output(
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator, 
    ///< the working match_collection_iterator -in
  SET_TYPE_T set_type  
    ///< what set of match collection are we creating? (TARGET, DECOY1~3) -in 
  )
{ 
  // prepare the match_collection
  MATCH_COLLECTION_T* match_collection = allocate_match_collection();

  // set this as a post_process match collection
  match_collection->post_process_collection = TRUE;
  
  // the protein counter size, create protein counter
  DATABASE_T* database = match_collection_iterator->database;
  match_collection->post_protein_counter_size 
   = get_database_num_proteins(database);
  match_collection->post_protein_counter 
   = (int*)mycalloc(match_collection->post_protein_counter_size, sizeof(int));
  match_collection->post_protein_peptide_counter 
   = (int*)mycalloc(match_collection->post_protein_counter_size, sizeof(int));

  // create hash table for peptides
  // Set initial capacity to protein count.
  match_collection->post_hash 
    = new_hash(match_collection->post_protein_counter_size);

  // get the list of files to open
  vector<string> file_names;
  get_target_decoy_filenames(file_names, 
                             match_collection_iterator->working_directory,
                             set_type);

  // open each file and add psms to match collection
  for(int file_idx = 0; file_idx < (int)file_names.size(); file_idx++){
    char* full_filename = 
      get_full_filename(match_collection_iterator->directory_name,
                        file_names[file_idx].c_str());
    MatchFileReader delimited_result_file(full_filename);
    carp(CARP_DEBUG, "Creating new match collection from '%s' file.",
         full_filename);
    free(full_filename);

    extend_match_collection_tab_delimited(match_collection, 
                                          database, 
                                          delimited_result_file);

    // for the first target file, set headers based on input files
    if( set_type == SET_TARGET && file_idx == 0 ){
      delimited_result_file.getMatchColumnsPresent(
                                  *match_collection_iterator->cols_in_file); 
    }
  } // next file

  return match_collection;
}

/**
 * parse all the match objects and add to match collection
 *\returns TRUE, if successfully parse all PSMs in result_file, else FALSE
 */
BOOLEAN_T extend_match_collection_tab_delimited(
  MATCH_COLLECTION_T* match_collection, ///< match collection to extend -out
  DATABASE_T* database, ///< the database holding the peptides -in
  MatchFileReader& result_file   ///< the result file to parse PSMs -in
  )
{


  MATCH_T* match = NULL;

  FLOAT_T delta_cn = 0;
  FLOAT_T ln_delta_cn = 0;
  FLOAT_T ln_experiment_size = 0;

  // only for post_process_collections
  if(!match_collection->post_process_collection){
    carp(CARP_ERROR, "Must be a post process match collection to extend.");
    return FALSE;
  }

  while (result_file.hasNext()) {

    /*** get spectrum specific features ***/
    match_collection-> zstate.setNeutralMass(
      result_file.getFloat(SPECTRUM_NEUTRAL_MASS_COL),
      result_file.getInteger(CHARGE_COL));
    delta_cn = result_file.getFloat(DELTA_CN_COL);
    if (delta_cn <= 0.0) {
      ln_delta_cn = 0;
    } else {
      ln_delta_cn = logf(delta_cn);
    }
    ln_experiment_size = log(result_file.getFloat(MATCHES_SPECTRUM_COL));
    

    //TODO: Parse all boolean indicators for scores
    match_collection -> 
      scored_type[SP] = !result_file.empty(SP_SCORE_COL);

    match_collection -> 
      scored_type[XCORR] = 
      !result_file.empty(XCORR_SCORE_COL);

    match_collection -> 
      scored_type[DECOY_XCORR_QVALUE] = 
      !result_file.empty(DECOY_XCORR_QVALUE_COL);

/* TODO
    match_collection -> 
      scored_type[LOGP_WEIBULL_XCORR] = 
      result_file.getString("logp weibull xcorr") != "";
*/
    match_collection -> 
      scored_type[LOGP_BONF_WEIBULL_XCORR] = 
      !result_file.empty(PVALUE_COL);

    match_collection -> 
      scored_type[PERCOLATOR_QVALUE] = 
      !result_file.empty(PERCOLATOR_QVALUE_COL);

    match_collection -> 
      scored_type[PERCOLATOR_SCORE] = 
      !result_file.empty(PERCOLATOR_SCORE_COL);

    match_collection -> 
      scored_type[LOGP_QVALUE_WEIBULL_XCORR] = 
      !result_file.empty(WEIBULL_QVALUE_COL);
  
    match_collection -> 
      scored_type[QRANKER_SCORE] = 
      !result_file.empty(QRANKER_SCORE_COL);
    
    match_collection -> 
      scored_type[QRANKER_QVALUE] = 
      !result_file.empty(QRANKER_QVALUE_COL);

    match_collection -> post_scored_type_set = TRUE;

    // parse match object
    match = parse_match_tab_delimited(result_file, database);
    if (match == NULL) {
      carp(CARP_ERROR, "Failed to parse tab-delimited PSM match");
      return FALSE;
    }

    //set all spectrum specific features to parsed match
    set_match_zstate(match, match_collection->zstate);
    set_match_delta_cn(match, delta_cn);
    set_match_ln_delta_cn(match, ln_delta_cn);
    set_match_ln_experiment_size(match, ln_experiment_size);

    //add match to match collection.
    add_match_to_post_match_collection(match_collection, match);
    //increment pointer.
    result_file.next();
  }

  return TRUE;
}



/**
 * \brief Adds the match to match_collection by copying the pointer.
 * 
 * No new match is allocated.  Match_collection total_matches must not
 * exceed the _MAX_NUMBER_PEPTIDES. 
 * \returns TRUE if successfully adds the match to the
 * match_collection, else FALSE 
 */
BOOLEAN_T add_match_to_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to free -out
  MATCH_T* match ///< the match to add -in
  )
{
  if( match_collection == NULL || match == NULL ){
    carp(CARP_FATAL, "Cannot add NULL match to NULL collection.");
  }

  // check if enough space for peptide match
  if(match_collection->match_total >= _MAX_NUMBER_PEPTIDES){
    carp(CARP_FATAL, "Cannot add to match collection; count exceeds limit: %d", 
         _MAX_NUMBER_PEPTIDES);
  }

  // add a new match to array
  match_collection->match[match_collection->match_total] = match;
  increment_match_pointer_count(match);
  
  // increment total rich match count
  ++match_collection->match_total;

  
  return TRUE;
}

/**
 * Adds the match object to match_collection
 * Must not exceed the _MAX_NUMBER_PEPTIDES to be match added
 * Only for post_process (i.e. post search) match_collections.  Keeps
 * track of all peptides in a hash table.
 * \returns TRUE if successfully adds the match to the
 * match_collection, else FALSE 
 */
// this method renamed so that a more general add_match_to_match_collection could be implemented
BOOLEAN_T add_match_to_post_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to free -out
  MATCH_T* match ///< the match to add -in
  )
{
  if( match_collection == NULL || match == NULL ){
    carp(CARP_FATAL, "Cannot add NULL match to NULL collection.");
  }
  PEPTIDE_T* peptide = NULL;

  // only for post_process_collections
  if(!match_collection->post_process_collection){
    carp(CARP_ERROR, "Must be a post process match collection to add a match.");
    return FALSE;
  }

  // check if enough space for peptide match
  if(match_collection->match_total >= _MAX_NUMBER_PEPTIDES){
    carp(CARP_ERROR, "Rich match count exceeds max match limit: %d", 
         _MAX_NUMBER_PEPTIDES);
    return FALSE;
  }

  // add a new match to array
  match_collection->match[match_collection->match_total] = match;
  increment_match_pointer_count(match);
  
  // increment total rich match count
  ++match_collection->match_total;
  
  // DEBUG, print total peptided scored so far
  if(match_collection->match_total % 1000 == 0){
    carp(CARP_INFO, "parsed PSM: %d", match_collection->match_total);
  }

  // match peptide
  peptide = get_match_peptide(match);
  
  // update protein counter, protein_peptide counter
  update_protein_counters(match_collection, peptide);
  
  // update hash table
  char* hash_value = get_peptide_hash_value(peptide); 
  add_hash(match_collection->post_hash, hash_value, NULL); 
  free(hash_value);
  
  return TRUE;
}

/**
 * updates the protein_counter and protein_peptide_counter for 
 * run specific features
 */
void update_protein_counters(
  MATCH_COLLECTION_T* match_collection, ///< working match collection -in
  PEPTIDE_T* peptide  ///< peptide information to update counters -in
  )
{
  PEPTIDE_SRC_ITERATOR_T* src_iterator = NULL;
  PEPTIDE_SRC_T* peptide_src = NULL;
  Protein* protein = NULL;
  unsigned int protein_idx = 0;
  int hash_count = 0;
  BOOLEAN_T unique = FALSE;
  
  // only for post_process_collections
  if(!match_collection->post_process_collection){
    carp(CARP_FATAL, 
         "Must be a post process match collection to update protein counter.");
  }
  
  // See if this peptide has been observed before?
  char* hash_value = get_peptide_hash_value(peptide);
  hash_count = get_hash_count(match_collection->post_hash, hash_value);
  free(hash_value);

  if(hash_count < 1){
    // yes this peptide is first time observed
    unique = TRUE;
  }
  // first update protein counter
  src_iterator = new_peptide_src_iterator(peptide);
  
  // iterate overall parent proteins
  while(peptide_src_iterator_has_next(src_iterator)){
    peptide_src = peptide_src_iterator_next(src_iterator);
    protein = get_peptide_src_parent_protein(peptide_src);
    protein_idx = protein->getProteinIdx();
    
    // update the number of PSM this protein matches
    ++match_collection->post_protein_counter[protein_idx];
    
    // number of peptides match this protein
    if(unique){
      ++match_collection->post_protein_peptide_counter[protein_idx];
    }
  }  
  free_peptide_src_iterator(src_iterator);
}

/**
 * Fill the match objects score with the given the array, and populate
 * the corresponding ranks.  The match object order must not have been
 * altered since scoring.  The result array size must equal the number
 * of matches in the given match collection.  After the function
 * completes, the match collection is sorted by the specified score,
 * unless preserve_order is set to TRUE.
 */
void fill_result_to_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< collection to iterate over -in/out
  double* results,           ///< array of scores -in
  SCORER_TYPE_T score_type,  ///< The score type of the results to fill -in
  BOOLEAN_T preserve_order   ///< preserve match order?
  )
{
  MATCH_T** match_array = NULL;
  SCORER_TYPE_T score_type_old = match_collection->last_sorted;

  // iterate over match object in collection, set scores
  int match_idx = 0;
  for(; match_idx < match_collection->match_total; ++match_idx){
    MATCH_T* match = match_collection->match[match_idx];
    set_match_score(match, score_type, results[match_idx]);    
  }
  
  // if need to preserve order store a copy of array in original order 
  if(preserve_order){
    match_array = (MATCH_T**)mycalloc(match_collection->match_total, sizeof(MATCH_T*));
    for(match_idx=0; match_idx < match_collection->match_total; ++match_idx){
      match_array[match_idx] = match_collection->match[match_idx];
    }
  }

  // populate the rank of match_collection
  if(!populate_match_rank_match_collection(match_collection, score_type)){
    free_match_collection(match_collection);
    carp(CARP_FATAL, "failed to populate match rank in match_collection");
  }
  
  // restore match order.
  if(preserve_order){
    for(match_idx=0; match_idx < match_collection->match_total; ++match_idx){
      match_collection->match[match_idx] = match_array[match_idx];
    }
    match_collection->last_sorted = score_type_old;
    free(match_array);
  }

  match_collection->scored_type[score_type] = TRUE;
}

/**
 * Process run specific features from all the PSMs
 */
void process_run_specific_features(
  MATCH_COLLECTION_T* match_collection ///< the match collection to free -out
  );

/**
 * \brief Calculate the delta_cn of each match and populate the field.
 * 
 * Delta_cn is the normalized difference between xcorrs of different
 * ranks.  For SEQUEST style searching
 * match[i] = (match[0] - match[i]) / match[0] 
 * For other searching
 * match[i] = (match[0] - match[i+1]) / match[0].  This function
 * defaults to the second case. Sorts match_collection by xcorr, if necessary.
 * 
 */
BOOLEAN_T calculate_delta_cn( MATCH_COLLECTION_T* match_collection,
                              COMMAND_T search_type ){

  if( match_collection == NULL ){
    carp(CARP_ERROR, "Cannot calculate deltaCn for NULL match collection");
    return FALSE;
  }

  if( match_collection->scored_type[XCORR] == FALSE ){
    carp(CARP_WARNING, 
      "Delta_cn not calculated because match collection not scored for xcorr");
    return FALSE;
  }

  // sort, if not already
  // N.B. Can't use sort_match_collection because iterator already exists!
  MATCH_T** matches = match_collection->match;
  int num_matches = match_collection->match_total;
  if( match_collection->last_sorted != XCORR ){
    qsort_match(matches, num_matches, (QSORT_COMPARE_METHOD)compare_match_xcorr);
    match_collection->last_sorted = XCORR;
  }

  // get xcorr of first match
  FLOAT_T max_xcorr = get_match_score(matches[0], XCORR);

  // for each match, calculate deltacn
  for(int match_idx = 0; match_idx < num_matches; match_idx++){
    FLOAT_T next_xcorr = 0;
    
    if( search_type == SEQUEST_COMMAND ){ // use this match's xcorr
      next_xcorr = get_match_score(matches[match_idx], XCORR);
    } else {                              // find next non-equal xcorr
      FLOAT_T this_xcorr = get_match_score(matches[match_idx], XCORR);
      int score_idx = match_idx + 1;
      
      while( score_idx < num_matches &&
             get_match_score(matches[score_idx], XCORR) == this_xcorr ){
        score_idx++;
      }
      
      if( score_idx < num_matches ){
        next_xcorr = get_match_score(matches[score_idx], XCORR);
      } else { // if this is the last match, set dcn to 0
        next_xcorr = max_xcorr;
      }
    }
    
    FLOAT_T delta_cn = (max_xcorr - next_xcorr) / max_xcorr;
    set_match_delta_cn(matches[match_idx], delta_cn);
  }

  return TRUE;
}


/**********************************
 * match_collection get, set methods
 **********************************/

/**
 * \returns TRUE if the match_collection only contains decoy matches,
 * else (all target or mixed) returns FALSE.
 */
BOOLEAN_T get_match_collection_is_decoy(
  MATCH_COLLECTION_T* match_collection
){
  return match_collection->null_peptide_collection;
}

/**
 *\returns the match_collection protein counter for the protein idx
 */
int get_match_collection_protein_counter(
  MATCH_COLLECTION_T* match_collection, ///< the working match collection -in
  unsigned int protein_idx ///< the protein index to return protein counter -in
  )
{
  // only for post_process_collections
  if(!match_collection->post_process_collection){
    carp(CARP_FATAL, "Must be a post process match collection to get protein counter.");
  }

  // number of PSMs match this protein
  return match_collection->post_protein_counter[protein_idx];
}

/**
 *\returns the match_collection protein peptide counter for the protein idx
 */
int get_match_collection_protein_peptide_counter(
  MATCH_COLLECTION_T* match_collection, ///< the working match collection -in
  unsigned int protein_idx ///< the protein index to return protein peptiide counter -in
  )
{
  // only for post_process_collections
  if(!match_collection->post_process_collection){
    carp(CARP_FATAL, "Must be a post process match collection to get peptide counter.");
  }
  
  // number of peptides match this protein
  return match_collection->post_protein_peptide_counter[protein_idx];
}

/**
 *\returns the match_collection hash value of PSMS for which this is the best scoring peptide
 */
int get_match_collection_hash(
  MATCH_COLLECTION_T* match_collection, ///< the working match collection -in
  PEPTIDE_T* peptide  ///< the peptide to check hash value -in
  )
{
  // only for post_process_collections
  if(!match_collection->post_process_collection){
    carp(CARP_FATAL, "Must be a post process match collection, to get match_collection_hash");
  }
  
  char* hash_value = get_peptide_hash_value(peptide);
  int count = get_hash_count(match_collection->post_hash, hash_value);
  free(hash_value);
  
  return count;
}

/**
 * \brief Get the number of proteins in the database associated with
 * this match collection.
 */
int get_match_collection_num_proteins(
  MATCH_COLLECTION_T* match_collection ///< the match collection of interest -
  ){

  return match_collection->post_protein_counter_size;
}


/******************************
 * match_collection_iterator
 ******************************/
     
/**
 * \brief Finds the next match_collection in directory and prepares
 * the iterator to hand it off when 'next' called.
 *
 * When no more match_collections (i.e. psm files) are available, set
 * match_collection_iterator->is_another_collection to FALSE
 * \returns void
 */
void setup_match_collection_iterator(
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator 
    ///< the match_collection_iterator to set up -in/out
  )
{
  // are there any more match_collections to return?
  if(match_collection_iterator->collection_idx 
      < match_collection_iterator->number_collections){
    // then go parse the match_collection
    match_collection_iterator->match_collection = 
      new_match_collection_psm_output(match_collection_iterator, 
         (SET_TYPE_T)match_collection_iterator->collection_idx);

    // we have another match_collection to return
    match_collection_iterator->is_another_collection = TRUE;
    
    // let's move on to the next one next time
    ++match_collection_iterator->collection_idx;

    // reset directory
    rewinddir(match_collection_iterator->working_directory);
  }
  else{
    // we're done, no more match_collections to return
    match_collection_iterator->is_another_collection = FALSE;
  }
}

/**
 * Create a match_collection iterator from a directory of serialized files.
 * Only handles up to one target and three decoy sets per folder.
 *\returns match_collection iterator instantiated from a result folder
 */
MATCH_COLLECTION_ITERATOR_T* new_match_collection_iterator(
  const char* output_file_directory, 
    ///< the directory path where the PSM output files are located -in
  const char* fasta_file, 
    ///< The name of the fasta file for peptides for match_collections. -in
  int* decoy_count
  )
{
  carp(CARP_DEBUG, 
       "Creating match collection iterator for dir %s and protein database %s",
       output_file_directory, fasta_file);

  // allocate match_collection
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator =
    (MATCH_COLLECTION_ITERATOR_T*)
      mycalloc(1, sizeof(MATCH_COLLECTION_ITERATOR_T));

  DIR* working_directory = NULL;
  struct dirent* directory_entry = NULL;
  DATABASE_T* database = NULL;
  BOOLEAN_T use_index = is_directory(fasta_file);

  /*
    BF: I think that this step is to count how many decoys there are
    per target file.  This is prone to errors as all it really does is
    check for the presence of a file with *decoy_1*, and one with
    *decoy_2* and *decoy_3*.  In fact, the three files could be from
    different targets.  Nothing was being done with the check for a
    target file.  There must be a better way to do this.
   */


  // do we have these files in the directory
  BOOLEAN_T boolean_result = FALSE;
  BOOLEAN_T decoy_1 = FALSE;
  BOOLEAN_T decoy_2 = FALSE;
  BOOLEAN_T decoy_3 = FALSE;

  // open PSM file directory
  working_directory = opendir(output_file_directory);
  
  if (working_directory == NULL) {
    carp(CARP_FATAL, "Failed to open PSM file directory: %s", 
        output_file_directory);
  }
  
  // determine how many decoy sets we have
  while((directory_entry = readdir(working_directory))){
    
    if(suffix_compare(directory_entry->d_name, "decoy-1.txt")) {
      carp(CARP_DEBUG, "Found decoy file %s", directory_entry->d_name);
      decoy_1 = TRUE;
    }
    else if(suffix_compare(directory_entry->d_name, "decoy.txt")) {
      decoy_1 = TRUE;
    }
    else if(suffix_compare(directory_entry->d_name, "decoy-2.txt")) {
      decoy_2 = TRUE;
    }
    else if(suffix_compare(directory_entry->d_name, "decoy-3.txt")) {
      decoy_3 = TRUE;
    }    
    else if(suffix_compare(directory_entry->d_name, ".txt")){
      carp(CARP_DEBUG, "Found target file %s", directory_entry->d_name);
      boolean_result = TRUE;
    }
    if (boolean_result && decoy_1 && decoy_2 && decoy_3) {
      break; // We've found all the files we can use.
    }
  }
  
  // set total_sets count
  int total_sets = 0;

  if(decoy_3){
    total_sets = 4; // 3 decoys + 1 target
    *decoy_count = 3;
  }
  else if(decoy_2){
    total_sets = 3; // 2 decoys + 1 target
    *decoy_count = 2;
  }
  else if(decoy_1){
    total_sets = 2; // 1 decoys + 1 target
    *decoy_count = 1;
  }
  else{
    total_sets = 1;
    *decoy_count = 0;
    carp(CARP_INFO, "No decoy sets exist in directory: %s", 
        output_file_directory);
  }
  if(!boolean_result){
    carp(CARP_FATAL, "No PSM files found in directory '%s'", 
         output_file_directory);
  }

  // get binary fasta file name with path to crux directory 
  char* binary_fasta  = NULL;
  if (use_index == TRUE){ 
    binary_fasta = get_index_binary_fasta_name(fasta_file);
  } else {
    binary_fasta = get_binary_fasta_name(fasta_file);
    carp(CARP_DEBUG, "Looking for binary fasta %s", binary_fasta);
    if (access(binary_fasta, F_OK)){
      carp(CARP_DEBUG, "Could not find binary fasta %s", binary_fasta);
      if (!create_binary_fasta_here(fasta_file, binary_fasta)){
       carp(CARP_FATAL, "Could not create binary fasta file %s", binary_fasta);
      };
    }
  }
  
  // check if input file exist
  if(access(binary_fasta, F_OK)){
    free(binary_fasta);
    carp(CARP_FATAL, "The file \"%s\" does not exist (or is not readable, "
        "or is empty) for crux index.", binary_fasta);
  }
  
  carp(CARP_DEBUG, "Creating a new database");
  // now create a database, 
  // using fasta file either binary_file(index) or fastafile
  database = new_database(binary_fasta, TRUE);
  
  // check if already parsed
  if(!get_database_is_parsed(database)){
    carp(CARP_DETAILED_DEBUG,"Parsing database");
    if(!parse_database(database)){
      carp(CARP_FATAL, "Failed to parse database, cannot create new index");
    }
  }
  
  free(binary_fasta);

  // reset directory
  rewinddir(working_directory);
  
  // set match_collection_iterator fields
  match_collection_iterator->working_directory = working_directory;
  match_collection_iterator->database = database;  
  match_collection_iterator->number_collections = total_sets;
  match_collection_iterator->directory_name = 
    my_copy_string(output_file_directory);
  match_collection_iterator->is_another_collection = FALSE;

  match_collection_iterator->cols_in_file = new vector<bool>();

  // setup the match collection iterator for iteration
  // here it will go parse files to construct match collections
  setup_match_collection_iterator(match_collection_iterator);

  if( match_collection_iterator == NULL ){
    carp(CARP_FATAL, "Failed to create a match collection iterator");
  }
  carp(CARP_DETAILED_DEBUG, "Created the match collection iterator");

  return match_collection_iterator;
}

/**
 *\returns TRUE, if there's another match_collection to return, else return FALSE
 */
BOOLEAN_T match_collection_iterator_has_next(
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator ///< the working match_collection_iterator -in
  )
{
  // Do we have another match_collection to return
  return match_collection_iterator->is_another_collection;
}

/**
 * free match_collection_iterator
 */
void free_match_collection_iterator(
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator ///< the working match_collection_iterator -in
  )
{
  // free unclaimed match_collection
  if(match_collection_iterator->match_collection != NULL){
    free_match_collection(match_collection_iterator->match_collection);
  }
  
  // if no index, remove the temp binary fasta file
  char* fasta_file = get_string_parameter("protein database");
  if( is_directory(fasta_file) == FALSE ){
    char* binary_fasta = get_binary_fasta_name(fasta_file);
    carp(CARP_DEBUG, "Protein source %s is not an index.  "
         "Removing temp binary fasta %s", fasta_file, binary_fasta);
    remove(binary_fasta);
  }
  free(fasta_file);

  // free up all match_collection_iteratory.
  free(match_collection_iterator->directory_name);
  free_database(match_collection_iterator->database);
  closedir(match_collection_iterator->working_directory); 
  delete match_collection_iterator->cols_in_file;

  free(match_collection_iterator);
}

/**
 * \brief Fetches the next match collection object and prepares for
 * the next iteration 
 *\returns The next match collection object
 */
MATCH_COLLECTION_T* match_collection_iterator_next(
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator 
    ///< the working match_collection_iterator -in
  )
{
  MATCH_COLLECTION_T* match_collection = NULL;
  
  if(match_collection_iterator->is_another_collection){
    match_collection = match_collection_iterator->match_collection;
    match_collection_iterator->match_collection = NULL;
    setup_match_collection_iterator(match_collection_iterator);
    return match_collection;
  }
  else{
    carp(CARP_ERROR, "No match_collection to return");
    return NULL;
  }
}

/**
 *\returns the total number of match_collections to return
 */
int get_match_collection_iterator_number_collections(
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator ///< the working match_collection_iterator -in
  )
{
  return match_collection_iterator->number_collections;
}

/**
 * \brief Get the name of the directory the match_collection_iterator
 * is working in.
 * \returns A heap allocated string (char*) of the directory name.
 */
char* get_match_collection_iterator_directory_name(
  MATCH_COLLECTION_ITERATOR_T* iterator ///< the match_collection_iterator -in
  ){

  char* dir_name = my_copy_string(iterator->directory_name);

  return dir_name;
}

/**
 * Try setting the match collection's charge.  Successful if the
 * current charge is 0 (i.e. hasn't yet been set) or if the current
 * charge is the same as the given value.  Otherwise, returns false
 *
 * \returns TRUE if the match_collection's charge state was changed.
 */

BOOLEAN_T set_match_collection_zstate(
  MATCH_COLLECTION_T* match_collection, ///< match collection to change
  SpectrumZState& zstate ///< new zstate
  ) {

  if (get_match_collection_charge(match_collection) == 0) {
    match_collection->zstate = zstate;
    return TRUE;
  } else {
    //error
    carp(CARP_WARNING, "Cannot change the zstate of a match collection "
        "once it has been set.");
    return FALSE;
  }
}
/*
BOOLEAN_T set_match_collection_charge(
  MATCH_COLLECTION_T* match_collection,  ///< match collection to change
  int charge){///< new charge value

  if( match_collection->charge == 0 ){
    match_collection->charge = charge;
    return TRUE;
  }// else error

  carp(CARP_WARNING, "Cannot change the charge state of a match collection "
       "once it has been set.");
  return FALSE;
}
*/

/**
 * Search the given database or index using shuffled peptides and the
 * spectrum/charge in the target psm match collection.  Add those
 * scores to the target psm match collection for use in weibull
 * parameter estimation but do not save the matches.
 */
void add_decoy_scores_match_collection(
  MATCH_COLLECTION_T* target_matches, ///< add scores to this collection
  Spectrum* spectrum, ///< search this spectrum
  int charge, ///< search spectrum at this charge state
  MODIFIED_PEPTIDES_ITERATOR_T* peptides ///< use these peptides to search
){

  // reuse these for scoring all matches
  IonConstraint* ion_constraint = 
    IonConstraint::newIonConstraintSmart(XCORR, charge);
 
  IonSeries* ion_series = new IonSeries(ion_constraint, charge);  
  SCORER_T* scorer = new_scorer(XCORR);
  
  // for each peptide in the iterator
  while( modified_peptides_iterator_has_next(peptides)){

    // get peptide and sequence
    PEPTIDE_T* peptide = modified_peptides_iterator_next(peptides);
    char* decoy_sequence = get_peptide_sequence(peptide);
    MODIFIED_AA_T* modified_seq = get_peptide_modified_aa_sequence(peptide);

    // create the ion series for this peptide
    ion_series->update(decoy_sequence, modified_seq);
    ion_series->predictIons();

    // get the score
    FLOAT_T score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);

    // add to collection's list of xcorrs
    target_matches->xcorrs[target_matches->num_xcorrs] = score;
    target_matches->num_xcorrs++;

    // clean up
    free(decoy_sequence);
    free(modified_seq);
    free_peptide(peptide);
  } // next peptide

  IonConstraint::free(ion_constraint);
  delete ion_series;
  free_scorer(scorer);

}

// cheater functions for testing

void force_scored_by(MATCH_COLLECTION_T* match_collection, SCORER_TYPE_T type){
  match_collection->scored_type[type] = TRUE;
}


/**
 * Extract a given type of score into an array.  The array is
 * allocated here and must be freed by the caller.
 */
FLOAT_T* extract_scores_match_collection(
  SCORER_TYPE_T       score_type, ///< Type of score to extract.
  MATCH_COLLECTION_T* all_matches ///< add scores to this collection
)
{
  FLOAT_T* return_value = (FLOAT_T*)mycalloc(all_matches->match_total,
                                             sizeof(FLOAT_T));

  MATCH_ITERATOR_T* match_iterator = 
    new_match_iterator(all_matches, XCORR, FALSE);
  int idx = 0;
  while(match_iterator_has_next(match_iterator)){
    MATCH_T* match = match_iterator_next(match_iterator); 
    return_value[idx] = get_match_score(match, score_type);
    idx++;
  }
  free_match_iterator(match_iterator);

  return(return_value);
}

/**
 * Given a hash table that maps from a score to its q-value, assign
 * q-values to all of the matches in a given collection.
 */
void assign_match_collection_qvalues(
  const map<FLOAT_T, FLOAT_T>* score_to_qvalue_hash,
  SCORER_TYPE_T score_type,
  MATCH_COLLECTION_T* all_matches
){

  // Iterate over the matches filling in the q-values
  MATCH_ITERATOR_T* match_iterator = new_match_iterator(all_matches, 
                                                        score_type, FALSE);
  while(match_iterator_has_next(match_iterator)){
    MATCH_T* match = match_iterator_next(match_iterator);
    FLOAT_T score = get_match_score(match, score_type);

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
    // Should never reach this point.
    case SP: 
    case LOGP_WEIBULL_XCORR: 
    case DECOY_XCORR_PEPTIDE_QVALUE:
    case LOGP_PEPTIDE_QVALUE_WEIBULL:
    case PERCOLATOR_PEPTIDE_QVALUE:
    case QRANKER_PEPTIDE_QVALUE:
    case NUMBER_SCORER_TYPES:
    case INVALID_SCORER_TYPE:
      carp(CARP_FATAL, "Something is terribly wrong!");
    }

    set_match_score(match, derived_score_type, qvalue);
    all_matches->scored_type[derived_score_type] = TRUE;

  }
  free_match_iterator(match_iterator);
}

const vector<bool>& get_match_collection_iterator_cols_in_file(
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator){

  return *match_collection_iterator->cols_in_file;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


