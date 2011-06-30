/*************************************************************************//**
 * \file q-value.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: Jan 03 2007
 * \brief  Given as input a directory containing binary psm files,
 * a protein database, and an optional parameter file, analyze the
 * matches (with percolator or q-value) and return scores indicating
 * how good the matches are. 
 *
 * Handles at most 4 files (target and decoy).  Expects psm files to
 * start with <fileroot>.se and 
 * end with the extension '.txt' and decoys to end with
 * '-decoy#.txt'.  Multiple target files in the given directory are
 * concatinated together and presumed to be non-overlaping parts of
 * the same ms2 file. 
 ****************************************************************************/
#include "q-value.h"
#include <map>

using namespace std;

static const int MAX_PSMS = 10000000;
// 14th decimal place
static const double EPSILON = 0.00000000000001;

/**
* Find the best-scoring match for each peptide in a given collection.
* Only consider the top-ranked PSM per spectrum.
*
* Results are stored in the given match collection.
*/
static void identify_best_psm_per_peptide
(MATCH_COLLECTION_T* all_matches,
 SCORER_TYPE_T score_type)
{
  /* Instantiate a hash table.  key = peptide; value = maximal xcorr
     for that peptide. */
  map<string, FLOAT_T> best_score_per_peptide;

  // Store in the hash the best score per peptide.
  MATCH_ITERATOR_T* match_iterator 
    = new_match_iterator(all_matches, score_type, FALSE);
  while(match_iterator_has_next(match_iterator)){
    Match* match = match_iterator_next(match_iterator);

    // Skip matches that are not top-ranked.
    if (match->getRank(score_type) == 1) {
      char *peptide = match->getModSequenceStrWithSymbols();
      FLOAT_T this_score = match->getScore(score_type);

      map<string, FLOAT_T>::iterator map_position 
        = best_score_per_peptide.find(peptide);

      if (map_position == best_score_per_peptide.end()) {
        best_score_per_peptide[peptide] = this_score;
      } else {
        // FIXME: Need a generic compare operator for score_type.
        if (map_position->second < this_score) {
          best_score_per_peptide[peptide] = this_score;
        }
      }
      free(peptide);
    }
  }
  free_match_iterator(match_iterator);


  // Set the best_per_peptide Boolean in the match, based on the hash.
  match_iterator = new_match_iterator(all_matches, score_type, FALSE);
  while(match_iterator_has_next(match_iterator)){
    Match* match = match_iterator_next(match_iterator);

     // Skip matches that are not top-ranked.
    if (match->getRank(score_type) == 1) {
      char* peptide = match->getModSequenceStrWithSymbols();
      FLOAT_T this_score = match->getScore(score_type);

      map<string, FLOAT_T>::iterator map_position 
        = best_score_per_peptide.find(peptide);

      if (map_position->second == this_score) {
        match->setBestPerPeptide();
        
        // Prevent ties from causing two peptides to be best.
        best_score_per_peptide[peptide] = HUGE_VAL;
      }
      
      free(peptide);
    }
  }
  free_match_iterator(match_iterator);
}


/**
 * The q-value is defined as the minimum FDR at which a given score is
 * deemed significant.  This function takes a list of FDRs and
 * converts them into q-values.  The FDRs should be ordered from
 * lowest to highest, sorted according to the underlying score.
 */
static void convert_fdr_to_qvalue 
  (FLOAT_T* qvalues,     ///< Come in as FDRs, go out as q-values.
   int      num_values)
{
  FLOAT_T prev_fdr = qvalues[num_values - 1];
  int idx;
  for (idx=num_values - 2; idx >= 0; idx--){
    carp(CARP_DETAILED_DEBUG, "fdr[%i] = %.10f", idx, qvalues[idx]);
    FLOAT_T this_fdr = qvalues[idx];
    if (prev_fdr < this_fdr) {
      qvalues[idx] = prev_fdr;
    }
    prev_fdr = qvalues[idx];
    carp(CARP_DETAILED_DEBUG, "qvalue[%i] = %.10f", idx, qvalues[idx]);
  }
}

/**
 * Store two parallel arrays of floats in a hash table.
 *
 * The new hash table must be freed by the caller.
 */
map<FLOAT_T, FLOAT_T>* store_arrays_as_hash
  (FLOAT_T* keys, 
   FLOAT_T* values,
   int      num_values
){

  map<FLOAT_T, FLOAT_T>* return_value = new map<FLOAT_T, FLOAT_T>();

  int idx;
  for (idx=0; idx < num_values; idx++){
    carp(CARP_DETAILED_DEBUG, "%g maps to %g", keys[idx], values[idx]);
    (*return_value)[keys[idx]] = values[idx];
  }
  return(return_value);
}

/**
 * Use the Benjamini-Hochberg procedure to convert a given set of
 * p-values into q-values.  
 *
 * Assumes that the input is an array of negative log p-values.  The
 * output q-values are not log-transformed.
 *
 * This function uses the command line parameter "pi-zero".
 */
FLOAT_T* compute_qvalues_from_pvalues(
  FLOAT_T* pvalues, 
  int      num_pvals,
  FLOAT_T  pi_zero
){

  // sort the - log p-values in descending order
  sort(pvalues, pvalues + num_pvals, compareDescending());

  // convert the p-values into FDRs using Benjamini-Hochberg
  FLOAT_T* qvalues = (FLOAT_T*)mycalloc(num_pvals, sizeof(FLOAT_T));
  int idx;
  for (idx=0; idx < num_pvals; idx++){
    carp(CARP_DETAILED_DEBUG, "pvalue[%i] = %.10f", idx, exp(-pvalues[idx]));
    double fdr = (exp(-pvalues[idx]) / (idx + 1)) 
      * (FLOAT_T)num_pvals * pi_zero;
    qvalues[idx] = fdr;
    carp(CARP_DETAILED_DEBUG, "FDR[%i] = %.10f", idx, qvalues[idx]);
  }

  // convert the FDRs into q-values
  convert_fdr_to_qvalue(qvalues, num_pvals);

  return(qvalues);
}

/**
 * \brief Compute q-values from a given set of scores, using a second
 * set of scores as an empirical null.  Sorts the incoming target
 * scores and returns a corresponding list of q-values.
 *
 * This function is only exported to allow unit testing.
 */
FLOAT_T* compute_decoy_qvalues(
  FLOAT_T* target_scores,
  int      num_targets,
  FLOAT_T* decoy_scores,
  int      num_decoys,
  FLOAT_T  pi_zero
){
  if ((num_targets == 0) || (num_decoys == 0)) {
    carp(CARP_FATAL, "Cannot compute q-values (%d targets, %d nulls).",
         num_targets, num_decoys);
  }
  carp(CARP_DEBUG, "Computing decoy q-values.");

  int target_idx;
  for (target_idx = 0; target_idx < num_targets; target_idx++) {
    carp(CARP_DEBUG, "target_scores[%d]=%g decoy_scores[%d]=%g",
         target_idx, target_scores[target_idx],
         target_idx, decoy_scores[target_idx]);
  }

  // Sort both sets of scores.
  sort(target_scores, target_scores + num_targets, compareDescending());
  sort(decoy_scores, decoy_scores + num_decoys, compareDescending());

  for (target_idx = 0; target_idx < num_targets; target_idx++) {
    carp(CARP_DEBUG, "target_scores[%d]=%g decoy_scores[%d]=%g",
         target_idx, target_scores[target_idx],
         target_idx, decoy_scores[target_idx]);
  }

  // Precompute the ratio of targets to decoys.
  FLOAT_T targets_to_decoys = (FLOAT_T)num_targets / (FLOAT_T)num_decoys;

  // Compute false discovery rate for each target score.
  FLOAT_T* qvalues = (FLOAT_T*)mycalloc(num_targets, sizeof(FLOAT_T));
  int decoy_idx = 0;
  for (target_idx = 0; target_idx < num_targets; target_idx++) {
    FLOAT_T target_score = target_scores[target_idx];

    // Find the index of the first decoy score greater than this target score.
    while ((decoy_idx < num_decoys) &&
           (decoy_scores[decoy_idx] > target_score)) {
      decoy_idx++;
    }

    // FDR = #decoys / #targets
    FLOAT_T fdr = pi_zero * targets_to_decoys * 
      ((FLOAT_T)decoy_idx / (FLOAT_T)(target_idx + 1));
    carp(CARP_DEBUG, "target_idx=%d target_score=%g decoy_idx=%d fdr=%g",
         target_idx, target_score, decoy_idx, fdr);
    
    if ( fdr > 1.0 ){
      fdr = 1.0;
    }
    
    qvalues[target_idx] = fdr;
  }
  
  // Convert the FDRs into q-values.
  convert_fdr_to_qvalue(qvalues, num_targets);

  return(qvalues);
}

/**
 * \brief Compute a q-values based on what is in the PSM files in the
 * directory.  Store q-values in the match collection returned.
 *
 * If p-values were computed, then perform Benjamini-Hochberg q-value
 * calculations. Otherwise, if decoys are present, then rank on xcorr
 * and compute empirical q-values based on the number of decoys and
 * targets above the score threshold.
 *
 * \returns a collection of target PSMs with one q-value in each
 * match.
 */
MATCH_COLLECTION_T* run_qvalue(
  char* input_directory, 
  char* fasta_file,
  OutputFiles& output 
  ){

  int num_decoys = 0; // to be set by match_collection_iterator
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator =
    new_match_collection_iterator(input_directory, fasta_file, &num_decoys);
  if( num_decoys > 1 ){
    carp(CARP_FATAL, "Only one decoy file per target can be processed "
         "but %d were found.  Please move extra decoy files.", num_decoys);
  }

  // Create two match collections, for targets and decoys.
  MATCH_COLLECTION_T* target_matches = new_empty_match_collection(FALSE);
  MATCH_COLLECTION_T* decoy_matches = new_empty_match_collection(TRUE);
  set_match_collection_scored_type(target_matches, XCORR, TRUE);
  set_match_collection_scored_type(decoy_matches, XCORR, TRUE);

  // Did we find something from which to get q-values?
  BOOLEAN_T have_pvalues = FALSE;
  BOOLEAN_T have_decoys = FALSE;

  // Iterate over all match collections in this directory.
  while(match_collection_iterator_has_next(match_collection_iterator)){
    MATCH_COLLECTION_T* match_collection = 
      match_collection_iterator_next(match_collection_iterator);

    // Keep track of whether we got p-values.
    // N.B. Assumes that if one collection has p-values, they all do.
    have_pvalues = get_match_collection_scored_type(match_collection,
                                                    LOGP_BONF_WEIBULL_XCORR); 
    
    // Iterate, gathering matches into one or two collections.
    MATCH_ITERATOR_T* match_iterator =
      new_match_iterator(match_collection, XCORR, FALSE);
    while(match_iterator_has_next(match_iterator)){
      Match* match = match_iterator_next(match_iterator);
      
      // Only use top-ranked matches.
      if( match->getRank(XCORR) != 1 ){
        continue;
      }
      
      if (match->getNullPeptide() == TRUE) {
        add_match_to_match_collection(decoy_matches, match);
        have_decoys = TRUE;
      } else {
        add_match_to_match_collection(target_matches, match);
      }
      Match::freeMatch(match);
    }
    free_match_iterator(match_iterator);
    free_match_collection(match_collection);
  }

  // get from the input files which columns to print in the output files
  const vector<bool>& cols_to_print =
    get_match_collection_iterator_cols_in_file(match_collection_iterator);
  output.writeHeaders(cols_to_print);

  // Compute q-values from p-values.
  FLOAT_T* pvalues = NULL; // N.B. Misnamed for decoy calculation.
  int num_pvals = get_match_collection_match_total(target_matches);
  FLOAT_T* qvalues = NULL;
  SCORER_TYPE_T score_type = INVALID_SCORER_TYPE;
  if (have_pvalues == TRUE) {
    carp(CARP_DEBUG, "There are %d PSMs for q-value computation.", num_pvals);
    set_match_collection_scored_type(target_matches, 
                                     LOGP_BONF_WEIBULL_XCORR, 
                                     TRUE);
    pvalues = extract_scores_match_collection(LOGP_BONF_WEIBULL_XCORR,
                                              target_matches);
    qvalues = compute_qvalues_from_pvalues(pvalues, num_pvals,
                                           get_double_parameter("pi-zero"));
    score_type = LOGP_BONF_WEIBULL_XCORR;
  }

  // Compute q-values from the XCorr decoy distribution.
  else if (have_decoys == TRUE) {
    int num_decoys = get_match_collection_match_total(decoy_matches);
    carp(CARP_DEBUG,
         "There are %d target and %d decoy PSMs for q-value computation.",
         num_pvals, num_decoys);
    pvalues = extract_scores_match_collection(XCORR, target_matches);
    FLOAT_T* decoy_xcorrs 
      = extract_scores_match_collection(XCORR, decoy_matches);
    qvalues = compute_decoy_qvalues(pvalues, num_pvals, 
                                    decoy_xcorrs, num_decoys,
                                    get_double_parameter("pi-zero"));
    free(decoy_xcorrs);
    score_type = XCORR;
  }

  // Fatal: Cannot compute q-values.
  else {
    carp(CARP_FATAL, "Cannot compute q-values without decoy PSMs or p-values.");
  }

  // Store p-values to q-values as a hash, and then assign them.
  map<FLOAT_T, FLOAT_T>* qvalue_hash 
    = store_arrays_as_hash(pvalues, qvalues, num_pvals);
  assign_match_collection_qvalues(qvalue_hash, score_type, target_matches);
  free(pvalues);
  free(qvalues);
  delete qvalue_hash;

  // Identify PSMs that are top-scoring per peptide.
  identify_best_psm_per_peptide(target_matches, score_type);

  // Compute peptide-level q-values.
  //  compute_decoy_q_values(all_matches/
  // TRUE); // Do peptide-level scoring.

  // Store targets by score.
  sort_match_collection(target_matches, score_type);
  output.writeMatches(target_matches);

  free_match_collection(decoy_matches);
  free_match_collection_iterator(match_collection_iterator);

  return(target_matches);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
