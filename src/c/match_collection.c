/*****************************************************************************
 * \file match_collection.c
 * AUTHOR: Chris Park
 * CREATE DATE: 11/27 2006
 * DESCRIPTION: Object for given a database and a spectrum, generate all match objects
 * REVISION: 
 ****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>
#include "carp.h"
#include "parse_arguments.h"
#include "spectrum.h"
#include "spectrum_collection.h"
#include "ion.h"
#include "ion_series.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "scorer.h" 
#include "generate_peptides_iterator.h" 
#include "match.h"
#include "match_collection.h"

static BOOLEAN_T is_first_spectrum = TRUE;
/**
 *\struct match_collection
 *\brief An object that contains match objects with a given spectrum and peptide database
 */
struct match_collection{
  MATCH_T* match[_MAX_NUMBER_PEPTIDES]; ///< array of match object
  BOOLEAN_T scored_type[_SCORE_TYPE_NUM]; ///< has the score type been computed in each match
  int experiment_size; ///< total peptide experiment sample size(peptide count form the database before any truncation
  int match_total; ///< total_match_count
  SCORER_TYPE_T last_sorted; ///< the last type the match has been sorted(if -1, then unsorted, if ever change the order must change to -1)
  BOOLEAN_T iterator_lock; ///< is there a iterator been curretly created?, if TRUE cannot manipulate match collection
  float sp_scores_mean;  ///< the mean value of the scored peptides sp score
  float mu; ///< EVD parameter Xcorr(characteristic value of extreme value distribution)
  float l_value; ///< EVD parameter Xcorr(decay constant of extreme value distribution)
};

/**
 *\struct match_iterator
 *\brief An object that iterates over the match objects in the specified score type (SP, XCORR)
 */
struct match_iterator{
  MATCH_COLLECTION_T* match_collection; ///< the match collection to iterate -out
  SCORER_TYPE_T match_mode; ///< the current score working mode (SP, XCORR)
  int match_idx; ///< current match to return
  int match_total; ///< total_match_count
};

/**
 * typedef, for descrition look below.
 */
BOOLEAN_T score_match_collection_sp(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  SPECTRUM_T* spectrum, ///< the spectrum to match peptides -in
  int charge,       ///< the charge of the spectrum -in
  GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator ///< peptide iteartor to use, must set it first before use
  );

BOOLEAN_T score_match_collection_xcorr(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  SPECTRUM_T* spectrum, ///< the spectrum to match peptides -in
  int charge       ///< the charge of the spectrum -in
  );

/**
 * The match collection must be scored under SP first
 * \returns TRUE, if successfully scores matches for LOGP_EXP_SP
 */
BOOLEAN_T score_match_collection_logp_exp_sp(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  int peptide_to_score ///< the number of top ranked sp scored peptides to score for logp_exp_sp -in
  );

/**
 * The match collection must be scored under SP first
 * \returns TRUE, if successfully scores matches for LOGP_BONF_EXP_SP
 */
BOOLEAN_T score_match_collection_logp_bonf_exp_sp(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  int peptide_to_score ///< the number of top ranked sp scored peptides to score for logp_bonf_exp_sp -in
  );

/**
 * The match collection must be scored under Xcorr first
 * \returns TRUE, if successfully scores matches for LOGP_EVD_XCORR
 */
BOOLEAN_T score_match_collection_logp_evd_xcorr(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  int peptide_to_score ///< the number of top ranked xcorr scored peptides to score for logp_evd_xcorr -in
  );

/**
 * The match collection must be scored under Xcorr first
 * \returns TRUE, if successfully scores matches for LOGP_BONF_EVD_XCORR
 */
BOOLEAN_T score_match_collection_logp_bonf_evd_xcorr(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  int peptide_to_score ///< the number of top ranked xcorr scored peptides to score for logp_bonf_evd_xcorr -in
  );


/**
 * Randomly samples max_count peptides from the peptide distribution and try to esitimate the Xcorr distribution of the the entire peptide distribution 
 * from the sampled peptide distribution. Populates the two EVD perameters mu, lambda in the match_collection.
 *
 * This function finds the location parameter, mu, and scale parameter, 1/L, 
 * that maximize the log likelihood of the data given an extreme value 
 * distribution.  It finds the parameters by using Newton-Raphson to find 
 * the zero of the constraint function.  The zero of the constraint function 
 * corresponds to the scale parameter giving the maximum log likelihood for the 
 * data.
 *
 * The parameter values contains the list of the data values.
 * The parameter starting_L contains a staring guess for L.
 * The parameter contains the tolerence for determining convergence.
 *
 * Returns the values of mu and L that maximize the log likelihood.
 * Throws an exception if Newton-Raphson fails to converge.
 *\returns TRUE, if successfully calculates the EVD parameters for the xcorr peptide distribution., else FALSE.
 */
BOOLEAN_T estimate_evd_perameters(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to estimate evd perameters -out
  int sample_count, ///< the number of peptides to sample from the match_collection -in
  SCORER_TYPE_T score_type, ///< score_type to estimate EVD distribution -in
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  int charge       ///< the charge of the spectrum -in
  );

/**
 * keeps the top max_rank number of matches and frees the rest
 * sorts by score_type(SP, XCORR, ...)
 */
void truncate_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to truncate -out
  int max_rank,     ///< max number of top rank matches to keep from SP -in
  SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
  );

/**
 * \returns An (empty) match_collection object.
 */
MATCH_COLLECTION_T* allocate_match_collection()
{
  MATCH_COLLECTION_T* match_collection =
    (MATCH_COLLECTION_T*)mycalloc(1, sizeof(MATCH_COLLECTION_T));
    
  //loop over to set all score type to FALSE
  int score_type_idx = 0;
  for(; score_type_idx < _SCORE_TYPE_NUM; ++score_type_idx){
    match_collection->scored_type[score_type_idx] = FALSE;
  }
  
  //set last score to -1, thus nothing has been done yet
  match_collection->last_sorted = -1;
  match_collection->iterator_lock = FALSE;

  return match_collection;
}

/**
 * free the memory allocated match collection
 */
void free_match_collection(
  MATCH_COLLECTION_T* match_collection ///< the match collection to free -out
  )
{
  //free all matches, actually we are only decrementing the pointer count in each match object
  while(match_collection->match_total > 0){
    --match_collection->match_total;
    free_match(match_collection->match[match_collection->match_total]);
  }
  free(match_collection);
}

/**
 * create a new match collection from spectrum
 * creates a peptide iterator for given mass window
 * return the top max_rank matches, first scored by prelim_score(SP), then by score_type(XCORR, LOGP_EXP_SP, LOGP_BONF_EXP_SP)
 *\returns a new match_collection object that is scored by score_type and contains the top max_rank matches
 */
MATCH_COLLECTION_T* new_match_collection_spectrum(
 SPECTRUM_T* spectrum, ///< the spectrum to match peptides -in
 int charge,       ///< the charge of the spectrum -in
 int max_rank,     ///< max number of top rank matches to keep from SP -in
 SCORER_TYPE_T prelim_score, ///< the preliminary score type (SP) -in
 SCORER_TYPE_T score_type ///< the score type (XCORR, LOGP_EXP_SP, LOGP_BONF_EXP_SP) -in
 )
{
  MATCH_COLLECTION_T* match_collection = allocate_match_collection();

  //top_rank_for_p_value is the amount of top ranked sp scored peptides to score for LOGP_EXP_SP
  //This parameter can only be set from crux_parameter file
  int top_rank_for_p_value = get_int_parameter("top-rank-p-value", 1);
  int sample_count = get_int_parameter("sample-count",500);
  
  //create a generate peptide iterator
  GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator =  //FIXME use neutral_mass, might chage to pick
    new_generate_peptides_iterator_sp(get_spectrum_neutral_mass(spectrum, charge));
  
  /***************Preliminary scoring**************************/
  //score SP match_collection
  if(prelim_score == SP){
    if(!score_match_collection_sp(match_collection, spectrum, charge, peptide_iterator)){
      carp(CARP_ERROR, "failed to score match collection for SP");
    }
  }

  //if scoring for EVD xcorr, sample from the original distribution of peptides 
  //for evd parameter estimation
  //Sample before truncate match collection so that the sampling will be from 
  //the entire peptide distribution.
  if(score_type == LOGP_EVD_XCORR || score_type == LOGP_BONF_EVD_XCORR){
    estimate_evd_perameters(match_collection, sample_count, XCORR, spectrum, charge);
  }

  //save only the top max_rank matches from prelim_scoring, sort and free the other matches
  truncate_match_collection(match_collection, max_rank, prelim_score);

  /***************Main scoring*******************************/
  //should we score for LOGP_EXP_SP?
  if(score_type == LOGP_EXP_SP){
    //score the top_rank_for_p_value amount of top ranked peptides their -log(p_value) match_collection
    if(!score_match_collection_logp_exp_sp(match_collection, top_rank_for_p_value)){
      carp(CARP_ERROR, "failed to score match collection for LOGP_EXP_SP");
    }
  }
  //should we score for LOGP_BONF_EXP_SP?
  else if(score_type == LOGP_BONF_EXP_SP){
    //score the top_rank_for_p_value amount of top ranked peptides their -log(p_value * number_of_peptides_scored) match_collection
    if(!score_match_collection_logp_bonf_exp_sp(match_collection, top_rank_for_p_value)){
      carp(CARP_ERROR, "failed to score match collection for LOGP_BONF_EXP_SP");
    }
  }
  //should we score for XCORR?
  else if(score_type == XCORR || score_type == LOGP_BONF_EVD_XCORR || score_type == LOGP_EVD_XCORR){
    //for all score types score for xcorr
    if(!score_match_collection_xcorr(match_collection, spectrum, charge)){
      carp(CARP_ERROR, "failed to score match collection for XCORR");
    }
    
    //Additional scoring? (EVD p_value based scoring)
    //should we score for LOGP_BONF_EVD_XCORR?
    if(score_type == LOGP_BONF_EVD_XCORR){
      //score the top_rank_for_p_value amount of top ranked peptides their -log(p_value * number_of_peptides_scored) match_collection
      if(!score_match_collection_logp_bonf_evd_xcorr(match_collection, top_rank_for_p_value)){
        carp(CARP_ERROR, "failed to score match collection for LOGP_BONF_EVD_XCORR");
      }
    }
    //should we score for LOGP_EVD_XCORR?
    else if(score_type == LOGP_EVD_XCORR){
      //score the top_rank_for_p_value amount of top ranked peptides their -log(p_value) match_collection
      if(!score_match_collection_logp_evd_xcorr(match_collection, top_rank_for_p_value)){
        carp(CARP_ERROR, "failed to score match collection for LOGP_EVD_XCORR");
      }
    }
  }
  
  //free generate_peptides_iterator
  free_generate_peptides_iterator(peptide_iterator);
  
  return match_collection;
}

/**
 * create a new match collection from spectrum
 * return the top max_rank matches, first scored by prelim_score(SP), then by score_type(XCORR, LOGP_EXP_SP, LOGP_BONF_EXP_SP);
 * uses a provided peptide iterator, MUST be a mutable iterator
 * Sets the iterator before useage.
 *\returns a new match_collection object that is scored by score_type and contains the top max_rank matches
 */
MATCH_COLLECTION_T* new_match_collection_spectrum_with_peptide_iterator(
 SPECTRUM_T* spectrum, ///< the spectrum to match peptides -in
 int charge,       ///< the charge of the spectrum -in
 int max_rank,     ///< max number of top rank matches to keep from SP -in 
 SCORER_TYPE_T prelim_score, ///< the preliminary score type (SP) -in
 SCORER_TYPE_T score_type ///< the score type (XCORR, LOGP_EXP_SP, LOGP_BONF_EXP_SP) -in
 //GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator ///< peptide iteartor to use, must set it first before use
 )
{
  MATCH_COLLECTION_T* match_collection = allocate_match_collection();
  //get perameters
  float neutral_mass = get_spectrum_neutral_mass(spectrum, charge);
  double mass_window = get_double_parameter("mass-window", 3);
  int sample_count = get_int_parameter("sample-count",500);
  double min_mass = neutral_mass - mass_window;
  double max_mass = neutral_mass + mass_window;

  //top_rank_for_p_value is the amount of top ranked sp scored peptides to score for LOGP_EXP_SP
  //This parameter can only be set from crux_parameter file
  int top_rank_for_p_value = get_int_parameter("top-rank-p-value", 1);

  carp(CARP_DEBUG,"searching peptide in %.2f ~ %.2f", min_mass, max_mass); 
  
  //free(peptide_iterator);

  //move out of crux index dir
  if(!is_first_spectrum){
    chdir("..");
  }else{
    is_first_spectrum = FALSE;
  }

  GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator =  //NULL;//FIXME use neutral_mass, might chage to pick
    new_generate_peptides_iterator_mutable();
  //chdir("human-NCBI-042104_crux_index/");

  //set the generate_peptides_iterator for the next round of peptides
  set_generate_peptides_mutable(peptide_iterator, max_mass, min_mass);
  
  /**********Preliminary scoring*****************/
  //score SP match_collection
  if(prelim_score == SP){
    if(!score_match_collection_sp(match_collection, spectrum, charge, peptide_iterator)){
      carp(CARP_ERROR, "failed to score match collection for SP");
    }
  }

  //if scoring for EVD xcorr, sample from the original distribution of peptides 
  //for evd parameter estimation
  //Sample before truncate match collection so that the sampling will be from 
  //the entire peptide distribution.
  if(score_type == LOGP_EVD_XCORR || score_type == LOGP_BONF_EVD_XCORR){
    estimate_evd_perameters(match_collection, sample_count, XCORR, spectrum, charge);
  }

  //save only the top max_rank matches from prelim_scoring, sort and free the other matches
  truncate_match_collection(match_collection, max_rank, prelim_score);
  
  /**********Main scoring***************/
  //should we score for LOGP_EXP_SP?
  if(score_type == LOGP_EXP_SP){
    //score the top_rank_for_p_value amount of top ranked peptides their -log(p_value) match_collection
    if(!score_match_collection_logp_exp_sp(match_collection, top_rank_for_p_value)){
      carp(CARP_ERROR, "failed to score match collection for LOGP_EXP_SP");
    }
  }
  //should we score for LOGP_BONF_EXP_SP?
  else if(score_type == LOGP_BONF_EXP_SP){
    //score the top_rank_for_p_value amount of top ranked peptides their -log(p_value * number_of_peptides_scored) match_collection
    if(!score_match_collection_logp_bonf_exp_sp(match_collection, top_rank_for_p_value)){
      carp(CARP_ERROR, "failed to score match collection for LOGP_BONF_EXP_SP");
    }
  }
  //should we score for XCORR?
  else if(score_type == XCORR || score_type == LOGP_BONF_EVD_XCORR || score_type == LOGP_EVD_XCORR){
    //for all score types score for xcorr
    if(!score_match_collection_xcorr(match_collection, spectrum, charge)){
      carp(CARP_ERROR, "failed to score match collection for XCORR");
    }
    
    //Additional scoring with Xcorr? (EVD p_value based scoring)
    //should we score for LOGP_BONF_EVD_XCORR?
    if(score_type == LOGP_BONF_EVD_XCORR){
      //score the top_rank_for_p_value amount of top ranked peptides their -log(p_value * number_of_peptides_scored) match_collection
      if(!score_match_collection_logp_bonf_evd_xcorr(match_collection, top_rank_for_p_value)){
        carp(CARP_ERROR, "failed to score match collection for LOGP_BONF_EVD_XCORR");
      }
    }
    //should we score for LOGP_EVD_XCORR?
    else if(score_type == LOGP_EVD_XCORR){
      //score the top_rank_for_p_value amount of top ranked peptides their -log(p_value) match_collection
      if(!score_match_collection_logp_evd_xcorr(match_collection, top_rank_for_p_value)){
        carp(CARP_ERROR, "failed to score match collection for LOGP_EVD_XCORR");
      }
    }
  }
  
  //free generate_peptides_iterator
  free_generate_peptides_iterator(peptide_iterator);
  
  return match_collection;
}

/**
 * sort the match collection by score_type(SP, XCORR, ... )
 *\returns TRUE, if successfully sorts the match_collection
 */
BOOLEAN_T sort_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
  )
{
  //check if we are allowed to alter match_collection
  if(match_collection->iterator_lock){
    carp(CARP_ERROR, "cannot alter match_collection when a match iterator is already been instantiated");
    return FALSE;
  }

  switch(score_type){
  case SP:
    //sort the match to decreasing SP order for the return
    qsort_match(match_collection->match, match_collection->match_total, (void *)compare_match_sp);
    match_collection->last_sorted = SP;
    return TRUE;
  case XCORR:
    //sort the match to decreasing XCORR order for the return
    qsort_match(match_collection->match, match_collection->match_total, (void *)compare_match_xcorr);
    match_collection->last_sorted = XCORR;
    return TRUE;
  case LOGP_EVD_XCORR:
    //LOGP_EVD_XCORR and XCORR have same order, sort the match to decreasing XCORR order for the return
    qsort_match(match_collection->match, match_collection->match_total, (void *)compare_match_xcorr);
    match_collection->last_sorted = XCORR;
    return TRUE;
  case LOGP_BONF_EVD_XCORR:
    //LOGP_BONF_EVD_XCORR and XCORR have same order, sort the match to decreasing XCORR order for the return
    qsort_match(match_collection->match, match_collection->match_total, (void *)compare_match_xcorr);
    match_collection->last_sorted = XCORR;
    return TRUE;
  case DOTP:
    //implement later
    return FALSE;
  case LOGP_EXP_SP:
    //LOGP_EXP_SP and SP have same order, thus sort the match to decreasing SP order for the return
    qsort_match(match_collection->match, match_collection->match_total, (void *)compare_match_sp);
    match_collection->last_sorted = SP;
    return TRUE;
  case LOGP_BONF_EXP_SP:
    //LOGP_EXP_SP and SP have same order, thus sort the match to decreasing SP order for the return
    qsort_match(match_collection->match, match_collection->match_total, (void *)compare_match_sp);
    match_collection->last_sorted = SP;
    return TRUE;
  }
  
  return FALSE;
}

/**
 * keeps the top max_rank number of matches and frees the rest
 * sorts by score_type(SP, XCORR, ...)
 */
void truncate_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to truncate -out
  int max_rank,     ///< max number of top rank matches to keep from SP -in
  SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
  )
{
  //sort match collection by score type
  if(!sort_match_collection(match_collection, score_type)){
    carp(CARP_ERROR, "failed to sort match collection");
    exit(-1);
  }

  //is there any matches to free?
  while(match_collection->match_total > max_rank){
    free_match(match_collection->match[match_collection->match_total - 1]);
    --match_collection->match_total;
  }
}

/**
 * Must provide a match_collection that is already scored and ranked in the score_type
 * Rank 1, means hight score
 *\returns TRUE, if successfully popluates the match rank in the match collection
 */
BOOLEAN_T populate_match_rank_match_collection(
 MATCH_COLLECTION_T* match_collection, ///< the match collection to populate match rank -out
 SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
 )
{
  //check if the match collection is in the correct sorted order
  if(match_collection->last_sorted != score_type){
    //sort match collection by score type
    if(!sort_match_collection(match_collection, score_type)){
      carp(CARP_ERROR, "failed to sort match collection");
      return FALSE;
    }
  }
  
  //set match rank for all match objects
  int match_index = 0;
  for(; match_index < match_collection->match_total; ++match_index){
    set_match_rank(match_collection->match[match_index], score_type, match_index+1);
  }
  
  return TRUE;
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
  srand(time(NULL));

  //ranomly select matches upto count_max
  for(; count_idx < count_max; ++count_idx){
    match_idx = ((double)rand()/((double)RAND_MAX + (double)1)) * match_collection->match_total;
    
    //match_idx = rand() % match_collection->match_total;
    sample_collection->match[count_idx] = match_collection->match[match_idx];
    //increment pointer count of the match object 
    increment_match_pointer_count(sample_collection->match[count_idx]);
  }
  
  //set total number of matches sampled
  sample_collection->match_total = count_idx;

  sample_collection->experiment_size = match_collection->experiment_size;

  //set scored types in the sampled matches
  for(; score_type_idx < _SCORE_TYPE_NUM;  ++score_type_idx){
    sample_collection->scored_type[score_type_idx] = match_collection->scored_type[score_type_idx];
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
  MATCH_COLLECTION_T* match_collection, ///< the match collection to estimate evd perameters -in
  SCORER_TYPE_T score_type, ///< score_type to estimate EVD distribution -in
  float l_value,  ///< L value -in
  float* function,  ///< the output function value -out
  float* derivative,  ///< the output derivative value -out
  float* exponential_sum ///< the final exponential array sum -out
  )
{
  int idx = 0;
  float* exponential = (float*)mycalloc(match_collection->match_total, sizeof(float));
  float numerator = 0;
  float second_numerator = 0;
  float score = 0;
  float denominator = 0;
  float score_sum = 0;
  MATCH_T** matches = match_collection->match;

  //iterate over the matches to calculate numerator, exponential value, denominator
  for(; idx < match_collection->match_total; ++idx){
    score = get_match_score(matches[idx], score_type);
    exponential[idx] = exp(-l_value * score);
    numerator += (exponential[idx] * score);
    denominator += exponential[idx];
    score_sum += score;
    second_numerator += (score * score * exponential[idx]);
  }

  //assign function value
  *function = (1.0 / l_value) - (score_sum / match_collection->match_total) 
    + (numerator / denominator);

  //assign derivative value
  *derivative =  ((numerator * numerator) / (denominator * denominator)) 
    - ((second_numerator / denominator)) - (1.0 / (l_value * l_value));

  //assign the total sum of the exponential values
  *exponential_sum = denominator;

  //free exponential array
  free(exponential);
}

/**
 * Randomly samples max_count peptides from the peptide distribution and try to esitimate the Xcorr distribution of the the entire peptide distribution 
 * from the sampled peptide distribution. Populates the two EVD perameters mu, lambda in the match_collection.
 *
 * This function finds the location parameter, mu, and scale parameter, 1/L, 
 * that maximize the log likelihood of the data given an extreme value 
 * distribution.  It finds the parameters by using Newton-Raphson to find 
 * the zero of the constraint function.  The zero of the constraint function 
 * corresponds to the scale parameter giving the maximum log likelihood for the 
 * data.
 *
 * The parameter values contains the list of the data values.
 * The parameter starting_L contains a staring guess for L.
 * The parameter contains the tolerence for determining convergence.
 *
 * Returns the values of mu and L that maximize the log likelihood.
 * Throws an exception if Newton-Raphson fails to converge.
 *\returns TRUE, if successfully calculates the EVD parameters for the xcorr peptide distribution., else FALSE.
 */
BOOLEAN_T estimate_evd_perameters(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to estimate evd perameters -out
  int sample_count, ///< the number of peptides to sample from the match_collection -in
  SCORER_TYPE_T score_type, ///< score_type to estimate EVD distribution -in
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  int charge       ///< the charge of the spectrum -in
  )
{
  //randomly sample from match collection
  MATCH_COLLECTION_T* sample_collection = random_sample_match_collection(match_collection, sample_count);
  float l_value = 1;
  float f = 0.0;
  float f_prime = 0.0;
  float epsilon = 0.001;
  float exponential_sum = 0;
  int max_iterations = 10000;
  int idx = 0;

  //print info
  carp(CARP_INFO, "Estimate EVD perameters, sample count: %d", sample_count);
  
  //first score the sample match_collection
  if(score_type == XCORR){
    if(!score_match_collection_xcorr(sample_collection, spectrum, charge)){
      carp(CARP_ERROR, "failed to score match collection for XCORR");
    }
  }
  //FIXME Add different scorin gif needed
  //...

  //estimate the EVD parameters
  for(; idx < max_iterations; ++idx){
    constraint_function(sample_collection, score_type, l_value, 
                                       &f, &f_prime, &exponential_sum);

    if(fabsf(f) < epsilon){
      break;
    }
    else{
      l_value = l_value - f / f_prime;
    }
    
    //failed to converge error..
    if(idx >= max_iterations){
      carp(CARP_ERROR, "Root finding failed to converge.");
      return FALSE;
    }
  }
  
  //Calculate best value of position parameter from best value of 
  //scale parameter.
  match_collection->mu = -1.0 / l_value * logf(1.0 / sample_count * exponential_sum);
  match_collection->l_value = l_value;
    
  //free up sampled match_collection 
  free_match_collection(sample_collection);
  
  //DEBUG
  //carp(CARP_DEBUG, "mu: %.5f, L: %.5f", match_collection->mu, match_collection->l_value);
  return TRUE;
}


/**
 * scores the match_collection, the score type SP
 * Assumes this is the first time scoring with this score_collection,
 * thus, prior number match object is 0.
 * the routine will use generate_peptides for each peptide will create a match
 * that maps the peptide to the spectrum.
 * If the score has already been computed simply returns TRUE 
 *\returns  TRUE, if successfully populates the sp score matches in the match_collection
 */
BOOLEAN_T score_match_collection_sp(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  SPECTRUM_T* spectrum, ///< the spectrum to match peptides -in
  int charge,       ///< the charge of the spectrum -in
  GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator ///< the peptide iterator to score
  )
{
  
  //is this a empty collection?
  if(match_collection->match_total != 0){
    carp(CARP_ERROR, "must start with empty match collection");
    return FALSE;
  }
  //FIXME, might need this..
  //has the match collection already been scored by sp?
  /*
  if(match_collection->scored_type[SP]){
    carp(CARP_INFO, "match collection has already been scored in SP");
    return TRUE;
  }
  */

  //create a generate peptide iterator
  //GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator =  //FIXME use neutral_mass, might chage to pick
  //  new_generate_peptides_iterator_sp(get_spectrum_neutral_mass(spectrum, charge));
  
  //set ion constraint to sequest settings
  ION_CONSTRAINT_T* ion_constraint = 
    new_ion_constraint_sequest_sp(charge); 
  
  //create new scorer
  SCORER_T* scorer = new_scorer(SP);  

  char* peptide_sequence = NULL;
  MATCH_T* match = NULL;
  float score = 0;
  PEPTIDE_T* peptide = NULL;
  ION_SERIES_T* ion_series = NULL;
  
  //iterate over all peptides
  while(generate_peptides_iterator_has_next(peptide_iterator)){
    peptide = generate_peptides_iterator_next(peptide_iterator);
    peptide_sequence = get_peptide_sequence(peptide);
    

    //create new ion series
    ion_series = new_ion_series(peptide_sequence, charge, ion_constraint);
    
    //now predict ions
    predict_ions(ion_series);
    
    //calculates the Sp score
    score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);

    //increment the total sp score
    match_collection->sp_scores_mean += score;

    //create a new match
    match = new_match();
    
    //set all fields in match
    set_match_score(match, SP, score);
    set_match_peptide(match, peptide);
    set_match_spectrum(match, spectrum);
    
    //check if enough space for peptide match
    if(match_collection->match_total >= _MAX_NUMBER_PEPTIDES){
      carp(CARP_ERROR, "peptide count exceed max match limit: %d", _MAX_NUMBER_PEPTIDES);
      //free heap
      free(peptide_sequence);
      free_ion_series(ion_series);
      free_scorer(scorer);
      free_ion_constraint(ion_constraint);

      return FALSE;
    }
    
    //add a new match to array
    match_collection->match[match_collection->match_total] = match;
    
    //increment total match count
    ++match_collection->match_total;

    //DEBUG, print total peptided scored so far
    if(match_collection->match_total % 1000 == 0){
      carp(CARP_INFO, "scored peptide for sp: %d", match_collection->match_total);
    }
    
    free(peptide_sequence);
    free_ion_series(ion_series);
  }

  //calculate the final sp score mean
  match_collection->sp_scores_mean /= match_collection->match_total;
  
  //total peptide experiment sample size
  match_collection->experiment_size = match_collection->match_total;

  //DEBUG, print total peptided scored so far
  carp(CARP_INFO, "total peptide scored for sp: %d", match_collection->match_total);
  
  //free heap
  free_scorer(scorer);
  free_ion_constraint(ion_constraint);
  //free_generate_peptides_iterator(peptide_iterator);

  //DEBUG
  /*
  MATCH_COLLECTION_T* sample_collection = random_sample_match_collection(match_collection, 300);
  
  //debug
  free_match_collection(match_collection);
  *match_collection1 = sample_collection;
  match_collection = sample_collection;
  */

  //now that the match_collection is sorted, populate the rank of each match object
  if(!populate_match_rank_match_collection(match_collection, SP)){
    carp(CARP_ERROR, "failed to populate match rank for SP in match_collection");
    free_match_collection(match_collection);
    exit(-1);
  }
  
  //yes, we have now scored for the match-mode: SP
  match_collection->scored_type[SP] = TRUE;
    
  return TRUE;
}

/**
 * The match collection must be scored under SP first
 * \returns TRUE, if successfully scores matches for LOGP_EXP_SP
 */
BOOLEAN_T score_match_collection_logp_exp_sp(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  int peptide_to_score ///< the number of top ranked sp scored peptides to score for logp_exp_sp -in
  )
{
  int match_idx = 0;
  float score = 0;
  MATCH_T* match = NULL;
  
  //has the score type been populated in match collection?
  if(!match_collection->scored_type[SP]){
    carp(CARP_ERROR, "the collection must be scored by SP first before LOGP_EXP_SP");
    exit(-1);
  }

  //sort by SP if not already sorted.
  //This enables to identify the top ranked SP scoring peptides
  if(match_collection->last_sorted != SP){
    //sort match collection by score type
    if(!sort_match_collection(match_collection, SP)){
      carp(CARP_ERROR, "failed to sort match collection by SP");
      free_match_collection(match_collection);
      exit(-1);
    }
  }
  
  //we are string xcorr!
  carp(CARP_INFO, "start scoring for LOGP_EXP_SP");

  //iterate over all matches to score for LOGP_EXP_SP
  while(match_idx < match_collection->match_total && match_idx < peptide_to_score){
    match = match_collection->match[match_idx];
    score = score_logp_exp_sp(get_match_score(match, SP), match_collection->sp_scores_mean);
    
    //set all fields in match
    set_match_score(match, LOGP_EXP_SP, score);
    ++match_idx;
  }
  
  //we are done
  carp(CARP_INFO, "total peptides scored for LOGP_EXP_SP: %d", match_idx);

  //match_collection is not populate with the rank of LOGP_EXP_SP, becuase the SP rank is  identical to the LOGP_EXP_SP rank
  
  //yes, we have now scored for the match-mode: LOGP_EXP_SP
  match_collection->scored_type[LOGP_EXP_SP] = TRUE;
  
  return TRUE;
}


/**
 * The match collection must be scored under SP first
 * \returns TRUE, if successfully scores matches for LOGP_BONF_EXP_SP
 */
BOOLEAN_T score_match_collection_logp_bonf_exp_sp(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  int peptide_to_score ///< the number of top ranked sp scored peptides to score for logp_bonf_exp_sp -in
  )
{
  int match_idx = 0;
  float score = 0;
  MATCH_T* match = NULL;
  
  //has the score type been populated in match collection?
  if(!match_collection->scored_type[SP]){
    carp(CARP_ERROR, "the collection must be scored by SP first before LOGP_EXP_SP");
    exit(-1);
  }

  //sort by SP if not already sorted.
  //This enables to identify the top ranked SP scoring peptides
  if(match_collection->last_sorted != SP){
    //sort match collection by score type
    if(!sort_match_collection(match_collection, SP)){
      carp(CARP_ERROR, "failed to sort match collection by SP");
      free_match_collection(match_collection);
      exit(-1);
    }
  }
  
  //we are string xcorr!
  carp(CARP_INFO, "start scoring for LOGP_BONF_EXP_SP");

  //iterate over all matches to score for LOGP_BONF_EXP_SP
  while(match_idx < match_collection->match_total && match_idx < peptide_to_score){
    match = match_collection->match[match_idx];
    score = score_logp_bonf_exp_sp(get_match_score(match, SP), match_collection->sp_scores_mean, match_collection->experiment_size);
    
    //set all fields in match
    set_match_score(match, LOGP_BONF_EXP_SP, score);
    ++match_idx;
  }
  
  //we are done
  carp(CARP_INFO, "total peptides scored for LOGP_BONF_EXP_SP: %d", match_idx);
    
  //match_collection is not populate with the rank of LOGP_BONF_EXP_SP, becuase the SP rank is  identical to the LOGP_EXP_SP rank
  
  //yes, we have now scored for the match-mode: LOGP_BONF_EXP_SP
  match_collection->scored_type[LOGP_BONF_EXP_SP] = TRUE;
  
  return TRUE;
}

/**
 * Assumes that match collection was scored under SP first
 * \returns TRUE, if successfully scores matches for xcorr
 */
BOOLEAN_T score_match_collection_xcorr(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  SPECTRUM_T* spectrum, ///< the spectrum to match peptides -in
  int charge       ///< the charge of the spectrum -in
  )
{
  int match_idx = 0;
  MATCH_T* match = NULL;
  char* peptide_sequence = NULL;
  ION_SERIES_T* ion_series = NULL;
  float score = 0;
  
  /*
  //is this a empty collection?
  if(match_collection->match_total > 0){
    carp(CARP_ERROR, "must start with SP scored match collection");
    return FALSE;
  }
  */
  
  //set ion constraint to sequest settings
  ION_CONSTRAINT_T* ion_constraint = 
    new_ion_constraint_sequest_xcorr(charge); 
  
  //create new scorer
  SCORER_T* scorer = new_scorer(XCORR);  

  //we are string xcorr!
  carp(CARP_INFO, "start scoring for XCORR");

  //iterate over all matches to score for xcorr
  for(; match_idx < match_collection->match_total; ++match_idx){
    match = match_collection->match[match_idx];
    peptide_sequence = get_peptide_sequence(get_match_peptide(match));
    
    //create new ion series
    ion_series = new_ion_series(peptide_sequence, charge, ion_constraint);
    
    //now predict ions
    predict_ions(ion_series);
    
    //calculates the Xcorr score
    score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);
    

    //set all fields in match
    set_match_score(match, XCORR, score);
    
    //free heap
    free(peptide_sequence);
    free_ion_series(ion_series);
  }

  //we are starting xcorr!
  carp(CARP_INFO, "total peptides scored for XCORR: %d", match_idx);

  //free heap
  free_scorer(scorer);
  free_ion_constraint(ion_constraint);

  //sort match collection by score type
  if(!sort_match_collection(match_collection, XCORR)){
    carp(CARP_ERROR, "failed to sort match collection");
    exit(-1);
  }
  
  //now that the match_collection is sorted, populate the rank of each match object
  if(!populate_match_rank_match_collection(match_collection, XCORR)){
    carp(CARP_ERROR, "failed to populate match rank for Xcorr in match_collection");
    free_match_collection(match_collection);
    exit(-1);
  }

  //yes, we have now scored for the match-mode: XCORR
  match_collection->scored_type[XCORR] = TRUE;

  return TRUE;
}

/**
 * The match collection must be scored under Xcorr first
 * \returns TRUE, if successfully scores matches for LOGP_EVD_XCORR
 */
BOOLEAN_T score_match_collection_logp_evd_xcorr(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  int peptide_to_score ///< the number of top ranked xcorr scored peptides to score for logp_evd_xcorr -in
  )
{
  int match_idx = 0;
  float score = 0;
  MATCH_T* match = NULL;
  
  //has the score type been populated in match collection?
  if(!match_collection->scored_type[XCORR]){
    carp(CARP_ERROR, "the collection must be scored by XCORR first before LOGP_EVD_XCORR");
    exit(-1);
  }

  //sort by XCORR if not already sorted.
  //This enables to identify the top ranked XCORR scoring peptides
  if(match_collection->last_sorted != XCORR){
    //sort match collection by score type
    if(!sort_match_collection(match_collection, XCORR)){
      carp(CARP_ERROR, "failed to sort match collection by XCORR");
      free_match_collection(match_collection);
      exit(-1);
    }
  }
  
  //we are starting LOGP_EVD_XCORR!
  carp(CARP_INFO, "start scoring for LOGP_EVD_XCORR");

  //iterate over all matches to score for LOGP_EVD_XCORR
  while(match_idx < match_collection->match_total && match_idx < peptide_to_score){
    match = match_collection->match[match_idx];
    score = score_logp_evd_xcorr(get_match_score(match, XCORR), match_collection->mu, match_collection->l_value);
    
    //set all fields in match
    set_match_score(match, LOGP_EVD_XCORR, score);
    ++match_idx;
  }
  
  //we are done
  carp(CARP_INFO, "total peptides scored for LOGP_EVD_XCORR: %d", match_idx);

  //match_collection is not populate with the rank of LOGP_EVD_XCORR, 
  //becuase the XCORR rank is  identical to the LOGP_EVD_XCORR rank
  
  //yes, we have now scored for the match-mode: LOGP_EVD_XCORR
  match_collection->scored_type[LOGP_EVD_XCORR] = TRUE;
  
  return TRUE;
}

/**
 * The match collection must be scored under Xcorr first
 * \returns TRUE, if successfully scores matches for LOGP_BONF_EVD_XCORR
 */
BOOLEAN_T score_match_collection_logp_bonf_evd_xcorr(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  int peptide_to_score ///< the number of top ranked xcorr scored peptides to score for logp_bonf_evd_xcorr -in
  )
{
  int match_idx = 0;
  float score = 0;
  MATCH_T* match = NULL;
  
  //has the score type been populated in match collection?
  if(!match_collection->scored_type[XCORR]){
    carp(CARP_ERROR, "the collection must be scored by XCORR first before LOGP_BONF_EVD_XCORR");
    exit(-1);
  }

  //sort by XCORR if not already sorted.
  //This enables to identify the top ranked XCORR scoring peptides
  if(match_collection->last_sorted != XCORR){
    //sort match collection by score type
    if(!sort_match_collection(match_collection, XCORR)){
      carp(CARP_ERROR, "failed to sort match collection by XCORR");
      free_match_collection(match_collection);
      exit(-1);
    }
  }
  
  //we are starting LOGP_BONF_EVD_XCORR!
  carp(CARP_INFO, "start scoring for LOGP_BONF_EVD_XCORR");

  //iterate over all matches to score for LOGP_BONF_EVD_XCORR
  while(match_idx < match_collection->match_total && match_idx < peptide_to_score){
    match = match_collection->match[match_idx];
    score = score_logp_bonf_evd_xcorr(get_match_score(match, XCORR), match_collection->mu, match_collection->l_value, match_collection->experiment_size);
    
    //set all fields in match
    set_match_score(match, LOGP_BONF_EVD_XCORR, score);
    ++match_idx;
  }
  
  //we are done
  carp(CARP_INFO, "total peptides scored for LOGP_BONF_EVD_XCORR: %d", match_idx);

  //match_collection is not populate with the rank of LOGP_BONF_EVD_XCORR, 
  //becuase the XCORR rank is  identical to the LOGP_BONF_EVD_XCORR rank
  
  //yes, we have now scored for the match-mode: LOGP_BONF_EVD_XCORR
  match_collection->scored_type[LOGP_BONF_EVD_XCORR] = TRUE;
  
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
 *\returns TRUE, if there is a  match_iterators instantiated by match collection 
 */
BOOLEAN_T get_match_collection_iterator_lock(
  MATCH_COLLECTION_T* match_collection ///< the match collection to iterate -in
  )
{
  return match_collection->iterator_lock;
}

/**
 *\returns the total match objects in match_collection
 */
int get_match_collection_match_total(
  MATCH_COLLECTION_T* match_collection ///< the match collection to iterate -in
  )
{
  return match_collection->match_total;
}


/**
 * match_iterator routines!
 *
 */

/**
 * create a new memory allocated match iterator, which iterates over match iterator
 * only one iterator is allowed to be instantiated per match collection at a time
 *\returns a new memory allocated match iterator
 */
MATCH_ITERATOR_T* new_match_iterator(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to iterate -out
  SCORER_TYPE_T score_type, ///< the score type to iterate (LOGP_EXP_SP, XCORR) -in
  BOOLEAN_T sort_match  ///< should I return the match in sorted order?
  )
{
  //is there any existing iterators?
  if(match_collection->iterator_lock){
    carp(CARP_ERROR, "can only have one match iterator instantiated at a time");
    exit(-1);
  }
  
  //has the score type been populated in match collection?
  if(!match_collection->scored_type[score_type]){
    carp(CARP_ERROR, "the collection has not been score for request score type");
    exit(-1);
  }
  
  //allocate a new match iterator
  MATCH_ITERATOR_T* match_iterator = (MATCH_ITERATOR_T*)mycalloc(1, sizeof(MATCH_ITERATOR_T));
  
  //set items
  match_iterator->match_collection = match_collection;
  match_iterator->match_mode = score_type;
  match_iterator->match_idx = 0;
  match_iterator->match_total = match_collection->match_total;

  //only sort if requested and match collection is not already sorted
  if(sort_match && (match_collection->last_sorted != score_type /*|| (match_collection->last_sorted == SP && score_type == LOGP_EXP_SP)*/)){

    if((score_type == LOGP_EXP_SP || score_type == LOGP_BONF_EXP_SP) &&
       match_collection->last_sorted == SP){
      //No need to sort, since the score_type has same rank as SP      
    }
    
    else if((score_type == LOGP_EVD_XCORR || score_type == LOGP_BONF_EVD_XCORR) &&
       match_collection->last_sorted == XCORR){
      //No need to sort, since the score_type has same rank as XCORR
    }
    //sort match collection by score type
    else if(!sort_match_collection(match_collection, score_type)){
      carp(CARP_ERROR, "failed to sort match collection");
      free_match_collection(match_collection);
      free(match_iterator);
      exit(-1);
    }
  }

  //ok lock up match collection
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
  //iterator lock now set to FALSE
  match_iterator->match_collection->iterator_lock = FALSE;

  //free iterator
  free(match_iterator);
}

/**********************************
 * match_collection get, set methods
 **********************************/



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

