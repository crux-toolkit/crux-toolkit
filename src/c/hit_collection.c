/*********************************************************************//**
 * FILE: \file hit_collection.c
 * AUTHOR: Aaron Klammer
 * DESCRIPTION: \brief A collection of hits.
 * CREATE DATE: 2008 March 11
 * REVISION: $Revision: 1.1 $
 ****************************************************************************/
#include "hit_collection.h"

//static BOOLEAN_T is_first_spectrum = TRUE;

/* Private data types (structs) */

/**
 * \struct hit_collection
 * \brief An object that contains a set of hit objects.
 *
 * May contain hites for one spectrum or many spectra. 
 * 
 * 
 */
struct hit_collection{
  HIT_T* hit[_MAX_NUMBER_HITS]; ///< array of hit object
  BOOLEAN_T scored_type[_SCORE_TYPE_NUM]; 
    ///< has the score type been computed in each hit
  int experiment_size; 
    ///< total peptide count from the database before any truncation
  int hit_total; ///< total_hit_count
  SCORER_TYPE_T last_sorted; 
    ///< the last type by which it's been sorted ( -1 if unsorted)
  BOOLEAN_T iterator_lock; 
    ///< has an itterator been created? if TRUE can't manipulate hites
  int charge; ///< charge of the associated spectrum
  BOOLEAN_T null_peptide_collection; ///< are the searched peptides null

  // values used for various scoring functions.
  float delta_cn; ///< the difference in top and second Xcorr scores
  float sp_scores_mean;  ///< the mean value of the scored peptides sp score
  float mu; 
  ///< EVD parameter Xcorr(characteristic value of extreme value distribution)
  float l_value; 
  ///< EVD parameter Xcorr(decay constant of extreme value distribution)
  int top_fit_sp; 
  ///< The top ranked sp scored peptides to use as EXP_SP parameter estimation
  float base_score_sp; 
 ///< The lowest sp score within top_fit_sp, used as the base to rescale sp
  float eta;  ///< The eta parameter for the Weibull distribution.i
  float beta; ///< The beta parameter for the Weibull distribution.
  float shift; ///< The location parameter for the Weibull distribution.

  // The following features (post_*) are only valid when
  // post_process_collection boolean is TRUE 
  BOOLEAN_T post_process_collection; 
  ///< Is this a post process hit_collection?
  int post_protein_counter_size; 
  ///< the size of the protein counter array, usually the number of proteins in database
  int* post_protein_counter; 
  ///< the counter for how many each protein has hites other PSMs
  int* post_protein_peptide_counter; 
  ///< the counter for how many each unique peptides each protein has hites other PSMs
  HASH_T* post_hash; ///< hash table that keeps tracks of the peptides
  BOOLEAN_T post_scored_type_set; 
  ///< has the scored type been confirmed for the hit collection,
  // set after the first hit collection is extended
};

/**
 *\struct hit_iterator
 *\brief An object that iterates over the hit objects in the
 * specified hit_collection for the specified score type (SP, XCORR)
 */
struct hit_iterator{
  MATCH_COLLECTION_T* hit_collection; 
                            ///< the hit collection to iterate -out
  SCORER_TYPE_T hit_mode; ///< the current working score (SP, XCORR)
  int hit_idx;            ///< current hit to return
  int hit_total;          ///< total_hit_count
};

/**
 * \struct hit_collection_iterator
 * \brief An object that iterates over the hit_collection objects in
 * the specified directory of serialized hit_collections 
 */
struct hit_collection_iterator{
  DIR* working_directory; 
  ///< the working directory for the iterator to find hit_collections
  char* directory_name; ///< the directory name in char
  DATABASE_T* database; ///< the database for which the hit_collection
  int number_collections; 
  ///< the total number of hit_collections in the directory (target+decoy)
  int collection_idx;  ///< the index of the current collection to return
  MATCH_COLLECTION_T* hit_collection; ///< the hit collection to return
  BOOLEAN_T is_another_collection; 
  ///< is there another hit_collection to return?
};

/******* Private function declarations, described in definintions below ***/

BOOLEAN_T score_hit_collection_sp(
  MATCH_COLLECTION_T* hit_collection, 
  SPECTRUM_T* spectrum, 
  int charge,
  GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator
  );

BOOLEAN_T score_hit_collection_xcorr(
  MATCH_COLLECTION_T* hit_collection,
  SPECTRUM_T* spectrum,
  int charge
  );

BOOLEAN_T score_hit_collection_logp_exp_sp(
  MATCH_COLLECTION_T* hit_collection, 
  int peptide_to_score 
  );

BOOLEAN_T score_hit_collection_logp_bonf_exp_sp(
  MATCH_COLLECTION_T* hit_collection, 
  int peptide_to_score 
  );

BOOLEAN_T score_hit_collection_logp_weibull_sp(
  MATCH_COLLECTION_T* hit_collection, 
  int peptide_to_score 
  );

BOOLEAN_T score_hit_collection_logp_bonf_weibull_sp(
  MATCH_COLLECTION_T* hit_collection, 
  int peptide_to_score 
  );

BOOLEAN_T score_hit_collection_logp_weibull_xcorr(
  MATCH_COLLECTION_T* hit_collection, 
  int peptide_to_score 
  );

BOOLEAN_T score_hit_collection_logp_bonf_weibull_xcorr(
  MATCH_COLLECTION_T* hit_collection, 
  int peptide_to_score 
  );

BOOLEAN_T score_hit_collection_logp_evd_xcorr(
  MATCH_COLLECTION_T* hit_collection, 
  int peptide_to_score 
  );

BOOLEAN_T score_hit_collection_logp_bonf_evd_xcorr(
  MATCH_COLLECTION_T* hit_collection, 
  int peptide_to_score 
  );

BOOLEAN_T estimate_evd_parameters(
  MATCH_COLLECTION_T* hit_collection, 
  int sample_count, 
  SCORER_TYPE_T score_type, 
  SPECTRUM_T* spectrum,    
  int charge       
  );

BOOLEAN_T estimate_exp_sp_parameters(
  MATCH_COLLECTION_T* hit_collection, 
  int top_count 
  );

BOOLEAN_T estimate_weibull_parameters(
  MATCH_COLLECTION_T* hit_collection, 
  SCORER_TYPE_T score_type,
  int sample_count, 
  SPECTRUM_T* spectrum,
  int charge
  );

void truncate_hit_collection(
  MATCH_COLLECTION_T* hit_collection, 
  int max_rank,     
  SCORER_TYPE_T score_type 
  );

BOOLEAN_T extend_hit_collection(
  MATCH_COLLECTION_T* hit_collection, 
  DATABASE_T* database, 
  FILE* result_file   
  );

BOOLEAN_T add_hit_to_hit_collection(
  MATCH_COLLECTION_T* hit_collection, 
  HIT_T* hit 
  );

void update_protein_counters(
  MATCH_COLLECTION_T* hit_collection, 
  PEPTIDE_T* peptide  
  );

/********* end of function declarations *******************/


/**
 * \returns An (empty) hit_collection object.
 */
MATCH_COLLECTION_T* allocate_hit_collection()
{
  MATCH_COLLECTION_T* hit_collection =
    (MATCH_COLLECTION_T*)mycalloc(1, sizeof(MATCH_COLLECTION_T));
    
  // loop over to set all score type to FALSE
  int score_type_idx = 0;
  for(; score_type_idx < _SCORE_TYPE_NUM; ++score_type_idx){
    hit_collection->scored_type[score_type_idx] = FALSE;
  }
  
  // set last score to -1, thus nothing has been done yet
  hit_collection->last_sorted = -1;
  hit_collection->iterator_lock = FALSE;
  hit_collection->post_process_collection = FALSE;
  hit_collection->null_peptide_collection = FALSE;
  carp(CARP_DETAILED_DEBUG, "Allocate hit coll"); // MEMLEAK
  
  return hit_collection;
}

/**
 * /brief Free the memory allocated for a hit collection
 */
void free_hit_collection(
  MATCH_COLLECTION_T* hit_collection ///< the hit collection to free -out
  )
{
  carp(CARP_DETAILED_DEBUG, "Free hit coll"); // MEMLEAK
  // decrement the pointer count in each hit object
  // MEMLEAK
  while(hit_collection->hit_total > 0){
    --hit_collection->hit_total;
    free_hit(hit_collection->hit[hit_collection->hit_total]);
    hit_collection->hit[hit_collection->hit_total] = NULL;
  }
  
  // free post_process_collection specific memory
  if(hit_collection->post_process_collection){
    // free protein counter
    free(hit_collection->post_protein_counter);
    
    // free protein peptide counter
    free(hit_collection->post_protein_peptide_counter);
  
    // free hash table
    free_hash(hit_collection->post_hash);
  }

  free(hit_collection);
}

// TODO (BF 1-28-08): max_rank, scores, offset can be taken from parameter.c

/**
 * \brief Creates a new hit collection by searching a database
 * for hites to a spectrum. in .c
 *
 * \details This is the main spectrum searching routine.  Allocates memory for
 * the hit collection. Creates a peptide iterator for given mass
 * window. Performs preliminary scoring on all candidate
 * peptides. Performs primary scoring on the <max_rank> best-scoring
 * peptides. Estimates EVD parameters. in .c
 *
 * \returns A new hit_collection object that is scored by score_type
 * and contains the top max_rank hites in .c
 * \callgraph
 */
MATCH_COLLECTION_T* new_hit_collection_from_spectrum(
 SPECTRUM_T* spectrum, 
    ///< the spectrum to hit peptides in -in
 int charge,       
   ///< the charge of the spectrum -in
 int max_rank,     
   ///< max number of top rank hites to keep from SP -in
 SCORER_TYPE_T prelim_score, 
   ///< the preliminary score type (SP) -in
 SCORER_TYPE_T score_type, 
   ///< the score type (XCORR, LOGP_EXP_SP, LOGP_BONF_EXP_SP) -in
 float mass_offset,  
   ///< the mass offset from neutral_mass to search for candidate peptides -in
 BOOLEAN_T null_peptide_collection,
   ///< is this hit_collection a null peptide collection? -in
 INDEX_T* index,      ///< the index source of peptides
 DATABASE_T* database ///< the database (fasta) soruce of peptides
 )
{
  MATCH_COLLECTION_T* hit_collection = allocate_hit_collection();
  
  // set charge of hit_collection creation
  hit_collection->charge = charge;
  hit_collection->null_peptide_collection = null_peptide_collection;

  int top_rank_for_p_value = get_int_parameter("top-rank-p-value");
  int sample_count = get_int_parameter("sample-count");
  int top_fit_sp = get_int_parameter("top-fit-sp");
  
  // move out of crux index dir
  /*  if(!is_first_spectrum){
    chdir("..");
  }else{
    is_first_spectrum = FALSE;
    }*/
  
  // create a generate peptide iterator
  // FIXME use neutral_mass for now, but should allow option to change
  GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator =  
    new_generate_peptides_iterator_from_mass(
        get_spectrum_neutral_mass(spectrum, charge) + mass_offset,
        index, database
        );
  
  /***************Preliminary scoring**************************/
  // When creating hit objects for first time, must set the
  // null peptide boolean parameter
  
  // score SP hit_collection
  if(prelim_score == SP){
    if(!score_hit_collection_sp(
          hit_collection, 
          spectrum, 
          charge, 
          peptide_iterator)){
      carp(CARP_ERROR, "Failed to score hit collection for SP");
      free_hit_collection(hit_collection);
      return NULL;
    }
    if (hit_collection->hit_total == 0){
      carp(CARP_WARNING, "No hites found for spectrum %i charge %i",
          get_spectrum_first_scan(spectrum), charge);
      free_hit_collection(hit_collection);
      return NULL;
    }
  }

  /*********** Estimate parameters *****************************/
  //if fitting an EVD xcorr, sample from the original distribution of peptides
  // for evd parameter estimation
  // Sample before truncate hit collection so that the sampling will be
  // from the entire peptide distribution.
  if(score_type == LOGP_EVD_XCORR || score_type == LOGP_BONF_EVD_XCORR){
    estimate_evd_parameters(
        hit_collection, 
        sample_count, 
        XCORR, 
        spectrum, 
        charge);
  }
  // if scoring for LOGP_EXP_SP, LOGP_BONF_EXP_SP estimate parameters
  else if(score_type == LOGP_EXP_SP || score_type == LOGP_BONF_EXP_SP){
    estimate_exp_sp_parameters(hit_collection, top_fit_sp);
  }
  // if scoring for LOGP_EXP_SP, LOGP_BONF_EXP_SP estimate parameters
  else if(score_type == LOGP_WEIBULL_SP || 
          score_type == LOGP_BONF_WEIBULL_SP){
    estimate_weibull_parameters(
        hit_collection, SP, sample_count, spectrum, charge);
  }
  else if(score_type == LOGP_WEIBULL_XCORR || 
          score_type == LOGP_BONF_WEIBULL_XCORR){
    estimate_weibull_parameters(
        hit_collection, XCORR, sample_count, spectrum, charge);
  }

  // save only the top max_rank hites from prelim_scoring
  truncate_hit_collection(hit_collection, max_rank, prelim_score);
  
  /***************Main scoring*******************************/
  //replace with switch(score_type){ case():success = score...;} if !success 
  if(score_type == LOGP_EXP_SP){
    if(!score_hit_collection_logp_exp_sp(
          hit_collection, top_rank_for_p_value)){
      carp(CARP_ERROR, "failed to score hit collection for LOGP_EXP_SP");
    }
  }
  else if(score_type == LOGP_BONF_EXP_SP){
    if(!score_hit_collection_logp_bonf_exp_sp(
          hit_collection, top_rank_for_p_value)){
      carp(CARP_ERROR, "failed to score hit collection for LOGP_BONF_EXP_SP");
    }
  }
  else if(score_type == LOGP_WEIBULL_SP){
    carp(CARP_DEBUG, "Scoring hit collection for LOGP_WEIBULL_SP");
    if(!score_hit_collection_logp_weibull_sp(
          hit_collection, top_rank_for_p_value)){
      carp(CARP_ERROR, "failed to score hit collection for LOGP_WEIBULL_SP");
    }
  }
  else if(score_type == LOGP_BONF_WEIBULL_SP){
    if(!score_hit_collection_logp_bonf_weibull_sp(
          hit_collection, top_rank_for_p_value)){
      carp(CARP_ERROR, 
          "failed to score hit collection for LOGP_BONF_WEIBULL_SP");
    }
  }
  else if(score_type == XCORR || 
          score_type == LOGP_BONF_EVD_XCORR || 
          score_type == LOGP_EVD_XCORR || 
          score_type == LOGP_BONF_WEIBULL_XCORR || 
          score_type == LOGP_WEIBULL_XCORR ){
    if(!score_hit_collection_xcorr(hit_collection, spectrum, charge)){
      carp(CARP_ERROR, 
      "Failed to score hit collection for XCORR for spectrum %d, charge %d",
           get_spectrum_first_scan(spectrum), charge);
    }
    
    if(score_type == LOGP_BONF_EVD_XCORR){
      if(!score_hit_collection_logp_bonf_evd_xcorr(hit_collection, top_rank_for_p_value)){
        carp(CARP_ERROR, "failed to score hit collection for LOGP_BONF_EVD_XCORR");
      }
    }
    else if(score_type == LOGP_EVD_XCORR){
      if(!score_hit_collection_logp_evd_xcorr(hit_collection, top_rank_for_p_value)){
        carp(CARP_ERROR, "failed to score hit collection for LOGP_EVD_XCORR");
      }
    }
    else if(score_type == LOGP_WEIBULL_XCORR){
      if(!score_hit_collection_logp_weibull_xcorr(hit_collection, top_rank_for_p_value)){
        carp(CARP_ERROR, "failed to score hit collection for LOGP_WEIBULL_XCORR");
      }
    }
    else if(score_type == LOGP_BONF_WEIBULL_XCORR){
      if(!score_hit_collection_logp_bonf_weibull_xcorr(hit_collection, top_rank_for_p_value)){
        carp(CARP_ERROR, "failed to score hit collection for LOGP_BONF_WEIBULL_XCORR");
      }
    }
  }
  
  // free generate_peptides_iterator
  free_generate_peptides_iterator(peptide_iterator);

  carp(CARP_DETAILED_DEBUG, 
       "Finished creating hit collection for spectrum %d, charge %d",
       get_spectrum_first_scan(spectrum), charge);  
  return hit_collection;
}

/**
 * sort the hit collection by score_type(SP, XCORR, ... )
 *\returns TRUE, if successfully sorts the hit_collection
 */
BOOLEAN_T sort_hit_collection(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to sort -out
  SCORER_TYPE_T score_type ///< the score type (SP, XCORR) to sort by -in
  )
{
  // check if we are allowed to alter hit_collection
  if(hit_collection->iterator_lock){
    carp(CARP_ERROR, 
         "Cannot alter hit_collection when a hit iterator is already"
         " instantiated");
    return FALSE;
  }

  switch(score_type){
  case DOTP:
    // implement later
    return FALSE;
  case XCORR:
  case LOGP_EVD_XCORR:
  case LOGP_BONF_EVD_XCORR:
  case LOGP_WEIBULL_XCORR: 
  case LOGP_BONF_WEIBULL_XCORR: 
    // LOGP_BONF_EVD_XCORR and XCORR have same order, 
    // sort the hit to decreasing XCORR order for the return
    qsort_hit(hit_collection->hit, hit_collection->hit_total, 
                (void *)compare_hit_xcorr);
    hit_collection->last_sorted = XCORR;
    return TRUE;
  case SP: 
  case LOGP_EXP_SP: 
  case LOGP_BONF_EXP_SP: 
  case LOGP_WEIBULL_SP: 
  case LOGP_BONF_WEIBULL_SP: 
  case LOGP_QVALUE_WEIBULL_XCORR: 
    // LOGP_EXP_SP and SP have same order, 
    // thus sort the hit to decreasing SP order for the return
    carp(CARP_DEBUG, "Sorting hit_collection %i", 
        hit_collection->hit_total);
    qsort_hit(hit_collection->hit, 
        hit_collection->hit_total, (void *)compare_hit_sp);
    carp(CARP_DEBUG, "Sorting hit_collection %i", 
        hit_collection->hit_total);
    hit_collection->last_sorted = SP;
    return TRUE;
  case Q_VALUE:
  case PERCOLATOR_SCORE:
    qsort_hit(hit_collection->hit, hit_collection->hit_total, 
        (void *)compare_hit_percolator_score);
    hit_collection->last_sorted = PERCOLATOR_SCORE;
    return TRUE;
  }
  return FALSE;
}

/**
 * \brief Sort a hit_collection by the given score type, grouping
 * hites by spectrum (if multiple spectra present).
 * \returns TRUE if sort is successful, else FALSE;
 */
BOOLEAN_T spectrum_sort_hit_collection(
  MATCH_COLLECTION_T* hit_collection, ///< hit collection to sort -out
  SCORER_TYPE_T score_type ///< the score type to sort by -in
  ){

  BOOLEAN_T success = FALSE;

  // check if we are allowed to alter hit_collection
  if(hit_collection->iterator_lock){
    carp(CARP_ERROR, 
         "Cannot alter hit_collection when a hit iterator is already"
         " instantiated");
    return FALSE;
  }

  switch(score_type){
  case DOTP:
    success = FALSE;
    break;

  case XCORR:
  case LOGP_EVD_XCORR:
  case LOGP_BONF_EVD_XCORR:
  case LOGP_WEIBULL_XCORR: 
  case LOGP_BONF_WEIBULL_XCORR: 
    qsort_hit(hit_collection->hit, hit_collection->hit_total,
                (void*)compare_hit_spectrum_xcorr);
    hit_collection->last_sorted = XCORR;
    success = TRUE;
    break;

  case SP: 
  case LOGP_EXP_SP: 
  case LOGP_BONF_EXP_SP: 
  case LOGP_WEIBULL_SP: 
  case LOGP_BONF_WEIBULL_SP: 
  case LOGP_QVALUE_WEIBULL_XCORR: 
    qsort_hit(hit_collection->hit, hit_collection->hit_total,
                (void*)compare_hit_spectrum_sp);
    hit_collection->last_sorted = SP;
    success = TRUE;
    break;

  case Q_VALUE:
    qsort_hit(hit_collection->hit, hit_collection->hit_total,
                (void*)compare_hit_spectrum_q_value);
    hit_collection->last_sorted = Q_VALUE;
    success = TRUE;
    break;

  case PERCOLATOR_SCORE:
    qsort_hit(hit_collection->hit, hit_collection->hit_total,
                (void*)compare_hit_spectrum_percolator_score);
    hit_collection->last_sorted = PERCOLATOR_SCORE;
    success = TRUE;
    break;


  }

  return success;
}


/**
 * keeps the top max_rank number of hites and frees the rest
 * sorts by score_type(SP, XCORR, ...)
 */
void truncate_hit_collection(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to truncate -out
  int max_rank,     ///< max number of top rank hites to keep from SP -in
  SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
  )
{
  // sort hit collection by score type
  // check if the hit collection is in the correct sorted order
  if(hit_collection->last_sorted != score_type){
    // sort hit collection by score type
    if(!sort_hit_collection(hit_collection, score_type)){
      carp(CARP_ERROR, "failed to sort hit collection");
      exit(1);
    }
  }

  // are there any hites to free?
  while(hit_collection->hit_total > max_rank){
    free_hit(hit_collection->hit[hit_collection->hit_total - 1]);
    --hit_collection->hit_total;
  }
}

/**
 * Must provide a hit_collection that is already scored, ranked for score_type
 * Rank 1, means highest score
 * \returns TRUE, if populates the hit rank in the hit collection
 */
BOOLEAN_T populate_hit_rank_hit_collection(
 MATCH_COLLECTION_T* hit_collection, ///< the hit collection to populate hit rank -out
 SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
 )
{
  // check if the hit collection is in the correct sorted order
  if(hit_collection->last_sorted != score_type){
    // sort hit collection by score type
    if(!sort_hit_collection(hit_collection, score_type)){
      carp(CARP_ERROR, "failed to sort hit collection");
      return FALSE;
    }
  }

  // set hit rank for all hit objects
  int hit_index;
  for(hit_index=0; hit_index<hit_collection->hit_total; ++hit_index){
    set_hit_rank(
        hit_collection->hit[hit_index], score_type, hit_index+1);
  }
  
  return TRUE;
}

/**
 * Create a new hit_collection by randomly sampling hites 
 * from hit_collection upto count number of hites
 * Must not free the hites
 * \returns a new hit_collection of randomly sampled hites 
 */
MATCH_COLLECTION_T* random_sample_hit_collection(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to sample -out
  int count_max ///< the number of hites to randomly select -in
  )
{
  int count_idx = 0;
  int hit_idx = 0;
  int score_type_idx = 0;
  MATCH_COLLECTION_T* sample_collection = allocate_hit_collection();
  srand(time(NULL));

  // make sure we don't sample more than the hites in the hit collection
  if (count_max < hit_collection->hit_total){
    free_hit_collection(sample_collection);
    return hit_collection;
  }

  // ranomly select hites upto count_max
  for(; count_idx < count_max; ++count_idx){
    hit_idx = ((double)rand()/((double)RAND_MAX + (double)1)) * hit_collection->hit_total;
    
    // hit_idx = rand() % hit_collection->hit_total;
    sample_collection->hit[count_idx] = hit_collection->hit[hit_idx];
    // increment pointer count of the hit object 
    increment_hit_pointer_count(sample_collection->hit[count_idx]);
  }
  
  // set total number of hites sampled
  sample_collection->hit_total = count_idx;

  sample_collection->experiment_size = hit_collection->experiment_size;

  // set scored types in the sampled hites
  for(; score_type_idx < _SCORE_TYPE_NUM;  ++score_type_idx){
    sample_collection->scored_type[score_type_idx] = hit_collection->scored_type[score_type_idx];
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
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to estimate evd parameters -in
  SCORER_TYPE_T score_type, ///< score_type to estimate EVD distribution -in
  float l_value,  ///< L value -in
  float* function,  ///< the output function value -out
  float* derivative,  ///< the output derivative value -out
  float* exponential_sum ///< the final exponential array sum -out
  )
{
  int idx = 0;
  float* exponential = (float*)mycalloc(hit_collection->hit_total, sizeof(float));
  float numerator = 0;
  float second_numerator = 0;
  float score = 0;
  float denominator = 0;
  float score_sum = 0;
  HIT_T** hites = hit_collection->hit;

  // iterate over the hites to calculate numerator, exponential value, denominator
  for(; idx < hit_collection->hit_total; ++idx){
    score = get_hit_score(hites[idx], score_type);
    exponential[idx] = exp(-l_value * score);
    numerator += (exponential[idx] * score);
    denominator += exponential[idx];
    score_sum += score;
    second_numerator += (score * score * exponential[idx]);
  }

  // assign function value
  *function = (1.0 / l_value) - (score_sum / hit_collection->hit_total) 
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
 * Randomly samples max_count peptides from the peptide distribution and try to esitimate the Xcorr distribution of the the entire peptide distribution 
 * from the sampled peptide distribution. Populates the two EVD parameters mu, lambda in the hit_collection.
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
BOOLEAN_T estimate_evd_parameters(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to estimate evd parameters -out
  int sample_count, ///< the number of peptides to sample from the hit_collection -in
  SCORER_TYPE_T score_type, ///< score_type to estimate EVD distribution -in
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  int charge       ///< the charge of the spectrum -in
  )
{
  // randomly sample from hit collection
  MATCH_COLLECTION_T* sample_collection = random_sample_hit_collection(hit_collection, sample_count);
  float l_value = 1;
  float f = 0.0;
  float f_prime = 0.0;
  float epsilon = 0.001;
  float exponential_sum = 0;
  int max_iterations = 10000;
  int idx = 0;

  // print info
  carp(CARP_INFO, "Estimate EVD parameters, sample count: %d", sample_count);
  
  // first score the sample hit_collection
  if(score_type == XCORR){
    if(!score_hit_collection_xcorr(sample_collection, spectrum, charge)){
      carp(CARP_ERROR, "failed to score hit collection for XCORR");
    }
  }
  // FIXME Add different scoring if needed
  // ...

  // estimate the EVD parameters
  for(; idx < max_iterations; ++idx){
    constraint_function(sample_collection, score_type, l_value, 
                                       &f, &f_prime, &exponential_sum);

    if(fabsf(f) < epsilon){
      break;
    }
    else{
      l_value = l_value - f / f_prime;
    }
    
    // failed to converge error..
    if(idx >= max_iterations){
      carp(CARP_ERROR, "Root finding failed to converge.");
      free_hit_collection(sample_collection);
      return FALSE;
    }
  }
  
  // Calculate best value of position parameter from best value of 
  // scale parameter.
  hit_collection->mu = -1.0 / l_value * logf(1.0 / sample_count * exponential_sum);
  hit_collection->l_value = l_value;
    
  // free up sampled hit_collection 
  free_hit_collection(sample_collection);
  
  // DEBUG
  // carp(CARP_DEBUG, "mu: %.5f, L: %.5f", hit_collection->mu, hit_collection->l_value);
  return TRUE;
}

// TODO why does xcorr need spectrum and charge but sp scoring function doesn't?
// TODO score_hit_collection_sp should probably not take an iterator?
// TODO change score_hit_collection* to single routine score_hit_collection
// TODO sample_count should probably not be an explicit parameter to the
// fitting code (should be like e.g. fraction-top-scores-to-fit)

/**
 * For the #top_count ranked peptides, calculate the Weibull parameters
 *\returns TRUE, if successfully calculates the Weibull parameters
 */
#define MIN_XCORR_SHIFT -5.0
#define MAX_XCORR_SHIFT  5.0
#define XCORR_SHIFT 0.05
#define MIN_SP_SHIFT -100.0
#define MAX_SP_SHIFT 300.0
#define SP_SHIFT 5.0

BOOLEAN_T estimate_weibull_parameters(
  MATCH_COLLECTION_T* hit_collection, 
  ///< the hit collection for which to estimate weibull parameters -out
  SCORER_TYPE_T score_type,
  int sample_count,
  SPECTRUM_T* spectrum,
  int charge
  )
{
  carp(CARP_INFO, "Estimating weibull params");
  MATCH_COLLECTION_T* sample_collection = hit_collection;

  if (sample_count != 0){
    sample_collection = 
      random_sample_hit_collection(hit_collection, sample_count);
  }

  // how many things are we going to fit. We may want to just fit to the
  // tail, thus the distinction between total* and fit*
  int total_data_points = sample_collection->hit_total;
  carp(CARP_INFO, "Stat: Total hites: %i\n", total_data_points);
  int fit_data_points = total_data_points;

  // less than 0.0 or 0 indicates use all peptides
  double fraction_to_fit = get_double_parameter("fraction-top-scores-to-fit");
  int number_to_fit = get_int_parameter("number-top-scores-to-fit");
  carp(CARP_INFO, "Stat: Number hites to fit: %i\n", number_to_fit);
  if (fraction_to_fit > -0.5){
    assert(fraction_to_fit <= 1.0);
    fit_data_points = (int)(total_data_points * fraction_to_fit);
  } else if (number_to_fit > -1 ){
    fit_data_points = number_to_fit < total_data_points ? 
        number_to_fit : total_data_points;
  }

  carp(CARP_INFO, "Estimate Weibull parameters, count: %d", fit_data_points);
  
  // first score the sample hit_collection
  if(score_type == XCORR){
    if(!score_hit_collection_xcorr(sample_collection, spectrum, charge)){
      carp(CARP_ERROR, "failed to score hit collection for XCORR");
    }
  } else if (score_type == SP){
    // FIXME assumes scored by SP already
    ;
  }

  // check if the hit collection is in the correct sorted order
  if(sample_collection->last_sorted != score_type){
    // sort hit collection by score type
    if(!sort_hit_collection(sample_collection, score_type)){
      carp(CARP_ERROR, "failed to sort hit collection");
      exit(1);
    }
  }

  // implementation of Weibull distribution parameter estimation from 
  // http:// www.chinarel.com/onlincebook/LifeDataWeb/rank_regression_on_y.htm
  
  int idx;
  float* data   = calloc(sizeof(float) , total_data_points);
  for(idx=0; idx < total_data_points; idx++){
    float score = get_hit_score(sample_collection->hit[idx], score_type);
    data[idx] = score;
  }

  float correlation = 0.0;
  if (score_type == XCORR){
    fit_three_parameter_weibull(data, fit_data_points, total_data_points,
      MIN_XCORR_SHIFT, MAX_XCORR_SHIFT, XCORR_SHIFT, 
      &(hit_collection->eta), &(hit_collection->beta),
      &(hit_collection->shift), &correlation);
  } else if (score_type == SP){
    fit_three_parameter_weibull(data, fit_data_points, total_data_points,
      MIN_SP_SHIFT, MAX_SP_SHIFT, SP_SHIFT, 
      &(hit_collection->eta), &(hit_collection->beta), 
      &(hit_collection->shift), &correlation);
  }
  carp(CARP_INFO, "Correlation: %.6f\nEta: %.6f\nBeta: %.6f\nShift: %.6f\n", 
      correlation, hit_collection->eta, hit_collection->beta,
      hit_collection->shift);
  
  return TRUE;
}


/**
 * For the #top_count SP ranked peptides, calculate the mean for which the
 * #top_ranked peptide score is set to 0, thus scaling the SP scores.
 *\returns TRUE, if successfully calculates the EXP_SP parameters
 */
BOOLEAN_T estimate_exp_sp_parameters(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to estimate evd parameters -out
  int top_count ///< the number of top SP peptides to use for the hit_collection -in
  )
{
  float top_sp_score = 0.0;
  float base_score = 0.0;
  int count  = 0;
  
  // sort hit collection by SP
  // check if the hit collection is in the correct sorted order
  if(hit_collection->last_sorted != SP){
    // sort hit collection by score type
    if(!sort_hit_collection(hit_collection, SP)){
      carp(CARP_ERROR, "failed to sort hit collection");
      exit(1);
    }
  }
  
  // adjust the number of top ranked peptides to sample
  // because the the total number of peptides are less than top_count
  if(top_count > hit_collection->hit_total){
    top_count = hit_collection->hit_total;
    carp(CARP_INFO, "");
  }
  
  // set the base score to which score is set to 0
  base_score = get_hit_score(hit_collection->hit[top_count-1], SP);
  
  // compile the scores
  while(count < top_count){
    top_sp_score += get_hit_score(hit_collection->hit[count], SP);
    ++count;
  }
  
  hit_collection->sp_scores_mean = ((top_sp_score) / count - base_score);
  hit_collection->base_score_sp = base_score;
  hit_collection->top_fit_sp = top_count;
  
  return TRUE;
}

/**
 * Preliminary scoring methods
 *
 * Preliminary scoring methods create new hit objects.
 * Also, they must set the hit object null peptide boolean value.
 * To get the peptide sequence the order should always be
 * 1. get peptide object from generate_peptides or some other sources.
 * 2. create new hit object.
 * 3. set hit object null_peptide value.
 * 4. set peptide as hit objects peptide.
 * 5. Finally get peptide sequence through get_hit_sequence method
 */

/**
 * Preliminary scoring method: 
 * creates new hit objects, and sets them as to which they are null peptide or not
 *
 * scores the hit_collection, the score type SP
 * Assumes this is the first time scoring with this score_collection,
 * thus, prior number hit object is 0.
 * the routine will use generate_peptides for each peptide will create a hit
 * that maps the peptide to the spectrum.
 * If the score has already been computed simply returns TRUE 
 *\returns  TRUE, if successfully populates the sp score hites in the hit_collection
 */
BOOLEAN_T score_hit_collection_sp(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to score -out
  SPECTRUM_T* spectrum, ///< the spectrum to hit peptides -in
  int charge,       ///< the charge of the spectrum -in
  GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator ///< the peptide iterator to score -in
  )
{
  
  // is this a empty collection?
  if(hit_collection->hit_total != 0){
    carp(CARP_ERROR, "must start with empty hit collection");
    return FALSE;
  }
  
  // set ion constraint to sequest settings
  ION_CONSTRAINT_T* ion_constraint = 
    new_ion_constraint_sequest_sp(charge); 
  
  // create new scorer
  SCORER_T* scorer = new_scorer(SP);  

  char* peptide_sequence = NULL;
  HIT_T* hit = NULL;
  float score = 0;
  PEPTIDE_T* peptide = NULL;  

  // create a generic ion_series that will be reused for each peptide sequence
  ION_SERIES_T* ion_series = new_ion_series_generic(ion_constraint, charge);  
  
  // iterate over all peptides
  carp(CARP_DEBUG, "Iterating over peptides to score Sp");
  while(generate_peptides_iterator_has_next(peptide_iterator)){
    peptide = generate_peptides_iterator_next(peptide_iterator);

    // create a new hit
    hit = new_hit();

    // set hit if it is to be generated as null peptide hit
    set_hit_null_peptide(hit, hit_collection->null_peptide_collection);
    
    // now set peptide and spectrum
    set_hit_peptide(hit, peptide);

    set_hit_spectrum(hit, spectrum);
    
    // get peptide sequence    
    peptide_sequence = get_hit_sequence(hit);
    
    // update ion_series for the peptide instance    
    update_ion_series(ion_series, peptide_sequence);

    // now predict ions for this peptide
    predict_ions(ion_series);
    
    // calculates the Sp score
    score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);

    // increment the total sp score
    hit_collection->sp_scores_mean += score;
        
    // set score in hit
    set_hit_score(hit, SP, score);

    // set b_y_ion_hit field
    set_hit_b_y_ion_info(hit, scorer);
    
    // check if enough space for peptide hit
    if(hit_collection->hit_total >= _MAX_NUMBER_PEPTIDES){
      carp(CARP_ERROR, "peptide count of %i exceeds max hit limit: %d", 
          hit_collection->hit_total, _MAX_NUMBER_PEPTIDES);
      // free heap
      free(peptide_sequence);
      free_ion_series(ion_series);
      free_scorer(scorer);
      free_ion_constraint(ion_constraint);

      return FALSE;
    }
    
    // add a new hit to array
    hit_collection->hit[hit_collection->hit_total] = hit;
    
    // increment total hit count
    ++hit_collection->hit_total;

    // print total peptides scored so far
    if(hit_collection->hit_total % 10000 == 0){
      carp(CARP_INFO, "scored peptide for sp: %d", 
          hit_collection->hit_total);
    }
    
    free(peptide_sequence);
  }
  // free ion_series now that we are done iterating over all peptides
  free_ion_series(ion_series);
  
  // calculate the final sp score mean
  hit_collection->sp_scores_mean /= hit_collection->hit_total;
  
  // total peptide experiment sample size
  hit_collection->experiment_size = hit_collection->hit_total;

  // print total peptides scored so far
  carp(CARP_DEBUG, "Total peptide scored for sp: %d", 
      hit_collection->hit_total);

  if (hit_collection->hit_total)
  
  // free heap
  free_scorer(scorer);
  free_ion_constraint(ion_constraint);
    
  // now hit_collection is sorted, populate the rank of each hit object
  if(!populate_hit_rank_hit_collection(hit_collection, SP)){
    carp(CARP_ERROR, "failed to populate rank for SP in hit_collection");
    free_hit_collection(hit_collection);
    exit(1);
  }
  
  // yes, we have now scored for the hit-mode: SP
  hit_collection->scored_type[SP] = TRUE;
  
  return TRUE;
}


/**
 * Main scoring methods
 * 
 * Unlike preliminary scoring methods only updates existing mathc objects with new scores.
 * In all cases should get peptide sequence only through get_hit_sequence method.
 */


/**
 * The hit collection must be scored under SP first
 * \returns TRUE, if successfully scores hites for LOGP_EXP_SP
 */
BOOLEAN_T score_hit_collection_logp_exp_sp(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to score -out
  int peptide_to_score ///< the number of top ranked sp scored peptides to score for logp_exp_sp -in
  )
{
  int hit_idx = 0;
  float score = 0;
  HIT_T* hit = NULL;
  
  // has the score type been populated in hit collection?
  if(!hit_collection->scored_type[SP]){
    carp(CARP_ERROR, "the collection must be scored by SP first before LOGP_EXP_SP");
    exit(1);
  }

  // sort by SP if not already sorted.
  // This enables to identify the top ranked SP scoring peptides
  if(hit_collection->last_sorted != SP){
    // sort hit collection by score type
    if(!sort_hit_collection(hit_collection, SP)){
      carp(CARP_ERROR, "failed to sort hit collection by SP");
      free_hit_collection(hit_collection);
      exit(1);
    }
  }
  
  // we are string xcorr!
  carp(CARP_INFO, "start scoring for LOGP_EXP_SP");

  // iterate over all hites to score for LOGP_EXP_SP
  while(hit_idx < hit_collection->hit_total && hit_idx < peptide_to_score){
    hit = hit_collection->hit[hit_idx];
    // scale the SP score by the base score found from estimate_exp_sp_parameters routine
    score = score_logp_exp_sp((get_hit_score(hit, SP) - hit_collection->base_score_sp), hit_collection->sp_scores_mean);
    
    // set all fields in hit
    set_hit_score(hit, LOGP_EXP_SP, score);
    ++hit_idx;
  }
  
  // we are done
  carp(CARP_INFO, "total peptides scored for LOGP_EXP_SP: %d", hit_idx);

  // hit_collection is not populate with the rank of LOGP_EXP_SP, becuase the SP rank is  identical to the LOGP_EXP_SP rank
  
  // yes, we have now scored for the hit-mode: LOGP_EXP_SP
  hit_collection->scored_type[LOGP_EXP_SP] = TRUE;
  
  return TRUE;
}


/**
 * The hit collection must be scored under SP first
 * \returns TRUE, if successfully scores hites for LOGP_BONF_EXP_SP
 */
BOOLEAN_T score_hit_collection_logp_bonf_exp_sp(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to score -out
  int peptide_to_score ///< the number of top ranked sp scored peptides to score for logp_bonf_exp_sp -in
  )
{
  int hit_idx = 0;
  double score = 0;
  HIT_T* hit = NULL;
  
  // has the score type been populated in hit collection?
  if(!hit_collection->scored_type[SP]){
    carp(CARP_ERROR, "the collection must be scored by SP first before LOGP_EXP_SP");
    exit(1);
  }

  // sort by SP if not already sorted.
  // This enables to identify the top ranked SP scoring peptides
  if(hit_collection->last_sorted != SP){
    // sort hit collection by score type
    if(!sort_hit_collection(hit_collection, SP)){
      carp(CARP_ERROR, "failed to sort hit collection by SP");
      free_hit_collection(hit_collection);
      exit(1);
    }
  }
  
  // we are string xcorr!
  carp(CARP_INFO, "start scoring for LOGP_BONF_EXP_SP");

  // iterate over all hites to score for LOGP_BONF_EXP_SP
  while(hit_idx < hit_collection->hit_total && hit_idx < peptide_to_score){
    hit = hit_collection->hit[hit_idx];
    // scale the SP score by the base score found from estimate_exp_sp_parameters routine
    score = score_logp_bonf_exp_sp((get_hit_score(hit, SP) - hit_collection->base_score_sp), hit_collection->sp_scores_mean, hit_collection->experiment_size);
    
    // set all fields in hit
    set_hit_score(hit, LOGP_BONF_EXP_SP, score);
    ++hit_idx;
  }
  
  // we are done
  carp(CARP_INFO, "total peptides scored for LOGP_BONF_EXP_SP: %d", hit_idx);
    
  // hit_collection is not populate with the rank of LOGP_BONF_EXP_SP, becuase the SP rank is  identical to the LOGP_EXP_SP rank
  
  // yes, we have now scored for the hit-mode: LOGP_BONF_EXP_SP
  hit_collection->scored_type[LOGP_BONF_EXP_SP] = TRUE;
  
  return TRUE;
}

/**
 * The hit collection must be scored under SP first
 * \returns TRUE, if successfully scores hites for LOGP_WEIBULL_SP
 */
BOOLEAN_T score_hit_collection_logp_weibull_sp(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to score -out
  int peptide_to_score ///< the number of top ranked sp scored peptides to score for logp_weibull_sp -in
  )
{
  int hit_idx = 0;
  float score = 0;
  HIT_T* hit = NULL;
  
  // has the score type been populated in hit collection?
  if(!hit_collection->scored_type[SP]){
    carp(CARP_ERROR, "the collection must be scored by SP first before LOGP_WEIBULL_SP");
    exit(1);
  }

  // sort by SP if not already sorted.
  // This enables to identify the top ranked SP scoring peptides
  if(hit_collection->last_sorted != SP){
    // sort hit collection by score type
    if(!sort_hit_collection(hit_collection, SP)){
      carp(CARP_ERROR, "failed to sort hit collection by SP");
      free_hit_collection(hit_collection);
      exit(1);
    }
  }
  
  // we are string xcorr!
  carp(CARP_INFO, "start scoring for LOGP_WEIBULL_SP");

  // iterate over all hites to score for LOGP_WEIBULL_SP
  while(hit_idx < hit_collection->hit_total && hit_idx < peptide_to_score){
    hit = hit_collection->hit[hit_idx];
    // scale the SP score by the base score found from estimate_weibull_sp_parameters routine
    score = score_logp_weibull(get_hit_score(hit, SP),
          hit_collection->eta, hit_collection->beta);
    
    // set all fields in hit
    set_hit_score(hit, LOGP_WEIBULL_SP, score);
    ++hit_idx;
  }
  
  // we are done
  carp(CARP_INFO, "total peptides scored for LOGP_WEIBULL_SP: %d", hit_idx);

  // hit_collection is not populate with the rank of LOGP_WEIBULL_SP, becuase the SP rank is  identical to the LOGP_WEIBULL_SP rank
  
  // yes, we have now scored for the hit-mode: LOGP_WEIBULL_SP
  hit_collection->scored_type[LOGP_WEIBULL_SP] = TRUE;
  
  return TRUE;
}

/**
 * The hit collection must be scored under XCORR first
 * \returns TRUE, if successfully scores hites for LOGP_WEIBULL_XCORR
 */
BOOLEAN_T score_hit_collection_logp_weibull_xcorr(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to score -out
  int peptide_to_score ///< the number of top ranked xcorr scored peptides to score for logp_weibull_xcorr -in
  )
{
  int hit_idx = 0;
  float score = 0;
  HIT_T* hit = NULL;
  
  // has the score type been populated in hit collection?
  if(!hit_collection->scored_type[XCORR]){
    carp(CARP_ERROR, "the collection must be scored by XCORR first before LOGP_WEIBULL_XCORR");
    exit(1);
  }

  // sort by XCORR if not already sorted.
  // This enables to identify the top ranked XCORR scoring peptides
  if(hit_collection->last_sorted != XCORR){
    // sort hit collection by score type
    if(!sort_hit_collection(hit_collection, XCORR)){
      carp(CARP_ERROR, "failed to sort hit collection by XCORR");
      free_hit_collection(hit_collection);
      exit(1);
    }
  }
  
  // we are string xcorr!
  carp(CARP_INFO, "start scoring for LOGP_WEIBULL_XCORR");

  // iterate over all hites to score for LOGP_WEIBULL_XCORR
  while(hit_idx < hit_collection->hit_total && hit_idx < peptide_to_score){
    hit = hit_collection->hit[hit_idx];
    // scale the XCORR score by the base score found from estimate_weibull_parameters routine
    score = score_logp_weibull(get_hit_score(hit, XCORR),
          hit_collection->eta, hit_collection->beta);
    // set all fields in hit
    set_hit_score(hit, LOGP_WEIBULL_XCORR, score);
    ++hit_idx;
  }
  
  // we are done
  carp(CARP_INFO, "total peptides scored for LOGP_WEIBULL_XCORR: %d", hit_idx);

  // hit_collection is not populate with the rank of LOGP_WEIBULL_XCORR, becuase the XCORR rank is  identical to the LOGP_WEIBULL_XCORR rank
  
  // yes, we have now scored for the hit-mode: LOGP_WEIBULL_XCORR
  hit_collection->scored_type[LOGP_WEIBULL_XCORR] = TRUE;
  
  return TRUE;
}


/**
 * The hit collection must be scored under XCORR first
 * \returns TRUE, if successfully scores hites for LOGP_BONF_WEIBULL_XCORR
 */
BOOLEAN_T score_hit_collection_logp_bonf_weibull_xcorr(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to score -out
  int peptide_to_score ///< the number of top ranked xcorr scored peptides to score for logp_bonf_weibull_xcorr -in
  )
{
  int hit_idx = 0;
  double score = 0;
  HIT_T* hit = NULL;
  
  // has the score type been populated in hit collection?
  if(!hit_collection->scored_type[XCORR]){
    carp(CARP_ERROR, "the collection must be scored by XCORR first before LOGP_WEIBULL_XCORR");
    exit(1);
  }

  // sort by XCORR if not already sorted.
  // This enables to identify the top ranked XCORR scoring peptides
  if(hit_collection->last_sorted != XCORR){
    // sort hit collection by score type
    if(!sort_hit_collection(hit_collection, XCORR)){
      carp(CARP_ERROR, "failed to sort hit collection by XCORR");
      free_hit_collection(hit_collection);
      exit(1);
    }
  }
  
  // we are string xcorr!
  carp(CARP_INFO, "start scoring for LOGP_BONF_WEIBULL_XCORR");

  // iterate over all hites to score for LOGP_BONF_WEIBULL_XCORR
  while(hit_idx < hit_collection->hit_total && hit_idx < peptide_to_score){
    hit = hit_collection->hit[hit_idx];
    // scale the XCORR score by the params from estimate_weibull_parameters
    score = score_logp_bonf_weibull(get_hit_score(hit, XCORR),
          hit_collection->eta, hit_collection->beta, 
          hit_collection->shift, hit_collection->experiment_size);
    
    // set all fields in hit
    set_hit_score(hit, LOGP_BONF_WEIBULL_XCORR, score);
    ++hit_idx;
  }
  
  // we are done
  carp(CARP_INFO, "total peptides scored for LOGP_BONF_WEIBULL_XCORR: %d", hit_idx);
    
  // hit_collection is not populate with the rank of LOGP_BONF_WEIBULL_XCORR, becuase the XCORR rank is  identical to the LOGP_WEIBULL_XCORR rank
  
  // yes, we have now scored for the hit-mode: LOGP_BONF_WEIBULL_XCORR
  hit_collection->scored_type[LOGP_BONF_WEIBULL_XCORR] = TRUE;
  
  return TRUE;
}



/**
 * The hit collection must be scored under SP first
 * \returns TRUE, if successfully scores hites for LOGP_BONF_WEIBULL_SP
 */
BOOLEAN_T score_hit_collection_logp_bonf_weibull_sp(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to score -out
  int peptide_to_score ///< the number of top ranked sp scored peptides to score for logp_bonf_weibull_sp -in
  )
{
  int hit_idx = 0;
  float score = 0;
  HIT_T* hit = NULL;
  
  // has the score type been populated in hit collection?
  if(!hit_collection->scored_type[SP]){
    carp(CARP_ERROR, "the collection must be scored by SP first before LOGP_WEIBULL_SP");
    exit(1);
  }

  // sort by SP if not already sorted.
  // This enables to identify the top ranked SP scoring peptides
  if(hit_collection->last_sorted != SP){
    // sort hit collection by score type
    if(!sort_hit_collection(hit_collection, SP)){
      carp(CARP_ERROR, "failed to sort hit collection by SP");
      free_hit_collection(hit_collection);
      exit(1);
    }
  }
  
  // we are string xcorr!
  carp(CARP_INFO, "start scoring for LOGP_BONF_WEIBULL_SP");

  // iterate over all hites to score for LOGP_BONF_WEIBULL_SP
  while(hit_idx < hit_collection->hit_total && hit_idx < peptide_to_score){
    hit = hit_collection->hit[hit_idx];
    // scale the SP score by the params from estimate_weibull_sp_parameters
    score = score_logp_bonf_weibull(get_hit_score(hit, SP),
          hit_collection->eta, hit_collection->beta, 
          hit_collection->shift, hit_collection->experiment_size);
    
    // set all fields in hit
    set_hit_score(hit, LOGP_BONF_WEIBULL_SP, score);
    ++hit_idx;
  }
  
  // we are done
  carp(CARP_INFO, "total peptides scored for LOGP_BONF_WEIBULL_SP: %d", hit_idx);
    
  // hit_collection is not populate with the rank of LOGP_BONF_WEIBULL_SP, becuase the SP rank is  identical to the LOGP_WEIBULL_SP rank
  
  // yes, we have now scored for the hit-mode: LOGP_BONF_WEIBULL_SP
  hit_collection->scored_type[LOGP_BONF_WEIBULL_SP] = TRUE;
  
  return TRUE;
}


/**
 * Assumes that hit collection was scored under SP first
 * Creates an ion constraint, a scorer, an ion series.  Modifies the
 * hites in the collection by setting the score.
 * \returns TRUE, if successfully scores hites for xcorr
 * \callgraph
 */
BOOLEAN_T score_hit_collection_xcorr(
  MATCH_COLLECTION_T* hit_collection, ///<the hit collection to score -out
  SPECTRUM_T* spectrum, ///< the spectrum to hit peptides -in
  int charge       ///< the charge of the spectrum -in
  )
{
  HIT_T* hit = NULL;
  char* peptide_sequence = NULL;  
  float score = 0;
  
  /*
  // is this a empty collection?
  if(hit_collection->hit_total > 0){
    carp(CARP_ERROR, "must start with SP scored hit collection");
    return FALSE;
  }
  */
  
  // set ion constraint to sequest settings
  ION_CONSTRAINT_T* ion_constraint = 
    new_ion_constraint_sequest_xcorr(charge); 
  
  // create new scorer
  SCORER_T* scorer = new_scorer(XCORR);  

  // create a generic ion_series that will be reused for each peptide sequence
  ION_SERIES_T* ion_series = new_ion_series_generic(ion_constraint, charge);  
  
  // we are scoring xcorr!
  carp(CARP_DEBUG, "Start scoring for XCORR");

  // iterate over all hites to score for xcorr
  int hit_idx;
  for(hit_idx=0; hit_idx < hit_collection->hit_total; ++hit_idx){
    hit = hit_collection->hit[hit_idx];
    peptide_sequence = get_hit_sequence(hit);
    
    // update ion_series for the peptide instance    
    update_ion_series(ion_series, peptide_sequence);
    
    // now predict ions
    predict_ions(ion_series);
    
    // calculates the Xcorr score
    score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);
    char* decoy = "target";
    if (get_hit_null_peptide(hit)==TRUE){
      decoy = "decoy";
    }
    carp(CARP_DETAILED_DEBUG, "Spectrum %i vs. %s peptide %s = %.6f", 
      get_spectrum_first_scan(spectrum), decoy, peptide_sequence, score);

    // set all fields in hit
    set_hit_score(hit, XCORR, score);
    
    // free heap
    free(peptide_sequence);   
  }  

  // free ion_series now that we are done iterating over all peptides
  free_ion_series(ion_series);

  // we scored xcorr!
  carp(CARP_DEBUG, "Total peptides scored for XCORR: %d", hit_idx);

  // free heap
  free_scorer(scorer);
  free_ion_constraint(ion_constraint);

  // sort hit collection by score type
  if(!sort_hit_collection(hit_collection, XCORR)){
    carp(CARP_FATAL, "Failed to sort hit collection by Xcorr");
    exit(1);
  }
  
  // now the hit_collection is sorted, update the rank of each hit object
  if(!populate_hit_rank_hit_collection(hit_collection, XCORR)){
    carp(CARP_FATAL, "Failed to populate hit rank in hit_collection");
    free_hit_collection(hit_collection);
    exit(1);
  }

  // calculate deltaCn value (difference between best and 2nd best score)
  if(hit_collection->hit_total > 1){
    hit_collection->delta_cn = 
      get_hit_score(hit_collection->hit[0], XCORR) -
      get_hit_score(hit_collection->hit[1], XCORR);
  }
  else{
    // set to very small number
    hit_collection->delta_cn = 0.000001;
  }
  
  // yes, we have now scored for the hit-mode: XCORR
  hit_collection->scored_type[XCORR] = TRUE;

  return TRUE;
}


/**
 * The hit collection must be scored under Xcorr first
 * \returns TRUE, if successfully scores hites for LOGP_EVD_XCORR
 */
BOOLEAN_T score_hit_collection_logp_evd_xcorr(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to score -out
  int peptide_to_score ///< the number of top ranked xcorr scored peptides to score for logp_evd_xcorr -in
  )
{
  int hit_idx = 0;
  float score = 0;
  HIT_T* hit = NULL;
  
  // has the score type been populated in hit collection?
  if(!hit_collection->scored_type[XCORR]){
    carp(CARP_ERROR, "the collection must be scored by XCORR first before LOGP_EVD_XCORR");
    exit(1);
  }

  // sort by XCORR if not already sorted.
  // This enables to identify the top ranked XCORR scoring peptides
  if(hit_collection->last_sorted != XCORR){
    // sort hit collection by score type
    if(!sort_hit_collection(hit_collection, XCORR)){
      carp(CARP_ERROR, "failed to sort hit collection by XCORR");
      free_hit_collection(hit_collection);
      exit(1);
    }
  }
  
  // we are starting LOGP_EVD_XCORR!
  carp(CARP_INFO, "start scoring for LOGP_EVD_XCORR");

  // iterate over all hites to score for LOGP_EVD_XCORR
  while(hit_idx < hit_collection->hit_total && hit_idx < peptide_to_score){
    hit = hit_collection->hit[hit_idx];
    score = score_logp_evd_xcorr(get_hit_score(hit, XCORR), hit_collection->mu, hit_collection->l_value);
    
    // set all fields in hit
    set_hit_score(hit, LOGP_EVD_XCORR, score);
    ++hit_idx;
  }
  
  // we are done
  carp(CARP_INFO, "total peptides scored for LOGP_EVD_XCORR: %d", hit_idx);

  // hit_collection is not populate with the rank of LOGP_EVD_XCORR, 
  // becuase the XCORR rank is  identical to the LOGP_EVD_XCORR rank
  
  // yes, we have now scored for the hit-mode: LOGP_EVD_XCORR
  hit_collection->scored_type[LOGP_EVD_XCORR] = TRUE;
  
  return TRUE;
}

/**
 * The hit collection must be scored under Xcorr first
 * \returns TRUE, if successfully scores hites for LOGP_BONF_EVD_XCORR
 */
BOOLEAN_T score_hit_collection_logp_bonf_evd_xcorr(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to score -out
  int peptide_to_score ///< the number of top ranked xcorr scored peptides to score for logp_bonf_evd_xcorr -in
  )
{
  int hit_idx = 0;
  float score = 0;
  HIT_T* hit = NULL;
  
  // has the score type been populated in hit collection?
  if(!hit_collection->scored_type[XCORR]){
    carp(CARP_ERROR, "the collection must be scored by XCORR first before LOGP_BONF_EVD_XCORR");
    exit(1);
  }

  // sort by XCORR if not already sorted.
  // This enables to identify the top ranked XCORR scoring peptides
  if(hit_collection->last_sorted != XCORR){
    // sort hit collection by score type
    if(!sort_hit_collection(hit_collection, XCORR)){
      carp(CARP_ERROR, "failed to sort hit collection by XCORR");
      free_hit_collection(hit_collection);
      exit(1);
    }
  }
  
  // we are starting LOGP_BONF_EVD_XCORR!
  carp(CARP_INFO, "start scoring for LOGP_BONF_EVD_XCORR");

  // iterate over all hites to score for LOGP_BONF_EVD_XCORR
  while(hit_idx < hit_collection->hit_total && hit_idx < peptide_to_score){
    hit = hit_collection->hit[hit_idx];
    score = score_logp_bonf_evd_xcorr(get_hit_score(hit, XCORR), hit_collection->mu, hit_collection->l_value, hit_collection->experiment_size);
    
    // set all fields in hit
    set_hit_score(hit, LOGP_BONF_EVD_XCORR, score);
    ++hit_idx;
  }
  
  // we are done
  carp(CARP_INFO, "total peptides scored for LOGP_BONF_EVD_XCORR: %d", hit_idx);

  // hit_collection is not populate with the rank of LOGP_BONF_EVD_XCORR, 
  // becuase the XCORR rank is  identical to the LOGP_BONF_EVD_XCORR rank
  
  // yes, we have now scored for the hit-mode: LOGP_BONF_EVD_XCORR
  hit_collection->scored_type[LOGP_BONF_EVD_XCORR] = TRUE;
  
  return TRUE;
}

/**
 * hit_collection get, set method
 */

/**
 *\returns TRUE, if the hit collection has been scored by score_type
 */
BOOLEAN_T get_hit_collection_scored_type(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to iterate -in
  SCORER_TYPE_T score_type ///< the score_type (MATCH_SP, MATCH_XCORR) -in
  )
{
  return hit_collection->scored_type[score_type];
}

/**
 * sets the score_type to value
 */
void set_hit_collection_scored_type(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to iterate -in
  SCORER_TYPE_T score_type, ///< the score_type (MATCH_SP, MATCH_XCORR) -in
  BOOLEAN_T value
  )
{
  hit_collection->scored_type[score_type] = value;
}

/**
 *\returns TRUE, if there is a  hit_iterators instantiated by hit collection 
 */
BOOLEAN_T get_hit_collection_iterator_lock(
  MATCH_COLLECTION_T* hit_collection ///< working hit collection -in
  )
{
  return hit_collection->iterator_lock;
}

/**
 *\returns the total hit objects avaliable in current hit_collection
 */
int get_hit_collection_hit_total(
  MATCH_COLLECTION_T* hit_collection ///< working hit collection -in
  )
{
  return hit_collection->hit_total;
}

/**
 *\returns the total peptides searched in the experiment in hit_collection
 */
int get_hit_collection_experimental_size(
  MATCH_COLLECTION_T* hit_collection ///< working hit collection -in
  )
{
  return hit_collection->experiment_size;
}

/**
 *\returns the top peptide count used in the logp_exp_sp in hit_collection
 */
int get_hit_collection_top_fit_sp(
  MATCH_COLLECTION_T* hit_collection ///< working hit collection -in
  )
{
  return hit_collection->top_fit_sp;
}

/**
 *\returns the charge of the spectrum that the hit collection was created
 */
int get_hit_collection_charge(
  MATCH_COLLECTION_T* hit_collection ///< working hit collection -in
  )
{
  return hit_collection->charge;
}

/**
 * Must have been scored by Xcorr, returns error if not scored by Xcorr
 *\returns the delta cn value(difference in top and second ranked Xcorr values)
 */
float get_hit_collection_delta_cn(
  MATCH_COLLECTION_T* hit_collection ///< working hit collection -in
  )
{
  // Check if xcorr value has been scored, thus delta cn value is valid
  if(hit_collection->scored_type[XCORR]){
    return hit_collection->delta_cn;
  }
  else{
    carp(CARP_ERROR, "must score hit_collection with XCORR to get delta cn value");
    return 0.0;
  }
}

/**
 * \brief Names and opens the correct number of binary psm files.
 *
 * Takes the values of hit-output-folder, ms2 filename (soon to be
 * named output file), overwrite, and number-decoy-set from parameter.c.
 * Exits with error if can't create new requested directory or if
 * can't create any of the psm files.
 * REPLACES: spectrum_collection::get_spectrum_collection_psm_result_filenames
 *
 * \returns An array of filehandles to the newly opened files
 */
FILE** create_psm_files(){

  int decoy_sets = get_int_parameter("number-decoy-set");
  int total_files = decoy_sets +1;
  // create FILE* array to return
  FILE** file_handle_array = (FILE**)mycalloc(total_files, sizeof(FILE*));
  int file_idx = 0;

  // Create null pointers if no binary output called for
  if( SQT_OUTPUT == get_output_type_parameter("output-mode") ){
    carp(CARP_DEBUG, "SQT mode: return empty array of file handles");
    return file_handle_array;
  }

  carp(CARP_DEBUG, "Opening %d new psm files", total_files);

  char* output_directory =get_string_parameter_pointer("hit-output-folder");

  // create the output folder if it doesn't exist
  if(access(output_directory, F_OK)){
    if(mkdir(output_directory, S_IRWXU+S_IRWXG+S_IRWXO) != 0){
      carp(CARP_FATAL, "Failed to create output directory %s", 
           output_directory);
      exit(1);
    }
  }

  // get ms2 file for naming result file
  //TODO change to output filename as argument, force .csm extension
  //     add _decoy1.csm
  //char* base_filename = get_string_parameter_pointer("ms2 file");
  char* ms2_filename = get_string_parameter_pointer("ms2 file");
  //char** filename_path_array = parse_filename_path(base_filename);
  char** filename_path_array = 
    parse_filename_path_extension(ms2_filename, ".ms2");
  if( filename_path_array[1] == NULL ){
    filename_path_array[1] = ".";
  }

  carp(CARP_DEBUG, "Base filename is %s and path is %s", 
       filename_path_array[0], filename_path_array[1]);

  char* filename_template = get_full_filename(output_directory, 
                                              filename_path_array[0]);

  //create target file
  //  int temp_file_descriptor = -1;
  //  char suffix[25];
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  //first file is target, remaining are decoys
  //  char* psm_filename = generate_name(filename_template, "_XXXXXX", 
  //                                     ".ms2", "crux_hit_target_");

  for(file_idx=0; file_idx<total_files; file_idx++){

    char* psm_filename = generate_psm_filename(filename_path_array[0], 
                                               file_idx);

    //get file descriptor for uniq version of name //TODO remove this
    //    temp_file_descriptor = mkstemp(psm_filename);
    //open the file from the descriptor
    //    file_handle_array[file_idx] = fdopen(temp_file_descriptor, "w");
    file_handle_array[file_idx] = create_file_in_path(psm_filename,
                                                      output_directory,
                                                      overwrite); 
    //check for error
    if( file_handle_array[file_idx] == NULL ){//||
      //temp_file_descriptor == -1 ){

      carp(CARP_FATAL, "Could not create psm file %s", psm_filename);
      exit(1);
    }
    //rename this, just for a quick fix
    filename_template = get_full_filename(output_directory, psm_filename);
    //chmod(psm_filename, 0664);
    chmod(filename_template, 0664);

    //get next decoy name
    free(psm_filename);
    //sprintf(suffix, "crux_hit_decoy_%d_", file_idx+1);
    //psm_filename = generate_name(filename_template, "_XXXXXX",
    //                             ".ms2", suffix);
  }
  return file_handle_array;

}

/**
 * \brief Serialize the PSM features to output file up to 'top_hit'
 * number of top peptides from the hit_collection.
 *
 * \details  First serialize the spectrum info of the hit collection
 * then  iterate over hites and serialize the structs
 *
 * <int: charge state of the spectrum>
 * <int: Total hit objects in the hit_collection>
 * <float: delta_cn>
 * <float: ln_delta_cn>
 * <float: ln_experiment_size>
 * <BOOLEAN_T: had the score type been scored?>* - for all score types
 * <MATCH: serialized hit struct>* <--serialize top_hit hit structs 
 *
 * \returns TRUE, if sucessfully serializes the PSMs, else FALSE 
 * \callgraph
 */
BOOLEAN_T serialize_psm_features(
  MATCH_COLLECTION_T* hit_collection, ///< working hit collection -in
  FILE* output,               ///< output file handle -out
  int top_hit,              ///< number of top hit to serialize -in
  SCORER_TYPE_T prelim_score, ///< the preliminary score to report -in
  SCORER_TYPE_T main_score    ///<  the main score to report -in
  )
{
  HIT_T* hit = NULL;
  
  // create hit iterator 
  // TRUE tells iterator to return hites in sorted order of main_score type
  MATCH_ITERATOR_T* hit_iterator = 
    new_hit_iterator(hit_collection, main_score, TRUE);
  
  float delta_cn =  get_hit_collection_delta_cn(hit_collection);
  float ln_delta_cn = logf(delta_cn);
  float ln_experiment_size = logf(hit_collection->experiment_size);

  // spectrum specific features
  // first, serialize the spectrum info of the hit collection  
  // the charge of the spectrum
  
  myfwrite(&(hit_collection->charge), sizeof(int), 1, output); 
  myfwrite(&(hit_collection->hit_total), sizeof(int), 1, output);
  myfwrite(&delta_cn, sizeof(float), 1, output);
  myfwrite(&ln_delta_cn, sizeof(float), 1, output);
  myfwrite(&ln_experiment_size, sizeof(float), 1, output);
  
  // serialize each boolean for scored type 
  int score_type_idx;
  for(score_type_idx=0; score_type_idx < _SCORE_TYPE_NUM; ++score_type_idx){
    myfwrite(&(hit_collection->scored_type[score_type_idx]), 
        sizeof(BOOLEAN_T), 1, output);
  }
  
  // second, iterate over hites and serialize them
  int hit_count = 0;
  while(hit_iterator_has_next(hit_iterator)){
    ++hit_count;
    hit = hit_iterator_next(hit_iterator);        
    
    // FIXME
    prelim_score = prelim_score;
    
    // serialize hites
    serialize_hit(hit, output); // FIXME main, preliminary type
    
    // print only up to max_rank_result of the hites
    if(hit_count >= top_hit){
      break;
    }
  }
  
  free_hit_iterator(hit_iterator);
  
  return TRUE;
}

void print_sqt_header(FILE* output, char* type, int num_proteins){
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

  char* database = get_string_parameter("protein input");
  //  fprintf(output, "H\tDatabase\t%s\n", database);

  if( get_boolean_parameter("use-index") == TRUE ){
    char* fasta_name  = get_index_binary_fasta_name(database);
    free(database);
    database = fasta_name;
  }
  fprintf(output, "H\tDatabase\t%s\n", database);

  if(decoy){
  fprintf(output, "H\tComment\tDatabase shuffled; these are decoy hites\n");
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

  double tol = get_double_parameter("mass-window");
  fprintf(output, "H\tAlg-PreMasTol\t%.1f\n",tol);
  fprintf(output, "H\tAlg-FragMassTol\t%.2f\n", 
          get_double_parameter("ion-tolerance"));
  fprintf(output, "H\tAlg-XCorrMode\t0\n");

  SCORER_TYPE_T score = get_scorer_type_parameter("prelim-score-type");
  scorer_type_to_string(score, temp_str);
  fprintf(output, "H\tComment\tpreliminary algorithm %s\n", temp_str);

  score = get_scorer_type_parameter("score-type");
  scorer_type_to_string(score, temp_str);
  fprintf(output, "H\tComment\tfinal algorithm %s\n", temp_str);

  /*
  int alphabet_size = get_alphabet_size(PROTEIN_ALPH);
  char* alphabet = get_alphabet(FALSE);
  */

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
  //for letters in alphabet
  //  double mod = get_double_parameter(letter);
  //  if mod != 0
  //     double mass = mod + getmass(letter);
  //     fprintf(output, "H\tStaticMod\t%s=%.3f\n", letter, mass);
  //  fprintf(output, "H\tStaticMod\tC=160.139\n");
  fprintf(output, "H\tAlg-DisplayTop\t%d\n", 
          get_int_parameter("max-sqt-result")); 
  // this is not correct for an sqt from analzyed hites

  PEPTIDE_TYPE_T cleavages = get_peptide_type_parameter("cleavages");
  peptide_type_to_string(cleavages, temp_str);
  fprintf(output, "H\tEnzymeSpec\t%s\n", temp_str);
  //  fprintf(output, "H\t\n");
}

/**
 * \brief Print the psm features to file in sqt format.
 *
 * Prints one S line, 'top_hit' M lines, and one locus line for each
 * peptide source of each M line.
 * Assumes one spectrum per hit collection.  Could get top_hit,
 * score types from parameter.c.  Could get spectrum from first hit.
 *\returns TRUE, if sucessfully print sqt format of the PSMs, else FALSE 
 */
BOOLEAN_T print_hit_collection_sqt(
  FILE* output,                  ///< the output file -out
  int top_hit,                 ///< the top hites to output -in
  MATCH_COLLECTION_T* hit_collection,
  ///< the hit_collection to print sqt -in
  SPECTRUM_T* spectrum,          ///< the spectrum to print sqt -in
  SCORER_TYPE_T prelim_score,    ///< the preliminary score to report -in
  SCORER_TYPE_T main_score       ///< the main score to report -in
  )
{

  if( output == NULL ){
    return FALSE;
  }
  time_t hold_time;
  hold_time = time(0);
  int charge = hit_collection->charge; 
  int num_hites = hit_collection->experiment_size;

  // First, print spectrum info
  print_spectrum_sqt(spectrum, output, num_hites, charge);
  
  HIT_T* hit = NULL;
  
  // create hit iterator
  // TRUE: return hit in sorted order of main_score type
  MATCH_ITERATOR_T* hit_iterator = 
    new_hit_iterator(hit_collection, main_score, TRUE);
  
  // Second, iterate over hites, prints M and L lines
  int hit_count = 0;
  while(hit_iterator_has_next(hit_iterator)){
    ++hit_count;
    hit = hit_iterator_next(hit_iterator);    

    print_hit_sqt(hit, output, main_score, prelim_score);

    // print only up to max_rank_result of the hites
    if(hit_count >= top_hit){
      break;
    }
  }// next hit
  
  free_hit_iterator(hit_iterator);
  
  return TRUE;
}

/**
 * hit_iterator routines!
 *
 */

/**
 * create a new memory allocated hit iterator, which iterates over
 * hit iterator only one iterator is allowed to be instantiated per
 * hit collection at a time 
 *\returns a new memory allocated hit iterator
 */
MATCH_ITERATOR_T* new_hit_iterator(
  MATCH_COLLECTION_T* hit_collection,
  ///< the hit collection to iterate -out
  SCORER_TYPE_T score_type,
  ///< the score type to iterate (LOGP_EXP_SP, XCORR) -in
  BOOLEAN_T sort_hit  ///< should I return the hit in sorted order?
  )
{
  // TODO (BF 06-Feb-08): Could we pass back an iterator with has_next==False
  if (hit_collection == NULL){
    die("Null hit collection passed to hit iterator");
  }
  // is there any existing iterators?
  if(hit_collection->iterator_lock){
    carp(CARP_FATAL, 
         "Can only have one hit iterator instantiated at a time");
    exit(1);
  }
  
  // has the score type been populated in hit collection?
  if(!hit_collection->scored_type[score_type]){
    carp(CARP_FATAL, 
         "The hit collection has not been scored for request score type.");
    exit(1);
  }
  
  // allocate a new hit iterator
  MATCH_ITERATOR_T* hit_iterator = 
    (MATCH_ITERATOR_T*)mycalloc(1, sizeof(MATCH_ITERATOR_T));
  
  // set items
  hit_iterator->hit_collection = hit_collection;
  hit_iterator->hit_mode = score_type;
  hit_iterator->hit_idx = 0;
  hit_iterator->hit_total = hit_collection->hit_total;

  // only sort if requested and hit collection is not already sorted
  if(sort_hit && (hit_collection->last_sorted != score_type 
  /*|| (hit_collection->last_sorted == SP && score_type == LOGP_EXP_SP)*/)){

    if((score_type == LOGP_EXP_SP || score_type == LOGP_BONF_EXP_SP ||
        score_type == LOGP_WEIBULL_SP || score_type == LOGP_BONF_WEIBULL_SP)
       &&
       hit_collection->last_sorted == SP){
      // No need to sort, since the score_type has same rank as SP      
    }
    
    else if((score_type == LOGP_EVD_XCORR || score_type ==LOGP_BONF_EVD_XCORR)
            && hit_collection->last_sorted == XCORR){
      // No need to sort, since the score_type has same rank as XCORR
    }
    else if((score_type == Q_VALUE) &&
            hit_collection->last_sorted == PERCOLATOR_SCORE){
      // No need to sort, the score_type has same rank as PERCOLATOR_SCORE
    }
    // sort hit collection by score type
    else if(!sort_hit_collection(hit_collection, score_type)){
      carp(CARP_FATAL, "failed to sort hit collection");
      free_hit_collection(hit_collection);
      free(hit_iterator);
      exit(1);
    }
  }

  // ok lock up hit collection
  hit_collection->iterator_lock = TRUE;
  
  return hit_iterator;
}

/**
 * \brief Create a hit iterator to return hites from a collection
 * grouped by spectrum and sorted by given score type.
 *
 * \returns A heap-allocated hit iterator.
 */
MATCH_ITERATOR_T* new_hit_iterator_spectrum_sorted(
  MATCH_COLLECTION_T* hit_collection,  ///< for iteration -in
  SCORER_TYPE_T scorer ///< the score type to sort by -in
){

  MATCH_ITERATOR_T* hit_iterator = 
    (MATCH_ITERATOR_T*)mycalloc(1, sizeof(MATCH_ITERATOR_T));

  // set up fields
  hit_iterator->hit_collection = hit_collection;
  hit_iterator->hit_mode = scorer;
  hit_iterator->hit_idx = 0;
  hit_iterator->hit_total = hit_collection->hit_total;

  spectrum_sort_hit_collection(hit_collection, scorer);

  hit_collection->iterator_lock = TRUE;

  return hit_iterator;
}

/**
 * Does the hit_iterator have another hit struct to return?
 *\returns TRUE, if hit iter has a next hit, else False
 */
BOOLEAN_T hit_iterator_has_next(
  MATCH_ITERATOR_T* hit_iterator ///< the working  hit iterator -in
  )
{
  return (hit_iterator->hit_idx < hit_iterator->hit_total);
}

/**
 * return the next hit struct!
 *\returns the hit in decreasing score order for the hit_mode(SCORER_TYPE_T)
 */
HIT_T* hit_iterator_next(
  MATCH_ITERATOR_T* hit_iterator ///< the working hit iterator -in
  )
{
  return hit_iterator->hit_collection->hit[hit_iterator->hit_idx++];
}

/**
 * free the memory allocated iterator
 */
void free_hit_iterator(
  MATCH_ITERATOR_T* hit_iterator ///< the hit iterator to free
  )
{
  // free iterator
  if (hit_iterator != NULL){
    if (hit_iterator->hit_collection != NULL){
      hit_iterator->hit_collection->iterator_lock = FALSE;
    }
    free(hit_iterator);
  }
}

/*
 * Copied from spectrum_collection::serialize_header
 * uses values from paramter.c rather than taking as arguments
 */
void serialize_headers(FILE** psm_file_array){

  if( *psm_file_array == NULL ){
    return;
  }
  int num_spectrum_features = 0; //obsolete?
  int num_charged_spectra = 0;  //this is set later
  int hites_per_spectrum = get_int_parameter("top-hit");
  char* filename = get_string_parameter_pointer("protein input");
  char* protein_file = parse_filename(filename);
  //  free(filename);
  filename = get_string_parameter_pointer("ms2 file");
  char* ms2_file = parse_filename(filename);
  //  free(filename);
           

  //write values to files
  int total_files = 1 + get_int_parameter("number-decoy-set");
  int i=0;
  for(i=0; i<total_files; i++){
    fwrite(&(num_charged_spectra), sizeof(int), 1, psm_file_array[i]);
    fwrite(&(num_spectrum_features), sizeof(int), 1, psm_file_array[i]);
    fwrite(&(hites_per_spectrum), sizeof(int), 1, psm_file_array[i]);
  }
  
  free(protein_file);
  free(ms2_file);

}

/**
 * \brief Writes the contents of a hit_collection to file(s)
 *
 * \details Takes information from parameter.c to decide which files
 * (binary, sqt) to write to, how many hites to write, etc.
 *
 */
void print_hites( 
                   MATCH_COLLECTION_T* hit_collection, ///< results to write
                   SPECTRUM_T* spectrum, ///< results for this spectrum
                   BOOLEAN_T is_decoy,   ///< peptides from target/decoy
                   FILE* psm_file,       ///< binary file -out
                   FILE* sqt_file,       ///< text file, target -out
                   FILE* decoy_file){    ///< text file, decoy -out

  carp(CARP_DETAILED_DEBUG, "Writing hites to file");
  // get parameters
  MATCH_SEARCH_OUTPUT_MODE_T output_type = get_output_type_parameter(
                                                            "output-mode");
  int max_sqt_hites = get_int_parameter("max-sqt-result");
  int max_psm_hites = get_int_parameter("top-hit");
  SCORER_TYPE_T main_score = get_scorer_type_parameter("score-type");
  SCORER_TYPE_T prelim_score = get_scorer_type_parameter("prelim-score-type");


  // write binary files
  if( output_type != SQT_OUTPUT ){ //i.e. binary or all
    carp(CARP_DETAILED_DEBUG, "Serializing psms");
    serialize_psm_features(hit_collection, psm_file, max_psm_hites,
                           prelim_score, main_score);
  }

  // write sqt files
  if( output_type != BINARY_OUTPUT ){ //i.e. sqt or all
    carp(CARP_DETAILED_DEBUG, "Writing sqt results");
    if( ! is_decoy ){
      print_hit_collection_sqt(sqt_file, max_sqt_hites,
                                 hit_collection, spectrum,
                                 prelim_score, main_score);
    }else{
      print_hit_collection_sqt(decoy_file, max_sqt_hites,
                                 hit_collection, spectrum,
                                 prelim_score, main_score);
    }
  }
}

/*******************************************
 * hit_collection post_process extension
 ******************************************/

/**
 * \brief Creates a new hit_collection from the PSM iterator.
 *
 * Used in the post_processing extension.  Also used by
 * setup_hit_collection_iterator which is called by next to find,
 * open, and parse the next psm file(s) to process.  If there are
 * multiple target psm files, it reads in all of them when set_type is
 * 0 and puts them all into one hit_collection. 
 *\returns A heap allocated hit_collection.
 */
MATCH_COLLECTION_T* new_hit_collection_psm_output(
  MATCH_COLLECTION_ITERATOR_T* hit_collection_iterator, 
    ///< the working hit_collection_iterator -in
  SET_TYPE_T set_type  
    ///< what set of hit collection are we creating? (TARGET, DECOY1~3) -in 
  )
{ 
  struct dirent* directory_entry = NULL;
  char* file_in_dir = NULL;
  FILE* result_file = NULL;
  char suffix[25];
  //  char prefix[25];

  carp(CARP_DEBUG, "Calling new_hit_collection_psm_output");
  DATABASE_T* database = hit_collection_iterator->database;
  
  // allocate hit_collection object
  MATCH_COLLECTION_T* hit_collection = allocate_hit_collection();

  // set this as a post_process hit collection
  hit_collection->post_process_collection = TRUE;
  
  // the protein counter size, create protein counter
  hit_collection->post_protein_counter_size 
   = get_database_num_proteins(database);
  hit_collection->post_protein_counter 
   = (int*)mycalloc(hit_collection->post_protein_counter_size, sizeof(int));
  hit_collection->post_protein_peptide_counter 
   = (int*)mycalloc(hit_collection->post_protein_counter_size, sizeof(int));

  // create hash table for peptides
  // Set initial capacity to protein count.
  hit_collection->post_hash 
    = new_hash(hit_collection->post_protein_counter_size);
  
  // set the prefix of the serialized file to parse
  // Also, tag if hit_collection type is null_peptide_collection

  /*char* filepath = hit_collection_iterator->directory_name;
  char** file_path_array = parse_filename_path_extension(filepath, ".csm");
  char* file_prefix = file_path_array[0];
  char* directory_name = file_path_array[1];
  if( directory_name == NULL ){
    directory_name = ".";
    }*/

  if(set_type == TARGET){
    //    sprintf(suffix, "crux_hit_target");
    //sprintf(prefix, "crux_hit_target");
    sprintf(suffix, ".csm");
    hit_collection->null_peptide_collection = FALSE;
  }
  else{
    //    sprintf(suffix, "crux_hit_decoy_%d", (int)set_type);
    //sprintf(prefix, "crux_hit_decoy_%d", (int)set_type);
    sprintf(suffix, "-decoy-%d.csm", (int)set_type);
    hit_collection->null_peptide_collection = TRUE;
  }
  
  carp(CARP_DEBUG, "Set type is %d and suffix is %s", (int)set_type, suffix);
  // iterate over all PSM files in directory to find the one to read
  while((directory_entry 
            = readdir(hit_collection_iterator->working_directory))){

    //carp(CARP_DETAILED_DEBUG, "Next file is %s", directory_entry->d_name);

    if (strcmp(directory_entry->d_name, ".") == 0 ||
        strcmp(directory_entry->d_name, "..") == 0 ||
        !suffix_compare(directory_entry->d_name, suffix)
        //!prefix_compare(directory_entry->d_name, prefix)
        ) {
      continue;
    }

    if( set_type == TARGET && name_is_decoy(directory_entry->d_name) ){
      continue;
    }
    file_in_dir =get_full_filename(hit_collection_iterator->directory_name, 
                                   directory_entry->d_name);

    carp(CARP_INFO, "Getting PSMs from %s", file_in_dir);
    result_file = fopen(file_in_dir, "r");
    if( access(file_in_dir, R_OK)){
      carp(CARP_FATAL, "Cannot read from psm file '%s'", file_in_dir);
      exit(1);
    }
    // add all the hit objects from result_file
    extend_hit_collection(hit_collection, database, result_file);
    carp(CARP_DETAILED_DEBUG, "Extended hit collection " );
    fclose(result_file);
    free(file_in_dir);
    carp(CARP_DETAILED_DEBUG, "Finished file.");
  }
  
  return hit_collection;
}


/**
 * parse all the hit objects and add to hit collection
 *\returns TRUE, if successfully parse all PSMs in result_file, else FALSE
 */
BOOLEAN_T extend_hit_collection(
  MATCH_COLLECTION_T* hit_collection, ///< hit collection to extend -out
  DATABASE_T* database, ///< the database holding the peptides -in
  FILE* result_file   ///< the result file to parse PSMs -in
  )
{
  int total_spectra = 0;
  int hit_idx = 0;
  int spectrum_idx = 0;
  int charge = 0;
  HIT_T* hit = NULL;
  int num_top_hit = 0;
  int num_spectrum_features = 0;
  float delta_cn =  0;
  float ln_delta_cn = 0;
  float ln_experiment_size = 0;
  int hit_total_of_serialized_collection = 0;
  int score_type_idx = 0;
  BOOLEAN_T type_scored = FALSE;

  // only for post_process_collections
  if(!hit_collection->post_process_collection){
    carp(CARP_ERROR, "Must be a post process hit collection");
    return FALSE;
  }
  
  // read in file specific info
  
  // get number of spectra serialized in the file
  if(fread(&total_spectra, (sizeof(int)), 1, result_file) != 1){
    carp(CARP_ERROR,"Serialized file corrupted, incorrect number of spectra");
    return FALSE;
  }
  carp(CARP_DETAILED_DEBUG, "There are %i spectra in the result file", 
       total_spectra);

  // FIXME unused feature, just set to 0
  // get number of spectra features serialized in the file
  if(fread(&num_spectrum_features, (sizeof(int)), 1, result_file) != 1){
    carp(CARP_ERROR, 
         "Serialized file corrupted, incorrect number of spectrum features");
    return FALSE;
  }
  
  carp(CARP_DETAILED_DEBUG, "There are %i spectrum features", 
       num_spectrum_features);

  // get number top ranked peptides serialized
  if(fread(&num_top_hit, (sizeof(int)), 1, result_file) != 1){
    carp(CARP_ERROR, 
         "Serialized file corrupted, incorrect number of top hit");  
    return FALSE;
  }
  carp(CARP_DETAILED_DEBUG, "There are %i top hites", num_top_hit);


  // FIXME
  // could parse fasta file and ms2 file
  
  
  // now iterate over all spectra serialized
  for(spectrum_idx = 0; spectrum_idx < total_spectra; ++spectrum_idx){
    /*** get all spectrum specific features ****/
    
    // get charge of the spectrum
    /*    if(fread(&charge, (sizeof(int)), 1, result_file) != 1){
      carp(CARP_ERROR, "Serialized file corrupted, charge value not read");  
      return FALSE;
    }
    */
    int chars_read = fread(&charge, (sizeof(int)), 1, result_file);
    carp(CARP_DETAILED_DEBUG, "Read %i characters, charge is %i",
         chars_read, charge);

    // get serialized hit_total
    /*if(fread(&hit_total_of_serialized_collection, (sizeof(int)),
             1, result_file) != 1){
      carp(CARP_ERROR, "Serialized file corrupted, "
          "incorrect hit_total_of_serialized_collection value");  
      return FALSE;
      }*/
    chars_read = fread(&hit_total_of_serialized_collection, (sizeof(int)),
                       1, result_file);
    carp(CARP_DETAILED_DEBUG, "Read %i characters, hit total is %i",
         chars_read, hit_total_of_serialized_collection);
      
    // get delta_cn value
    if(fread(&delta_cn, (sizeof(float)), 1, result_file) != 1){
      carp(CARP_ERROR, 
       "Serialized file corrupted, incorrect delta cn value for top hit");  
      return FALSE;
    }
    
    // get ln_delta_cn value
    if(fread(&ln_delta_cn, (sizeof(float)), 1, result_file) != 1){
      carp(CARP_ERROR, 
    "Serialized file corrupted, incorrect ln_delta cn value for top hit");  
      return FALSE;
    }
    
    // get ln_experiment_size
    if(fread(&ln_experiment_size, (sizeof(float)), 1, result_file) != 1){
      carp(CARP_ERROR, "Serialized file corrupted, incorrect "
           "ln_experiment_size cn value for top hit");  
      return FALSE;
    }
    
    // Read each boolean for scored type 
    // parse all boolean indicators for scored hit object
    for(score_type_idx=0; score_type_idx < _SCORE_TYPE_NUM; ++score_type_idx){
      fread(&(type_scored), sizeof(BOOLEAN_T), 1, result_file);
      
      // if this is the first time extending the hit collection
      // set scored boolean values
      if(!hit_collection->post_scored_type_set){
        hit_collection->scored_type[score_type_idx] = type_scored;
      }
      else{
        // if boolean values already set compare if no
        // conflicting scored types 
        if(hit_collection->scored_type[score_type_idx] != type_scored){
          carp(CARP_ERROR, "Serialized hit objects has not been scored "
               "as other hit objects");
        }
      }
      
      // now once we are done with setting scored type
      // set hit collection status as set!
      if(!hit_collection->post_scored_type_set &&
         score_type_idx == (_SCORE_TYPE_NUM-1)){
        hit_collection->post_scored_type_set = TRUE;
      }
    }
    
    // now iterate over all 
    for(hit_idx = 0; hit_idx < num_top_hit; ++hit_idx){
      // break if there are no hit objects serialized
      //      if(hit_total_of_serialized_collection <= 0){
      if(hit_total_of_serialized_collection <= hit_idx){
        break;
      }
      
      carp(CARP_DETAILED_DEBUG, "Reading hit %i", hit_idx);
      // parse hit object
      if((hit = parse_hit(result_file, database))==NULL){
        carp(CARP_ERROR, "Failed to parse serialized PSM hit");
        return FALSE;
      }
      
      // set all spectrum specific features to parsed hit
      set_hit_charge(hit, charge);
      set_hit_delta_cn(hit, delta_cn);
      set_hit_ln_delta_cn(hit, ln_delta_cn);
      set_hit_ln_experiment_size(hit, ln_experiment_size);
      
      // now add hit to hit collection
      add_hit_to_hit_collection(hit_collection, hit);
    }// next hit for this spectrum

  }// next spectrum
  
  return TRUE;
}

/**
 * Adds the hit object to hit_collection
 * Must not exceed the _MAX_NUMBER_PEPTIDES to be hit added
 * \returns TRUE if successfully adds the hit to the
 * hit_collection, else FALSE 
 */
BOOLEAN_T add_hit_to_hit_collection(
  MATCH_COLLECTION_T* hit_collection, ///< the hit collection to free -out
  HIT_T* hit ///< the hit to add -in
  )
{
  PEPTIDE_T* peptide = NULL;

  // only for post_process_collections
  if(!hit_collection->post_process_collection){
    carp(CARP_ERROR, "Must be a post process hit collection");
    return FALSE;
  }
  
  // check if enough space for peptide hit
  if(hit_collection->hit_total >= _MAX_NUMBER_PEPTIDES){
    carp(CARP_ERROR, "Rich hit count exceeds max hit limit: %d", 
         _MAX_NUMBER_PEPTIDES);
    return FALSE;
  }
  
  // add a new hit to array
  hit_collection->hit[hit_collection->hit_total] = hit;
  
  // increment total rich hit count
  ++hit_collection->hit_total;
  
  // DEBUG, print total peptided scored so far
  if(hit_collection->hit_total % 1000 == 0){
    carp(CARP_INFO, "parsed PSM: %d", hit_collection->hit_total);
  }
  
  // hit peptide
  peptide = get_hit_peptide(hit);
  
  // update protein counter, protein_peptide counter
  update_protein_counters(hit_collection, peptide);
  
  // update hash table
  char* hash_value = get_peptide_hash_value(peptide); 
  add_hash(hit_collection->post_hash, hash_value, NULL); 
  
  return TRUE;
}

/**
 * updates the protein_counter and protein_peptide_counter for 
 * run specific features
 */
void update_protein_counters(
  MATCH_COLLECTION_T* hit_collection, ///< working hit collection -in
  PEPTIDE_T* peptide  ///< peptide information to update counters -in
  )
{
  PEPTIDE_SRC_ITERATOR_T* src_iterator = NULL;
  PEPTIDE_SRC_T* peptide_src = NULL;
  PROTEIN_T* protein = NULL;
  unsigned int protein_idx = 0;
  int hash_count = 0;
  BOOLEAN_T unique = FALSE;
  
  // only for post_process_collections
  if(!hit_collection->post_process_collection){
    carp(CARP_ERROR, "Must be a post process hit collection");
    exit(1);
  }
  
  // See if this peptide has been observed before?
  char* hash_value = get_peptide_hash_value(peptide);
  hash_count = get_hash_count(hit_collection->post_hash, hash_value);
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
    protein_idx = get_protein_protein_idx(protein);
    
    // update the number of PSM this protein hites
    ++hit_collection->post_protein_counter[protein_idx];
    
    // number of peptides hit this protein
    if(unique){
      ++hit_collection->post_protein_peptide_counter[protein_idx];
    }
  }  
  
  free_peptide_src_iterator(src_iterator);
}

/**
 * Fill the hit objects score with the given the float array. 
 * The hit object order must not have been altered since scoring.
 * The result array size must hit the hit_total count.
 * Match ranks are also populated to preserve the original order of the
 * hit input TRUE for preserve_order.
 *\returns TRUE, if successfully fills the scores into hit object, else FALSE.
 */
BOOLEAN_T fill_result_to_hit_collection(
  MATCH_COLLECTION_T* hit_collection, 
    ///< the hit collection to iterate -out
  double* results,  
    ///< The result score array to fill the hit objects -in
  SCORER_TYPE_T score_type,  
    ///< The score type of the results to fill (XCORR, Q_VALUE, ...) -in
  BOOLEAN_T preserve_order 
    ///< preserve hit order?
  )
{
  int hit_idx = 0;
  HIT_T* hit = NULL;
  HIT_T** hit_array = NULL;
  SCORER_TYPE_T score_type_old = hit_collection->last_sorted;

  // iterate over hit object in collection, set scores
  for(; hit_idx < hit_collection->hit_total; ++hit_idx){
    hit = hit_collection->hit[hit_idx];
    set_hit_score(hit, score_type, results[hit_idx]);    
  }
  
  // if need to preserve order store a copy of array in original order 
  if(preserve_order){
    hit_array = (HIT_T**)mycalloc(hit_collection->hit_total, sizeof(HIT_T*));
    for(hit_idx=0; hit_idx < hit_collection->hit_total; ++hit_idx){
      hit_array[hit_idx] = hit_collection->hit[hit_idx];
    }
  }

  // populate the rank of hit_collection
  if(!populate_hit_rank_hit_collection(hit_collection, score_type)){
    carp(CARP_ERROR, "failed to populate hit rank in hit_collection");
    free_hit_collection(hit_collection);
    exit(1);
  }
  
  // restore hit order if needed
  if(preserve_order){
    for(hit_idx=0; hit_idx < hit_collection->hit_total; ++hit_idx){
      hit_collection->hit[hit_idx] = hit_array[hit_idx];
    }
    hit_collection->last_sorted = score_type_old;
    free(hit_array);
  }

  hit_collection->scored_type[score_type] = TRUE;
  
  return TRUE;
}

/**
 * Process run specific features from all the PSMs
 */
void process_run_specific_features(
  MATCH_COLLECTION_T* hit_collection ///< the hit collection to free -out
  );


/**********************************
 * hit_collection get, set methods
 **********************************/

/**
 *\returns the hit_collection protein counter for the protein idx
 */
int get_hit_collection_protein_counter(
  MATCH_COLLECTION_T* hit_collection, ///< the working hit collection -in
  unsigned int protein_idx ///< the protein index to return protein counter -in
  )
{
  // only for post_process_collections
  if(!hit_collection->post_process_collection){
    carp(CARP_ERROR, "Must be a post process hit collection");
    exit(1);
  }

  // number of PSMs hit this protein
  return hit_collection->post_protein_counter[protein_idx];
}

/**
 *\returns the hit_collection protein peptide counter for the protein idx
 */
int get_hit_collection_protein_peptide_counter(
  MATCH_COLLECTION_T* hit_collection, ///< the working hit collection -in
  unsigned int protein_idx ///< the protein index to return protein peptiide counter -in
  )
{
  // only for post_process_collections
  if(!hit_collection->post_process_collection){
    carp(CARP_ERROR, "Must be a post process hit collection");
    exit(1);
  }
  
  // number of peptides hit this protein
  return hit_collection->post_protein_peptide_counter[protein_idx];
}

/**
 *\returns the hit_collection hash value of PSMS for which this is the best scoring peptide
 */
int get_hit_collection_hash(
  MATCH_COLLECTION_T* hit_collection, ///< the working hit collection -in
  PEPTIDE_T* peptide  ///< the peptide to check hash value -in
  )
{
  // only for post_process_collections
  if(!hit_collection->post_process_collection){
    carp(CARP_ERROR, "Must be a post process hit collection");
    exit(1);
  }
  
  char* hash_value = get_peptide_hash_value(peptide);
  int count = get_hash_count(hit_collection->post_hash, hash_value);
  free(hash_value);
  
  return count;
}

/**
 * \brief Get the number of proteins in the database associated with
 * this hit collection.
 */
int get_hit_collection_num_proteins(
  MATCH_COLLECTION_T* hit_collection ///< the hit collection of interest -
  ){

  return hit_collection->post_protein_counter_size;
}


/******************************
 * hit_collection_iterator
 ******************************/
     
/**
 * \brief Finds the next hit_collection in directory and prepares
 * the iterator to hand it off when 'next' called.
 *
 * When no more hit_collections (i.e. psm files) are available, set
 * hit_collection_iterator->is_another_collection to FALSE
 * \returns void
 */
void setup_hit_collection_iterator(
  MATCH_COLLECTION_ITERATOR_T* hit_collection_iterator 
    ///< the hit_collection_iterator to set up -in/out
  )
{
  // are there any more hit_collections to return?
  if(hit_collection_iterator->collection_idx 
      < hit_collection_iterator->number_collections){

    // then go parse the hit_collection
    hit_collection_iterator->hit_collection = 
      new_hit_collection_psm_output(hit_collection_iterator, 
         (SET_TYPE_T)hit_collection_iterator->collection_idx);

    // we have another hit_collection to return
    hit_collection_iterator->is_another_collection = TRUE;
    
    // let's move on to the next one next time
    ++hit_collection_iterator->collection_idx;

    // reset directory
    rewinddir(hit_collection_iterator->working_directory);
  }
  else{
    // we're done, no more hit_collections to return
    hit_collection_iterator->is_another_collection = FALSE;
  }
}

/**
 * Create a hit_collection iterator from a directory of serialized files
 * Only hadles up to one target and three decoy sets per folder
 *\returns hit_collection iterator instantiated from a result folder
 */
MATCH_COLLECTION_ITERATOR_T* new_hit_collection_iterator(
  char* output_file_directory, 
    ///< the directory path where the PSM output files are located -in
  char* fasta_file 
    ///< The name of the fasta file for peptides for hit_collections. -in
  )
{
  carp(CARP_DEBUG, 
       "Creating hit collection iterator for dir %s and protein input %s",
       output_file_directory, fasta_file);

  // allocate hit_collection
  MATCH_COLLECTION_ITERATOR_T* hit_collection_iterator =
    (MATCH_COLLECTION_ITERATOR_T*)
      mycalloc(1, sizeof(MATCH_COLLECTION_ITERATOR_T));

  DIR* working_directory = NULL;
  struct dirent* directory_entry = NULL;
  DATABASE_T* database = NULL;
  BOOLEAN_T use_index_boolean = get_boolean_parameter("use-index");  
  if( use_index_boolean == FALSE ){
    carp(CARP_ERROR, 
         "Analyzing hites without index not implemented.  Using index.");
    use_index_boolean = TRUE;
  }

  // get directory from path name and prefix from filename
  /*
  char** file_path_array = parse_filename_path_extension(  // free me
                                     output_file_directory, ".csm");
  char* filename = parse_filename(output_file_directory);  // free me
  char* filename_prefix = file_path_array[0];
  char* decoy_prefix = cat_string(filename_prefix, "-decoy-");// free me
  char* psm_dir_name = file_path_array[1];
  if( psm_dir_name == NULL ){
    psm_dir_name = ".";
  }
  */

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
  //working_directory = opendir(psm_dir_name);
  
  if(working_directory == NULL){
    carp(CARP_FATAL, "Failed to open PSM file directory: %s", 
        output_file_directory);
    exit(1);
  }
  
  // determine how many decoy sets we have
  while((directory_entry = readdir(working_directory))){
    //if(prefix_compare(directory_entry->d_name, "crux_hit_target")){
    /*if(suffix_compare(directory_entry->d_name, ".csm")){
      carp(CARP_DEBUG, "Found target file %s", directory_entry->d_name);
      boolean_result = TRUE;
    }
    //else if(prefix_compare(directory_entry->d_name, "crux_hit_decoy_1")) {
    else*/ if(suffix_compare(directory_entry->d_name, "decoy-1.csm")) {
      carp(CARP_DEBUG, "Found decoy file %s", directory_entry->d_name);
      decoy_1 = TRUE;
      //num_decoys++;
    }
    //else if(prefix_compare(directory_entry->d_name, "crux_hit_decoy_2")) {
    else if(suffix_compare(directory_entry->d_name, "decoy-2.csm")) {
      decoy_2 = TRUE;
    }
    //else if(prefix_compare(directory_entry->d_name, "crux_hit_decoy_3")) {
    else if(suffix_compare(directory_entry->d_name, "decoy-3.csm")) {
      decoy_3 = TRUE;
      break;
    }    
    else if(suffix_compare(directory_entry->d_name, ".csm")){
      carp(CARP_DEBUG, "Found target file %s", directory_entry->d_name);
      boolean_result = TRUE;
    }
  }
  
  // set total_sets count
  int total_sets = 0;

  if(decoy_3){
    total_sets = 4; // 3 decoys + 1 target
  }
  else if(decoy_2){
    total_sets = 3; // 2 decoys + 1 target
  }
  else if(decoy_1){
    total_sets = 2; // 1 decoys + 1 target
  }
  else{
    total_sets = 1;
    carp(CARP_INFO, "No decoy sets exist in directory: %s", 
        output_file_directory);
  }
  if(!boolean_result){
    carp(CARP_FATAL, "No PSM files found in directory '%s'", 
         output_file_directory);
    exit(1);
  }

  // get binary fasta file name with path to crux directory 
  //  char* binary_fasta = get_binary_fasta_name_in_crux_dir(fasta_file);
  char* binary_fasta = get_index_binary_fasta_name(fasta_file);
  
  // check if input file exist
  if(access(binary_fasta, F_OK)){
    carp(CARP_FATAL, "The file \"%s\" does not exist (or is not readable, "
        "or is empty) for crux index.", binary_fasta);
    free(binary_fasta);
    exit(1);
  }
  
  carp(CARP_DEBUG, "Creating a new database");
  // now create a database, 
  // using fasta file either binary_file(index) or fastafile
  database = new_database(binary_fasta, use_index_boolean);
  
  // check if already parsed
  if(!get_database_is_parsed(database)){
    carp(CARP_DETAILED_DEBUG,"Parsing database");
    if(!parse_database(database)){
      carp(CARP_FATAL, "Failed to parse database, cannot create new index");
      free_database(database);
      exit(1);
    }
  }
  
  free(binary_fasta);

  // reset directory
  rewinddir(working_directory);
  
  // set hit_collection_iterator fields
  hit_collection_iterator->working_directory = working_directory;
  hit_collection_iterator->database = database;  
  hit_collection_iterator->number_collections = total_sets;
  hit_collection_iterator->directory_name = 
    my_copy_string(output_file_directory);
  hit_collection_iterator->is_another_collection = FALSE;

  // setup the hit collection iterator for iteration
  // here it will go parse files to construct hit collections
  setup_hit_collection_iterator(hit_collection_iterator);

  // clean up strings
  //free(file_path_array);
  //free(filename);
  //free(decoy_prefix);

  return hit_collection_iterator;
}

/**
 *\returns TRUE, if there's another hit_collection to return, else return FALSE
 */
BOOLEAN_T hit_collection_iterator_has_next(
  MATCH_COLLECTION_ITERATOR_T* hit_collection_iterator ///< the working hit_collection_iterator -in
  )
{
  // Do we have another hit_collection to return
  return hit_collection_iterator->is_another_collection;
}

/**
 * free hit_collection_iterator
 */
void free_hit_collection_iterator(
  MATCH_COLLECTION_ITERATOR_T* hit_collection_iterator ///< the working hit_collection_iterator -in
  )
{
  // free unclaimed hit_collection
  if(hit_collection_iterator->hit_collection != NULL){
    free_hit_collection(hit_collection_iterator->hit_collection);
  }
  
  // free up all hit_collection_iterator 
  free(hit_collection_iterator->directory_name);
  free_database(hit_collection_iterator->database);
  closedir(hit_collection_iterator->working_directory); 
  free(hit_collection_iterator);
}

/**
 * \brief Fetches the next hit collection object and prepares for
 * the next iteration 
 *\returns The next hit collection object
 */
MATCH_COLLECTION_T* hit_collection_iterator_next(
  MATCH_COLLECTION_ITERATOR_T* hit_collection_iterator 
    ///< the working hit_collection_iterator -in
  )
{
  MATCH_COLLECTION_T* hit_collection = NULL;
  
  if(hit_collection_iterator->is_another_collection){
    hit_collection = hit_collection_iterator->hit_collection;
    hit_collection_iterator->hit_collection = NULL;
    setup_hit_collection_iterator(hit_collection_iterator);
    return hit_collection;
  }
  else{
    carp(CARP_ERROR, "No hit_collection to return");
    return NULL;
  }
}

/**
 *\returns the total number of hit_collections to return
 */
int get_hit_collection_iterator_number_collections(
  MATCH_COLLECTION_ITERATOR_T* hit_collection_iterator ///< the working hit_collection_iterator -in
  )
{
  return hit_collection_iterator->number_collections;
}

/**
 * \brief Get the name of the directory the hit_collection_iterator
 * is working in.
 * \returns A heap allocated string (char*) of the directory name.
 */
char* get_hit_collection_iterator_directory_name(
  MATCH_COLLECTION_ITERATOR_T* iterator ///< the hit_collection_iterator -in
  ){

  char* dir_name = my_copy_string(iterator->directory_name);

  return dir_name;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

