/*****************************************************************************
 * \file match.c
 * AUTHOR: Chris Park
 * CREATE DATE: 11/27 2006
 * DESCRIPTION: Object for matching a peptide and a spectrum, generate a preliminary score(ex, Sp)
 *
 * REVISION: 
 ****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
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
#include "match.h" 
#include "match_collection.h" 
#include "generate_peptides_iterator.h" 
#include "peptide.h"

/**
 * README!
 * <Issues on the overall_type field in match struct>
 * 
 * Outstanding question: How do you determine the
 * peptide trypticity for multiple protein sources?
 * 
 * For example, if one protein src is tryptic while the other is 
 * not tryptic, what is the peptide trypticity used for feature and 
 * shuffling peptide sequence?
 *
 * Currently, we use the "Tryptic wins all" approach, where
 * if N-terminus is tryptic in any of the src proteins we claim it
 * tryptic on the N terminus. Same applies for C-terminus.
 * Thus, even if no peptide is tryptic in any of its src protein,
 * if there are one src protein where it is N tryptic and another where
 * it is C tryptic, overall we will call the peptide in the match tryptic.
 * 
 * There are for sure other methods, for example to randomly sample the src protein
 * and consider that protein as its src or to shuffle the flanking sequence of the peptide
 * in each of peptide and randomly sample from the shuffled flanking sequence to determine
 * the shuffled peptide's trypticity.
 *
 */

/**
 *\struct match
 *\brief An object that stores the score & rank for each match of spectrum & peptide
 */
struct match{
  SPECTRUM_T* spectrum; ///< the spectrum we are scoring with
  PEPTIDE_T* peptide;  ///< the peptide we are scoring
  float match_scores[_SCORE_TYPE_NUM]; ///< the scoring result array (use enum_type SCORER_TYPE_T to index)
  int match_rank[_SCORE_TYPE_NUM];  ///< the rank of scoring result (use enum_type SCORER_TYPE_T to index)
  int pointer_count; ///< the number of pointers to this match object (when reach 0, free memory)
  float b_y_ion_match; ///< the fraction of the b, y ion matched while scoring for SP
  BOOLEAN_T null_peptide; ///< Is the match a null peptide match?
  char* peptide_sequence; ///< cached peptide sequence, if not called before set as NULL
  PEPTIDE_TYPE_T overall_type; ///< the overall peptide trypticity, this is set in set_match_peptide routine, go to README top
  int charge; ///< the charge state of the match collection created
  //post_process match object features
  //only valid when post_process_match is TRUE
  BOOLEAN_T post_process_match; ///< Is this a post process match object?
  float delta_cn; ///< the difference in top and second Xcorr scores
  float ln_delta_cn; ///< the natural log of delta_cn
  float ln_experiment_size; ///< the natural log of total number of candidate peptides evaluated
};

/**
 * \returns a new memory allocated match
 */
MATCH_T* new_match(void){
  MATCH_T* match = (MATCH_T*)mycalloc(1, sizeof(MATCH_T));
  
  //initialize score, rank !!!!DEBUG
  int index = 0;
  for(; index < _SCORE_TYPE_NUM; ++index){
    match->match_rank[index] = 0;
    match->match_scores[index] = 0;
  }
  
  ++match->pointer_count;

  //default is not a null peptide match
  match->null_peptide = FALSE;

  //set default as not tryptic
  //a full evaluation is done when set peptide
  match->overall_type = NOT_TRYPTIC;

  return match;
}

/**
 * free the memory allocated match
 * spectrum is not freed by match
 */
void free_match(
  MATCH_T* match ///< the match to free -in
  )
{
  --match->pointer_count;
  
  //only free match when pointer count reaches
  if(match->pointer_count == 0){
    //free peptide
    free_peptide(match->peptide);
    
    if(match->post_process_match){
      //free spectrum
      free_spectrum(match->spectrum);
    }

    //free cached peptide sequence if not NULL
    free(match->peptide_sequence);

    free(match);  
  }
}

/**
 * compare two matches, used for qsort
 * \returns the difference between sp score in match_a and match_b
 */
int compare_match_sp(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
  )
{
  //might have to worry about cases below 1 and -1
  //return (int)((*match_b)->match_scores[SP] - (*match_a)->match_scores[SP]);


  if((*match_b)->match_scores[SP] > (*match_a)->match_scores[SP]){
    return 1;
  }
  else if((*match_b)->match_scores[SP] < (*match_a)->match_scores[SP]){
    return -1;
  }
  return 0;

}

/**
 * compare two matches, used for qsort
 * \returns the difference between xcorr score in match_a and match_b
 */
int compare_match_xcorr(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
)
{

  if((*match_b)->match_scores[XCORR] > (*match_a)->match_scores[XCORR]){
    return 1;
  }
  else if((*match_b)->match_scores[XCORR] < (*match_a)->match_scores[XCORR]){
    return -1;
  }
  return 0;

}


/**
 * compare two matches, used for qsort
 * \returns the difference between xcorr score in match_a and match_b
 */
int compare_match_q_value(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
)
{

  if((*match_b)->match_scores[Q_VALUE] > (*match_a)->match_scores[Q_VALUE]){
    return 1;
  }
  else if((*match_b)->match_scores[Q_VALUE] < (*match_a)->match_scores[Q_VALUE]){
    return -1;
  }
  return 0;
}

/**
 * print the information of the match
 */
void print_match(
  MATCH_T* match, ///< the match to print -in  
  FILE* file, ///< output stream -out
  BOOLEAN_T output_sequence, ///< should I output peptide sequence -in
  SCORER_TYPE_T output_mode  ///< the output mode -in
  )
{
  char* peptide_sequence = NULL;
  
  //print according to the output mode
  //

  SCORER_TYPE_T primary_score = output_mode;
  SCORER_TYPE_T secondary_score = SP;
  switch (output_mode) {
    case DOTP:
      //FIXME fill in once implemented
      break;
    case SP:
    case XCORR:
    case LOGP_EXP_SP:
    case LOGP_BONF_EXP_SP:
    case LOGP_WEIBULL_SP:
    case LOGP_BONF_WEIBULL_SP:
    case LOGP_EVD_XCORR:
    case LOGP_BONF_EVD_XCORR:
    case LOGP_WEIBULL_XCORR:
    case LOGP_BONF_WEIBULL_XCORR:
    case Q_VALUE:      
    case PERCOLATOR_SCORE:  
      secondary_score = SP;
      break;
  }

  fprintf(file, "P %d\t%d\t%.2f\t%.2f\t%.2f\t", 
      match->match_rank[primary_score], 
      match->match_rank[secondary_score], 
      get_peptide_peptide_mass(match->peptide), 
      match->match_scores[primary_score], 
      match->match_scores[secondary_score]);

  // TODO resolve spectrum_header output and above to not be coupled
 
  //should I print sequence?
  if(output_sequence){
    peptide_sequence = get_match_sequence(match);
    fprintf(file, "%s\n", peptide_sequence);
    free(peptide_sequence);
  }
}

/**
 * sort the match array with the corresponding compare method
 */
void qsort_match(
  MATCH_T** match_array, ///< the match array to sort -in  
  int match_total,  ///< the total number of match objects -in
  void* compare_method ///< the compare method to use -in
  )
{
  qsort(match_array, match_total, sizeof(MATCH_T*), compare_method);
}

/**
 * serializes the match in binary
 *
 */
void serialize_match(
  MATCH_T* match, ///< the match to print -in
  FILE* file ///< output stream -out
  )
{
  //first serialize peptide
  serialize_peptide(match->peptide, file);
  
  //Serialize each score and rank
  int score_type_idx = 0;
  for(; score_type_idx < _SCORE_TYPE_NUM; ++score_type_idx){
    fwrite(&(match->match_scores[score_type_idx]), sizeof(float), 1, file);
    fwrite(&(match->match_rank[score_type_idx]), sizeof(int), 1, file);
  }
  
  //serialize spectrum in binary
  serialize_spectrum(match->spectrum, file);
  
  //b/y ion matches ratio
  fwrite(&(match->b_y_ion_match), sizeof(float), 1, file);

  //serialize match peptide overall trypticity
  fwrite(&(match->overall_type), sizeof(PEPTIDE_TYPE_T), 1, file);
  
  //serialize match is it null_peptide?
  fwrite(&(match->null_peptide), sizeof(BOOLEAN_T), 1, file);
}

/*******************************************
 * match post_process extension
 ******************************************/

/**
 * Constructs the 20 feature array that pass over to percolator registration
 * Go to top README for N,C terminus tryptic feature info.
 *\returns the feature float array
 */
double* get_match_percolator_features(
  MATCH_T* match, ///< the match to work -in                                          
  MATCH_COLLECTION_T* match_collection ///< the match collection to iterate -in
  )
{
  int feature_count = 20;
  PEPTIDE_SRC_ITERATOR_T* src_iterator = NULL;
  PEPTIDE_SRC_T* peptide_src = NULL;
  PROTEIN_T* protein = NULL;
  unsigned int protein_idx = 0;
  double* feature_array = (double*)mycalloc(feature_count, sizeof(double));
  float weight_diff = get_peptide_peptide_mass(match->peptide) -
    get_spectrum_neutral_mass(match->spectrum, match->charge);

  //Xcorr
  feature_array[0] = get_match_score(match, XCORR);
  //DeltCN
  feature_array[1] = match->delta_cn;
  //DeltLCN
  feature_array[2] = match->ln_delta_cn;
  //SP
  feature_array[3] = get_match_score(match, SP);
  //lnrSP
  feature_array[4] = logf(get_match_rank(match, SP));
  //dM
  feature_array[5] = weight_diff;
  //absdM
  feature_array[6] = fabsf(weight_diff);
  //Mass
  feature_array[7] = get_spectrum_neutral_mass(match->spectrum, match->charge);
  //ionFrac
  feature_array[8] = match->b_y_ion_match;
  //lnSM
  feature_array[9] = match->ln_experiment_size;
  
  //peptide cleavage info.
  if(match->overall_type == TRYPTIC){
    feature_array[10] = TRUE;
    feature_array[11] = TRUE;
  }
  else if(match->overall_type == N_TRYPTIC){
    feature_array[10] = TRUE;
  }
  else if(match->overall_type == C_TRYPTIC){
    feature_array[11] = TRUE;
  }
  
  //get the missed cleave sites
  feature_array[12] = get_peptide_missed_cleavage_sites(match->peptide);
  
  //pepLen
  feature_array[13] = get_peptide_length(match->peptide);
  
  //set charge
  if(match->charge == 1){
    feature_array[14] = TRUE;
  }
  else if(match->charge == 2){
    feature_array[15] = TRUE;
  }
  else if(match->charge == 3){
    feature_array[16] = TRUE;
  }
  
  //run specific features
  feature_array[17] = get_match_collection_hash(match_collection, match->peptide);
  
  src_iterator = new_peptide_src_iterator(match->peptide);
  
  //iterate overall parent proteins
  //find largest numProt and pepSite among the parent proteins
  while(peptide_src_iterator_has_next(src_iterator)){
    peptide_src = peptide_src_iterator_next(src_iterator);
    protein = get_peptide_src_parent_protein(peptide_src);
    protein_idx = get_protein_protein_idx(protein);
    
    //numProt
    if(feature_array[18] < get_match_collection_protein_counter(match_collection,
                                                                protein_idx)){
      feature_array[18] = get_match_collection_protein_counter(match_collection, 
                                                               protein_idx);
    }
    
    //pepSite
    if(feature_array[19] < get_match_collection_protein_peptide_counter(match_collection,
                                                                             protein_idx)){
      feature_array[19] = get_match_collection_protein_peptide_counter(match_collection,
                                                                       protein_idx);      
    }
  }
  
  //now check that no value is with in infinity
  int check_idx = 0;
  for(; check_idx < 20; ++check_idx){
    if(feature_array[check_idx] <= -INFINITY ||
       feature_array[check_idx] >= INFINITY){
      carp(CARP_ERROR, "Percolator feature out of bounds: %d, with value %.2f", check_idx, feature_array[check_idx]);
    }
  }
    
  free_peptide_src_iterator(src_iterator);
  
  return feature_array;
}

/**
 *
 *\returns a match object that is parsed from the serialized result file
 */
MATCH_T* parse_match(
  FILE* result_file,  ///< the result file to parse PSMs -in
  DATABASE_T* database ///< the database to which the peptides are created -in
  //int num_top_match  ///< number of top PSMs serialized per spectrum -in
  )
{
  MATCH_T* match = new_match();
  SPECTRUM_T* spectrum = NULL;
  PEPTIDE_T* peptide = NULL;
  
  //this is a post_process match object
  match->post_process_match = TRUE;
  int score_type_idx = 0;
  
  //parse score, ranks of the match    
  if((peptide = parse_peptide(result_file, database, TRUE))== NULL){
    carp(CARP_ERROR, "failed to parse peptide");
    exit(1);
  }
  
  //parse each score and rank of match
  for(; score_type_idx < _SCORE_TYPE_NUM; ++score_type_idx){
    fread(&(match->match_scores[score_type_idx]), sizeof(float), 1, result_file);
    fread(&(match->match_rank[score_type_idx]), sizeof(int), 1, result_file);
  }
  
  //parse spectrum
  if((spectrum = parse_spectrum_binary(result_file))== NULL){
    carp(CARP_ERROR, "failed to");
  }
  
  //spectrum specific features
  fread(&(match->b_y_ion_match), sizeof(float), 1, result_file);

  //parse match peptide overall trypticity
  fread(&(match->overall_type), sizeof(PEPTIDE_TYPE_T), 1, result_file);
  
  //parse if match is it null_peptide?
  fread(&(match->null_peptide), sizeof(BOOLEAN_T), 1, result_file);

  //assign fields
  match->peptide_sequence = NULL;
  match->spectrum = spectrum;
  match->peptide = peptide;
  
  return match;  
}


/****************************
 * match get, set methods
 ***************************/


/**
 * Returns a heap allocaated peptide sequence of the PSM
 * User must free the sequence.
 *
 * Go to top README for N,C terminus tryptic feature info.
 *
 *\returns the match peptide sequence, returns NULL if no sequence avaliable
 */
char* get_match_sequence(
  MATCH_T* match ///< the match to work -in
  )
{
  //if the match is post_process_match and a null peptide you cannot get sequence
  if(match->post_process_match && match->null_peptide){
    carp(CARP_ERROR, "cannot retrieve null peptide sequence for post_process_match");
    return NULL;
  }
  
  // if peptide sequence is cached
  // return copy of cached peptide sequence
  if(match->peptide_sequence != NULL){
    return my_copy_string(match->peptide_sequence);
  }
  
  //if not cached go generate the sequence

  //First, is this a null peptide?
  //Then must use the shuffled sequence
  if(match->null_peptide){
    //generate the shuffled peptide sequence
    match->peptide_sequence = 
      generate_shuffled_sequence(match->peptide, match->overall_type);
  }
  else{
    // just go parse it out from protein, no need to shuffle
    match->peptide_sequence = get_peptide_sequence(match->peptide);
  }
  
  return my_copy_string(match->peptide_sequence);
}


/**
 * Must ask for score that has been computed
 *\returns the match_mode score in the match object
 */
float get_match_score(
  MATCH_T* match, ///< the match to work -in  
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  )
{
  return match->match_scores[match_mode];
}

/**
 * sets the match score
 */
void set_match_score(
  MATCH_T* match, ///< the match to work -out
  SCORER_TYPE_T match_mode, ///< the working mode (SP, XCORR) -in
  float match_score ///< the score of the match -in
  )
{
  match->match_scores[match_mode] = match_score;
}

/**
 * Must ask for score that has been computed
 *\returns the match_mode rank in the match object
 */
int get_match_rank(
  MATCH_T* match, ///< the match to work -in  
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  )
{
  return match->match_rank[match_mode];
}

/**
 * sets the rank of the match
 */
void set_match_rank(
  MATCH_T* match, ///< the match to work -in  
  SCORER_TYPE_T match_mode, ///< the working mode (SP, XCORR) -in
  int match_rank ///< the rank of the match -in
  )
{
  match->match_rank[match_mode] = match_rank;
}

/**
 *\returns the spectrum in the match object
 */
SPECTRUM_T* get_match_spectrum(
  MATCH_T* match ///< the match to work -in  
  )
{
  return match->spectrum;
}

/**
 * sets the match spectrum
 */
void set_match_spectrum(
  MATCH_T* match, ///< the match to work -out
  SPECTRUM_T* spectrum  ///< the working spectrum -in
  )
{

  match->spectrum = spectrum;
}

/**
 *\returns the peptide in the match object
 */
PEPTIDE_T* get_match_peptide(
  MATCH_T* match ///< the match to work -in  
  )
{
  return match->peptide;
}

/**
 * sets the match peptide, and also determines the overall trypticity of the peptide
 * Go to top README for N,C terminus tryptic feature info.
 */
void set_match_peptide(
  MATCH_T* match, ///< the match to work -out
  PEPTIDE_T* peptide  ///< the working peptide -in
  )
{
  //first set peptide
  match->peptide = peptide;
  
  //now set peptide overall trypticity
  PEPTIDE_SRC_ITERATOR_T* src_iterator = 
    new_peptide_src_iterator(peptide);

  PEPTIDE_SRC_T* peptide_src = NULL;

  //iterate overall parent proteins
  //determine the match overal trypticity
  //for more detail look at README at top
  while(peptide_src_iterator_has_next(src_iterator)){
    peptide_src = peptide_src_iterator_next(src_iterator);

    //now if its tryptic we are done
    if(get_peptide_src_peptide_type(peptide_src) == TRYPTIC){
      match->overall_type = TRYPTIC;
      break;
    }
    else if(get_peptide_src_peptide_type(peptide_src) == N_TRYPTIC){
      //now there're at least tryptic on both side
      // thus now we set as Tryptic and are done
      if(match->overall_type == C_TRYPTIC){
        match->overall_type = TRYPTIC;
        break;
      }
      else{
        match->overall_type = N_TRYPTIC;
      }
    }
    else if(get_peptide_src_peptide_type(peptide_src) == C_TRYPTIC){
      // now there're at least tryptic on both side
      // thus now we set as Tryptic and are done
      if(match->overall_type == N_TRYPTIC){
        match->overall_type = TRYPTIC;
        break;
      }
      else{
        match->overall_type = C_TRYPTIC;
      }
    } 
    else{
      match->overall_type = NOT_TRYPTIC;
    }
  }
  
  free_peptide_src_iterator(src_iterator);
  
}

/**
 * sets the match if it is a null_peptide match
 */
void set_match_null_peptide(
  MATCH_T* match, ///< the match to work -out
  BOOLEAN_T is_null_peptid  ///< is the match a null peptide? -in
  )
{
  match->null_peptide = is_null_peptid;  
}

/**
 * gets the match if it is a null_peptide match
 *\returns TRUE if match is null peptide, else FALSE
 */
BOOLEAN_T get_match_null_peptide(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->null_peptide;
}

/**
 * sets the match charge
 */
void set_match_charge(
  MATCH_T* match, ///< the match to work -out
  int charge  ///< the charge of spectrum -in
  )
{
  match->charge = charge;
}

/**
 * gets the match charge
 */
int get_match_charge(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->charge;
}

/**
 * sets the match delta_cn
 */
void set_match_delta_cn(
  MATCH_T* match, ///< the match to work -out
  float delta_cn  ///< the delta cn value of PSM -in
  )
{
  match->delta_cn = delta_cn;
}

/**
 * gets the match delta_cn
 */
float get_match_delta_cn(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->delta_cn;
}

/**
 * sets the match ln_delta_cn
 */
void set_match_ln_delta_cn(
  MATCH_T* match, ///< the match to work -out
  float ln_delta_cn  ///< the ln delta cn value of PSM -in
  )
{
  match->ln_delta_cn = ln_delta_cn;
}

/**
 * gets the match ln_delta_cn
 */
float get_match_ln_delta_cn(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->ln_delta_cn;
}

/**
 * sets the match ln_experiment_size
 */
void set_match_ln_experiment_size(
  MATCH_T* match, ///< the match to work -out
  float ln_experiment_size ///< the ln_experiment_size value of PSM -in
  )
{
  match->ln_experiment_size = ln_experiment_size;
}

/**
 * gets the match ln_experiment_size
 */
float get_match_ln_experiment_size(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->ln_experiment_size;
}

/**
 * sets the match b_y_ion_match
 */
void set_match_b_y_ion_match(
  MATCH_T* match, ///< the match to work -out
  float b_y_ion_match ///< the b_y_ion_match value of PSM -in
  )
{
  match->b_y_ion_match = b_y_ion_match; 
}

/**
 * gets the match b_y_ion_match
 */
float get_match_b_y_ion_match(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->b_y_ion_match;
}



/**
 *Increments the pointer count to the match object
 */
void increment_match_pointer_count(
  MATCH_T* match ///< the match to work -in  
  )
{
  ++match->pointer_count;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

