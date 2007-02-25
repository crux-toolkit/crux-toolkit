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

#define MAX_NUMBER_PEPTIDES 1000000 //What to set?

/**
 *\struct match_collection
 *\brief An object that contains match objects with a given spectrum and peptide database
 */
struct match_collection{
  MATCH_T* match[MAX_NUMBER_PEPTIDES]; ///< array of match object
  BOOLEAN_T scored_type[_SCORE_TYPE_NUM]; ///< has the score type been computed in each match
  int match_total; ///< total_match_count
  SCORER_TYPE_T last_sorted; ///< the last type the match has been sorted(if -1, then unsorted, if ever change the order must change to -1)
  BOOLEAN_T iterator_lock; ///< is there a iterator been curretly created?, if TRUE cannot manipulate match collection
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
  int max_rank,     ///< max number of top rank matches to keep from SP -in
  GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator ///< peptide iteartor to use, must set it first before use
  );

BOOLEAN_T score_match_collection_xcorr(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  SPECTRUM_T* spectrum, ///< the spectrum to match peptides -in
  int charge       ///< the charge of the spectrum -in
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
 * return the top max_rank matches, by score_type(SP, XCORR);
 *\returns a new match_collection object that is scored by score_type and contains the top max_rank matches
 */
MATCH_COLLECTION_T* new_match_collection_spectrum(
 SPECTRUM_T* spectrum, ///< the spectrum to match peptides -in
 int charge,       ///< the charge of the spectrum -in
 int max_rank,     ///< max number of top rank matches to keep from SP -in
 SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
 )
{
  MATCH_COLLECTION_T* match_collection = allocate_match_collection();
  
  //create a generate peptide iterator
  GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator =  //FIXME use neutral_mass, might chage to pick
    new_generate_peptides_iterator_sp(get_spectrum_neutral_mass(spectrum, charge));
  
  //score SP match_collection
  if(!score_match_collection_sp(match_collection, spectrum, charge, max_rank, peptide_iterator)){
    carp(CARP_ERROR, "failed to score match collection for SP");
  }
  
  //should we score for XCORR?
  if(score_type == XCORR){
    /* implement later
    if(!score_match_collection_xcorr(match_collection, spectrum, charge)){
      carp(CARP_ERROR, "failed to score match collection for XCORR");
    }
    */
  }

  //free generate_peptides_iterator
  free_generate_peptides_iterator(peptide_iterator);
  
  return match_collection;
}

/**
 * create a new match collection from spectrum
 * return the top max_rank matches, by score_type(SP, XCORR);
 * uses a provided peptide iterator, MUST be a mutable iterator
 * Sets the iterator before useage.
 *\returns a new match_collection object that is scored by score_type and contains the top max_rank matches
 */
MATCH_COLLECTION_T* new_match_collection_spectrum_with_peptide_iterator(
 SPECTRUM_T* spectrum, ///< the spectrum to match peptides -in
 int charge,       ///< the charge of the spectrum -in
 int max_rank,     ///< max number of top rank matches to keep from SP -in
 SCORER_TYPE_T score_type, ///< the score type (SP, XCORR) -in
 GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator ///< peptide iteartor to use, must set it first before use
 )
{
  MATCH_COLLECTION_T* match_collection = allocate_match_collection();
  
  //get perameters
  float neutral_mass = get_spectrum_neutral_mass(spectrum, charge);
  double mass_window = get_double_parameter("mass-window", 3);
  double min_mass = neutral_mass - mass_window;
  double max_mass = neutral_mass + mass_window;
    
  carp(CARP_DEBUG,"searching peptide in %.2f ~ %.2f", min_mass, max_mass); 
  
  //set the generate_peptides_iterator for the next round of peptides
  set_generate_peptides_mutable(peptide_iterator, max_mass, min_mass);
  
  //score SP match_collection
  if(!score_match_collection_sp(match_collection, spectrum, charge, max_rank, peptide_iterator)){
    carp(CARP_ERROR, "failed to score match collection for SP");
  }
  
  //should we score for XCORR?
  if(score_type == XCORR){
    /* implement later
    if(!score_match_collection_xcorr(match_collection, spectrum, charge)){
      carp(CARP_ERROR, "failed to score match collection for XCORR");
    }
    */
  }

  //free generate_peptides_iterator
  //free_generate_peptides_iterator(peptide_iterator);
  
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
    //qsort_match(match_collection->match, match_collection->match_total, (void *)compare_match_xcorr);
    //match_collection->last_sorted = XCORR;
    return TRUE;
  case DOTP:
    //implement later
    return FALSE;
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
    carp(CARP_ERROR, "match collection has not been set to populate match rank");
    return FALSE;
  }
  
  //set match rank for all match objects
  int match_index = 0;
  for(; match_index < match_collection->match_total; ++match_index){
    set_match_rank(match_collection->match[match_index], score_type, match_index+1);
  }
  
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
  int max_rank,     ///< max number of top rank matches to keep from SP -in
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
    
    //create a new match
    match = new_match();
    
    //set all fields in match
    set_match_score(match, SP, score);
    set_match_peptide(match, peptide);
    set_match_spectrum(match, spectrum);
    
    //check if enough space for peptide match
    if(match_collection->match_total >= MAX_NUMBER_PEPTIDES){
      carp(CARP_ERROR, "peptide count exceed max match limit: %d", MAX_NUMBER_PEPTIDES);
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

  //DEBUG, print total peptided scored so far
  carp(CARP_INFO, "total peptide scored for sp: %d", match_collection->match_total);
  
  //free heap
  free_scorer(scorer);
  free_ion_constraint(ion_constraint);
  //free_generate_peptides_iterator(peptide_iterator);

  //save only the top max_rank matches, sort and free the other matches
  truncate_match_collection(match_collection, max_rank, SP);
  
  //now that the match_collection is sorted, populate the rank of each match object
  if(!populate_match_rank_match_collection(match_collection, SP)){
    carp(CARP_ERROR, "failed to populate match rank in match_collection");
    free_match_collection(match_collection);
    exit(-1);
  }
  
  //yes, we have now scored for the match-mode: SP
  match_collection->scored_type[SP] = TRUE;
    
  return TRUE;
}

/**
 * \returns TRUE, if successfully scores matches for xcorr
 */
/*
BOOLEAN_T score_match_collection_xcorr(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  SPECTRUM_T* spectrum, ///< the spectrum to match peptides -in
  int charge       ///< the charge of the spectrum -in
  )
{
  //implement later
  return TRUE;
}
*/

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
  SCORER_TYPE_T score_type, ///< the score type to iter6ate (SP, XCORR) -in
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
  
  //ok lock up match collection
  match_collection->iterator_lock = TRUE;
  
  //allocate a new match iterator
  MATCH_ITERATOR_T* match_iterator = (MATCH_ITERATOR_T*)mycalloc(1, sizeof(MATCH_ITERATOR_T));
  
  //set items
  match_iterator->match_collection = match_collection;
  match_iterator->match_mode = SP;
  match_iterator->match_idx = 0;
  match_iterator->match_total = match_collection->match_total;

  //only sort if requested and match collection is not already sorted
  if(sort_match && match_collection->last_sorted != score_type){
    //sort match collection by score type
    if(!sort_match_collection(match_collection, score_type)){
      carp(CARP_ERROR, "failed to sort match collection");
      free_match_collection(match_collection);
      free(match_iterator);
      exit(-1);
    }
  }
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

