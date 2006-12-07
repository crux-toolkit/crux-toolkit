/*****************************************************************************
 * \file score_peptide_iterator
 * AUTHOR: Chris Park
 * CREATE DATE: 11/27 2006
 * DESCRIPTION: Object for given a peptide and a spectrum, generate a perliminary score(ex, Sp)
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
#include "generate_peptides_iterator.h" 

#define MAX_NUMBER_PEPTIDES 500000 //What to set?

/**
 *\struct match
 *\brief An object that stores the score & rank for each match of spectrum & peptide
 */
struct match{
  PEPTIDE_T* peptide;  ///< the peptide we are scoring
  float match_scores[_SCORE_TYPE_NUM]; ///< the scoring result array (use enum_type SCORER_TYPE_T to index)
  int match_rank[_SCORE_TYPE_NUM];  ///< the rank of scoring result (use enum_type SCORER_TYPE_T to index)
};

/**
 *\struct match iterator
 *\brief An object that generates score match objects with a given spectrum and peptide database
 */
struct match_iterator{
  MATCH_T* match[MAX_NUMBER_PEPTIDES]; ///< array of match object
  SCORER_TYPE_T match_mode; ///< the current score working mode (SP, XCORR)
  BOOLEAN_T scored_type[_SCORE_TYPE_NUM]; ///< has the score type been computed in each match
  SPECTRUM_T* spectrum; ///< the spectrum to score
  GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator; ///< the peptide database to score against spectrum
  int match_idx; ///< current match to return
  int match_total; ///< total_match_count
};

/**
 * \returns a new memory allocated match
 */
MATCH_T* new_match(void){
  MATCH_T* match = (MATCH_T*)mycalloc(1, sizeof(MATCH_T));
  return match;
}

/**
 * free the memory allocated match
 * must provide the match iterator to use the correct method to free peptide contained in match
 */
void free_match(
  MATCH_T* match, ///< the match to free -in
  MATCH_ITERATOR_T* match_iterator ///< the working match iterator -in
  )
{
  //free peptide, must free with the through the free method provided by iterator
  free_peptide_produced_by_iterator(match_iterator->peptide_iterator, match->peptide);
  free(match);  
}

/**
 * create a new memory allocated match iterator
 * creates a new the generate_peptides_iterator inside the match_iterator
 *\returns a new memory allocated match iterator
 */
MATCH_ITERATOR_T* new_match_iterator(
  SPECTRUM_T* spectrum ///< the spectrum to match peptides -in
  )
{
  //allocate a new match iterator
  MATCH_ITERATOR_T* match_iterator = (MATCH_ITERATOR_T*)mycalloc(1, sizeof(MATCH_ITERATOR_T));
  
  //loop over to set all score type to FALSE
  int score_type_idx = 0;
  for(; score_type_idx < _SCORE_TYPE_NUM; ++score_type_idx){
    match_iterator->scored_type[score_type_idx] = FALSE;
  }

  //set working spectrum
  match_iterator->spectrum = spectrum;

  //set to nothing
  match_iterator->match_mode = -1;

  //create peptide iterator
  match_iterator->peptide_iterator = new_generate_peptides_iterator();
  match_iterator->match_idx = 0;
  match_iterator->match_total = 0;

  return match_iterator;
}


/**
 * compare two matches, used for qsort
 * \returns the difference between sp score in match_a and match_b
 */
int compare_match_sp(
  MATCH_T* match_a, ///< the first match -in  
  MATCH_T* match_b  ///< the scond match -in
)
{
  //might have to worry about cases below 1 and -1
  return (int)(match_b->match_scores[SP] - match_a->match_scores[SP]);
}

/**
 * given a match_iterator, generates sp scores in the match objects
 * Also, since this is the first socring method required,
 * it creates new match object for each scoring pair
 * \returns TRUE if all scores been calculated in match iterator, else FALSE
 */
BOOLEAN_T score_match_iterator_sp(
  MATCH_ITERATOR_T* match_iterator ///< the match iterator to score -out                               
  )
{
  int total_peptides = 0;
  GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator = match_iterator->peptide_iterator;
  int peptide_charge = get_int_parameter("charge", 1);
  char* peptide_sequence = NULL;
  MATCH_T* match = NULL;
  float score = 0;
  SCORER_T* scorer = NULL;
  PEPTIDE_T* peptide = NULL;
  ION_SERIES_T* ion_series = NULL;

  //set ion constraint to sequest settings
  ION_CONSTRAINT_T* ion_constraint = 
    new_ion_constraint_sequest_sp(peptide_charge); 
  
  //create new scorer
  scorer = new_scorer(SP);  
    
  //iterate over all peptides
  while(generate_peptides_iterator_has_next(peptide_iterator)){
    ++total_peptides;
    peptide = generate_peptides_iterator_next(peptide_iterator);
    peptide_sequence = get_peptide_sequence_pointer(peptide);
    
    //create new ion series
    ion_series = new_ion_series(peptide_sequence, peptide_charge, ion_constraint);
    
    //now predict ions
    predict_ions(ion_series);
    
    //calculates the Sp score
    score = score_spectrum_v_ion_series(scorer, match_iterator->spectrum, ion_series);
    
    //create a new match
    match = new_match();
    
    //set all fields in match
    match->match_scores[SP] = score;
    match->peptide = peptide;
    
    //check if enough space for peptide match
    if(match_iterator->match_total >= MAX_NUMBER_PEPTIDES){
      carp(CARP_ERROR, "peptide count exceed max match limit: %d", MAX_NUMBER_PEPTIDES);
      //free heap
      free_ion_series(ion_series);
      free_scorer(scorer);
      free_ion_constraint(ion_constraint);
      return FALSE;
    }
    
    //add a new match to array
    match_iterator->match[match_iterator->match_total] = match;
    
    //increment total match count
    ++match_iterator->match_total;

    free_ion_series(ion_series);
  }
  
  //free heap
  free_scorer(scorer);
  free_ion_constraint(ion_constraint);

  //yes, we have now scored for the match-mode: SP
  match_iterator->scored_type[SP] = TRUE;
  
  return TRUE;
}

/**
 * given a match_iterator, generates scores in the match objects for the given score type
 * scores all the matches with correct match_mode(score type)
 * \returns TRUE if all scores been calculated in match iterator, else FALSE
 */
BOOLEAN_T score_match_iterator(
  MATCH_ITERATOR_T* match_iterator, ///< the match iterator to score -out                               
  SCORER_TYPE_T match_mode ///< the mode to score (MATCH_SP, MATCH_XCORR) -in
  )
{
  if(match_mode == SP){
    return score_match_iterator_sp(match_iterator);
  }
  else if(match_mode == XCORR){
    /* add this later
      return score_match_iterator_xcorr(match_iterator);
    */
  }
  
  //can't find correct score type
  carp(CARP_ERROR, "cannot find correct scoring type");
  return FALSE;
}

/**
 * sets the match_iterator mode.
 * If not already scored in iterator, creates the scorer for the correct mode and claculates the scores for,
 * each match object, then sort the match structs so that return in decreasing score order.
 * MUST, use this method to set match_iterator ready for a given score type before iterating through scores
 * \returns TRUE if match iterator is successfully set to the correct mode
 */
BOOLEAN_T set_mode_match_iterator(
  MATCH_ITERATOR_T* match_iterator, ///< the match iterator to set -out
  SCORER_TYPE_T match_mode ///< the mode to set (MATCH_SP, MATCH_XCORR) -in
  )
{
  //the match mode has already been set
  if(match_iterator->match_mode == match_mode){
    carp(CARP_INFO, "the match iterator is already set to correct match mode");
    carp(CARP_INFO, "starting from begining of match list");
    //set the match index to start
    match_iterator->match_idx = 0;
    
    return TRUE;
  }
  
  //set the match mode
  match_iterator->match_mode = match_mode;
  
  //have we already scored this match mode?
  if(!match_iterator->scored_type[match_mode]){
    //score the match iterator, only socres are filled, rankings are still balnk
    if(!score_match_iterator(match_iterator, match_mode)){
      carp(CARP_ERROR, "failed to score peptide, spectrum for match iterator");
      return FALSE;
    }
  }
  
  if(match_mode == SP){
    int sp_match_idx = 0;

    //sort the match to the correct order for the return
    qsort(match_iterator->match, match_iterator->match_total, sizeof(MATCH_T), (void *)compare_match_sp);

    //set SP ranking in match
    for(; sp_match_idx < match_iterator->match_total; ++sp_match_idx){
      match_iterator->match[sp_match_idx]->match_rank[SP] = sp_match_idx + 1;
    }
  }
  else{
    //add other sorting methods...
  }
  
  //set match index to start
  match_iterator->match_idx = 0;

  //good ready to go!
  return TRUE;
}

/**
 * Does the match_iterator have another match struct to return?
 * MUST set the iterator to correct mode before calling this method
 *\returns TRUE, if match iter has a next match, else False
 */
BOOLEAN_T mode_match_iterator_has_next(
  MATCH_ITERATOR_T* match_iterator, ///< the working  match iterator -in
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  )
{
  //Is the match mode correctly set?
  if(match_iterator->match_mode != match_mode){
    carp(CARP_ERROR, "The match iterator is not set correctly");
    exit(1);
  }
  
  return (match_iterator->match_idx < match_iterator->match_total);
}

/**
 * return the next match struct!
 * MUST set the iterator to correct mode before initialially calling this method
 *\returns the match in decreasing score order for the match_mode(SCORER_TYPE_T)
 */
MATCH_T* mode_match_iterator_next(
  MATCH_ITERATOR_T* match_iterator, ///< the working match iterator -in
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  )
{
  //Is the match mode correctly set?
  if(match_iterator->match_mode != match_mode){
    carp(CARP_ERROR, "The match iterator is not set correctly");
    exit(1);
  }
  
  return match_iterator->match[match_iterator->match_idx++];
}

/**
 * free the memory allocated iterator
 */
void free_match_iterator(
  MATCH_ITERATOR_T* match_iterator ///< the match iterator to free
  )
{
  int free_idx = 0; 
  //free all matches
  for(; free_idx < match_iterator->match_total; ++free_idx){
    free_match(match_iterator->match[free_idx], match_iterator);
  }
  
  //free peptide iterator, free matches first, before freeing up iterator!
  free_generate_peptides_iterator(match_iterator->peptide_iterator);
  
  //free iterator
  free(match_iterator);
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
  //print according to the output mode
  switch (output_mode) {
  case SP:
    fprintf(file, "%d\t%.2f\t", match->match_rank[SP], match->match_scores[SP]);
    
    //should I print sequence
    if(output_sequence){
      fprintf(file, "%s\n", get_peptide_sequence_pointer(match->peptide));
    }
    
    break;
  case XCORR:
    //fill in
    break;
  case DOTP:
    //fill in
    break;
  }
}

