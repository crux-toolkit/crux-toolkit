/*****************************************************************************
 * \file match.c
 * AUTHOR: Chris Park
 * CREATE DATE: 11/27 2006
 * DESCRIPTION: Object for matching a peptide and a spectrum, generate a perliminary score(ex, Sp)
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
  SPECTRUM_T* spectrum; ///< the spectrum we are scoring with
  PEPTIDE_T* peptide;  ///< the peptide we are scoring
  float match_scores[_SCORE_TYPE_NUM]; ///< the scoring result array (use enum_type SCORER_TYPE_T to index)
  int match_rank[_SCORE_TYPE_NUM];  ///< the rank of scoring result (use enum_type SCORER_TYPE_T to index)
  int pointer_count; ///< the number of pointers to this match object (when reach 0, free memory)
};

/**
 * \returns a new memory allocated match
 */
MATCH_T* new_match(void){
  MATCH_T* match = (MATCH_T*)mycalloc(1, sizeof(MATCH_T));
  ++match->pointer_count;
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
  //only free match when pointer count reaches
  if(--match->pointer_count == 0){
    //free peptide
    free_peptide(match->peptide);
    free(match);  
  }
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


/****************************
 * match get, set methods
 ***************************/

/**
 * Must ask for score that has been computed
 *\returns the match_mode score in the match object
 */
float get_match_score(
  MATCH_T* match, ///< the match to print -in  
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  )
{
  return match->match_scores[match_mode];
}

/**
 * sets the match score
 */
void set_match_score(
  MATCH_T* match, ///< the match to print -out
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
float get_match_rank(
  MATCH_T* match, ///< the match to print -in  
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  )
{
  return match->match_rank[match_mode];
}

/**
 * sets the rank of the match
 */
void set_match_rank(
  MATCH_T* match, ///< the match to print -in  
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
  MATCH_T* match ///< the match to print -in  
  )
{
  return match->spectrum;
}

/**
 * sets the match spectrum
 */
void set_match_spectrum(
  MATCH_T* match, ///< the match to print -out
  SPECTRUM_T* spectrum  ///< the working spectrum -in
  )
{
  match->spectrum = spectrum;
}

/**
 *\returns the peptide in the match object
 */
PEPTIDE_T* get_match_peptide(
  MATCH_T* match ///< the match to print -in  
  )
{
  return match->peptide;
}

/**
 * sets the match peptide
 */
void set_match_peptide(
  MATCH_T* match, ///< the match to print -out
  PEPTIDE_T* peptide  ///< the working peptide -in
  )
{
  match->peptide = peptide;
}
