/**
 * \file match.h
 * $Revision: 1.7 $ 
 * \brief Object for given a peptide and a spectrum, generate a preliminary score(ex, Sp)
 ****************************************************************************/
#ifndef MATCH_H
#define MATCH_H

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

/**
 * \returns a new memory allocated match
 */
MATCH_T* new_match(void);

/**
 * free the memory allocated match
 */
void free_match(
  MATCH_T* match ///< the match to free -in
  );

/**
 * sort the match array with the corresponding compare method
 */
void qsort_match(
  MATCH_T** match_array, ///< the match array to sort -in  
  int match_total,  ///< the total number of match objects -in
  void* compare_method ///< the compare method to use -in
  );

/**
 * compare two matches, used for qsort
 * \returns the difference between sp score in match_a and match_b
 */
int compare_match_sp(
  MATCH_T* match_a, ///< the first match -in  
  MATCH_T* match_b  ///< the scond match -in
  );

/**
 * compare two matches, used for qsort
 * \returns the difference between xcorr score in match_a and match_b
 */
int compare_match_xcorr(
  MATCH_T* match_a, ///< the first match -in  
  MATCH_T* match_b  ///< the scond match -in
  );

/**
 * print the information of the match
 */
void print_match(
  MATCH_T* match, ///< the match to print -in  
  FILE* file, ///< output stream -out
  BOOLEAN_T output_sequence, ///< should I output peptide sequence -in
  SCORER_TYPE_T output_mode  ///< the output mode -in
);

/****************************
 * match get, set methods
 ***************************/

/**
 * Must ask for score that has been computed
 *\returns the match_mode score in the match object
 */
float get_match_score(
  MATCH_T* match, ///< the match to work -in  
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  );

/**
 * sets the match score
 */
void set_match_score(
  MATCH_T* match, ///< the match to work -out
  SCORER_TYPE_T match_mode, ///< the working mode (SP, XCORR) -in
  float match_score ///< the score of the match -in
  );

/**
 * Must ask for score that has been computed
 *\returns the match_mode rank in the match object
 */
float get_match_rank(
  MATCH_T* match, ///< the match to work -in  
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  );

/**
 * sets the rank of the match
 */
void set_match_rank(
  MATCH_T* match, ///< the match to work -in  
  SCORER_TYPE_T match_mode, ///< the working mode (SP, XCORR) -in
  int match_rank ///< the rank of the match -in
  );

/**
 *\returns the spectrum in the match object
 */
SPECTRUM_T* get_match_spectrum(
  MATCH_T* match ///< the match to work -in  
  );

/**
 * sets the match spectrum
 */
void set_match_spectrum(
  MATCH_T* match, ///< the match to work -out
  SPECTRUM_T* spectrum  ///< the working spectrum -in
  );

/**
 *\returns the peptide in the match object
 */
PEPTIDE_T* get_match_peptide(
  MATCH_T* match ///< the match to work -in  
  );

/**
 * sets the match peptide
 */
void set_match_peptide(
  MATCH_T* match, ///< the match to work -out
  PEPTIDE_T* peptide  ///< the working peptide -in
  );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
