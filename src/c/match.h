/**
 * \file match.h
 * $Revision: 1.4 $ 
 * \brief Object for given a peptide and a spectrum, generate a perliminary score(ex, Sp)
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

/**
 * \returns a new memory allocated match
 */
MATCH_T* new_match(void);

/**
 * free the memory allocated match
 */
void free_match(
  MATCH_T* match, ///< the match to free -in
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
  MATCH_T* match, ///< the match to print -in  
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  );

/**
 * sets the match score
 */
void set_match_score(
  MATCH_T* match, ///< the match to print -out
  SCORER_TYPE_T match_mode, ///< the working mode (SP, XCORR) -in
  float match_score ///< the score of the match -in
  );

/**
 * Must ask for score that has been computed
 *\returns the match_mode rank in the match object
 */
float get_match_rank(
  MATCH_T* match, ///< the match to print -in  
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  );

/**
 * sets the rank of the match
 */
void set_match_rank(
  MATCH_T* match, ///< the match to print -in  
  SCORER_TYPE_T match_mode, ///< the working mode (SP, XCORR) -in
  int match_rank ///< the rank of the match -in
  );

/**
 *\returns the spectrum in the match object
 */
SPECTRUM_T* get_match_spectrum(
  MATCH_T* match ///< the match to print -in  
  );

/**
 * sets the match spectrum
 */
void set_match_spectrum(
  MATCH_T* match, ///< the match to print -out
  SPECTRUM_T* spectrum  ///< the working spectrum -in
  );

/**
 *\returns the peptide in the match object
 */
PEPTIDE_T* get_match_peptide(
  MATCH_T* match ///< the match to print -in  
  );

/**
 * sets the match peptide
 */
void set_match_peptide(
  MATCH_T* match, ///< the match to print -out
  PEPTIDE_T* peptide  ///< the working peptide -in
  );
