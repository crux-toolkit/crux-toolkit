/*****************************************************************************
 * \file protein.c
 * $Revision: 1.6 $
 * \brief: Object for representing a single protein.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"
#include "protein.h"

/*
 * Constants
 */
#define PROTEIN_ID_LENGTH 100
#define PROTEIN_SEQUENCE_LENGTH 10000
#define PROTEIN_ANNOTATION_LENGTH 100

/**
 * \struct protein 
 * \brief A protein sequence.
 */
struct protein {
  char*         id; ///< The protein sequence id.
  char*   sequence; ///< The protein sequence.
  int       length; ///< The length of the protein sequence.
  char* annotation; ///< Optional protein annotation. 
};    

/**
 * \struct peptide_iterator
 * \brief Object to iterate over the peptides within a protein.
 */
struct peptide_iterator {
  PROTEIN_T* protein; ///< The spectrum whose peaks to iterate over. 
  unsigned short int cur_start; ///< Start in protein of the current peptide.
  unsigned char cur_length; ///< The length of the current peptide.
};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

