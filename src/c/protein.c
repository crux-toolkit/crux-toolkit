/*****************************************************************************
 * \file protein.c
 * $Revision: 1.10 $
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
 * \struct protein_peptide_iterator
 * \brief Object to iterate over the peptides within a protein in an
 * unspecified order. The peptides should satisfy the constraints specified
 * in the peptide_constraint object.
 * 
 */
struct protein_peptide_iterator {
  PROTEIN_T* protein; ///< The protein whose peptides to iterate over. 
  unsigned short int cur_start; ///< Start in protein of the current peptide.
  unsigned char cur_length; ///< The length of the current peptide.
  unsigned int peptide_idx; ///< The index of the current peptide.
  PEPTIDE_CONSTRAINT_T* peptide_constraint; ///< The type of peptide to iterate over.
};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

