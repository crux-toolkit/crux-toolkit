/*****************************************************************************
 * \file peptide.c
 * $Revision: 1.3 $
 * \brief: Object for representing a single peptide.
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "protein.h"
#include "peptide.h"


/**
 * \struct peptide
 * A subsequence of a protein.
 */
struct peptide {
  char* sequence;       ///< The protein sequence that contains the peptide.
  unsigned short int start; ///< The starting idx of the peptide 
  unsigned char length;     ///< The length of the peptide
};

/**
 * \struct residue_iterator
 * Object to iterate over the residues in a peptide.
 */
struct residue_iterator {
  PEPTIDE_T*  peptide; ///< The peptide whose residues to iterate over.
  int     residue_idx; ///< The index of the current peak
};

  
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

