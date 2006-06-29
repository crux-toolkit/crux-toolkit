/*****************************************************************************
 * \file peptide.c
 * $Revision: 1.1 $
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
 */
struct peptide {
  PROTEIN_T* protein;       ///< The protein from whence the peptide came
  unsigned short int start; ///< The starting idx of the peptide 
  unsigned char length;     ///< The length of the peptide
};

/**
 * \struct residue_iterator
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

