/*****************************************************************************
 * \file peptide.c
 * $Revision: 1.9 $
 * \brief: Object for representing a single peptide.
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "mass.h"
#include "protein.h"
#include "peptide.h"


/**
 * \struct peptide
 * \brief A subsequence of a protein.
 */
struct peptide {
  char* sequence;       ///< A pointer to the peptide sequence.
  unsigned char length; ///< The length of the peptide
};

/**
 * \struct peptide_constraint
 * \brief Object to represent constraints which a peptide may or may not
 * satisfy.
 */
struct peptide_constraint {
  PEPTIDE_TYPE_T*  peptide_type; ///< The type of peptides.
  float min_mass; ///< The minimum mass of the peptide
  float max_mass; ///< The maximum mass of the peptide
  int min_length; ///< The minimum length of the peptide
  int max_length; ///< The maximum length of the peptide
};
 
/**
 * \struct residue_iterator
 * \brief Object to iterate over the residues in a peptide, starting at the
 * first residue of the peptide, and proceeding in order.
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

