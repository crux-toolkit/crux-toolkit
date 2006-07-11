/*****************************************************************************
 * \file protein_peptide_association.c
 * $Revision: 1.1 $
 * \brief: Object for mapping a peptide to it's parent protein.
 ****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "mass.h"
#include "objects.h"

/**
 * \struct protein_peptide_association
 * \brief object for mapping a peptide to it's parent protein.
 */
struct protein_peptide_association{
  PEPTIDE_TYPE_T peptide_type; ///< the peptide type for the corresponding protein
  PROTEIN_T* parent_protein; ///< the parent of this preptide
  int start; ///< start index of the peptide in the protein sequence
  PROTEIN_PEPTIDE_ASSOCIATION_T* next_association; ///< the next protein_peptide_association 
};



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
