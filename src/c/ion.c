/*****************************************************************************
 * \file ion.c
 * $Revision: 1.2 $
 * \brief: Object for representing a single ion.
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "objects.h"
#include "ion.h"

#define MAX_MODIFICATIONS 5

/**
 * \struct ion
 * \brief An object for representing a (fragment) ion of a peptide.
 */
struct ion {
  ION_TYPE_T type;  ///< type of the ion 
  int cleavage_idx; ///< index of peptide amide that fragments to form this ion, starting from the N-term end 
  // N.b. this is different than the b1,y1 index, in that it always starts
  // from the N-term
  int charge; ///< the ion charge
  char* peptide; ///< the peptide sequence that fragments to form this ion
  ION_MODIFICATION_T modification_counts[MAX_MODIFICATIONS]; ///< an array of the number of different ion modifications
  float ion_mass;   ///< The (neutral) mass of the ion. 
};


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

