/*****************************************************************************
 * \file ion.c
 * $Revision: 1.1 $
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
  int cleavage_idx; ///< index of peptide amide that fragments to form this ion
  int charge; ///< the ion charge
  PEPTIDE_T* peptide; ///< the peptide that fragments to form this ion
  ION_MODIFICATION_T modifications[MAX_MODIFICATIONS]; ///< an array of ion modifications
  float ion_mass;   ///< The (neutral) mass of the ion. 
};


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

