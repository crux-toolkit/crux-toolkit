/*****************************************************************************
 * \file ion_series.c
 * AUTHOR: Chris Park
 * CREATE DATE: 21 Sep 2006
 * DESCRIPTION: code to support working with a series of ions
 * REVISION: $Revision: 1.3 $
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "objects.h"
#include "ion.h"
#include "utils.h"

#define MAX_IONS 10000

/**
 * \struct ion_series
 * \brief An object to represent a series of ions, and organize them!
 * For which additional data structures will be created as needed 
 */
struct ion_series {
  char* peptide; ///< The peptide for this ion series
  int charge; ///< The charge state of the peptide for this ion series
  ION_CONSTRAINT_T* constraint; ///< The constraints which the ions in this series obey
  ION_T* ions[MAX_IONS]; ///< The ions in this series
  int num_ions; ///< the number of ions in this series

};

/**
 * \struct ion_constraint
 * \brief An object to represent the contraints which the ions in this
 * series obey.
 * CHRIS you can add BOOLEANS and other things as needed for implementation
 * of the predict-peptide-ions executable
 */
struct ion_constraint {
  BOOLEAN_T* use_neutral_losses; ///< A boolean to determine if the ions series should include neutral losses
};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
