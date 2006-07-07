/**
 * \file mass.h 
 * $Revision: 1.2 $
 * \brief Provides constants and methods for calculating mass
 *****************************************************************************/
#ifndef _MASS_H
#define _MASS_H

#include "peptide.h"
#define MASS_NH3 17.0306
#define MASS_H2O 18.0153
#define MASS_H 1.007
#define MASS_O 16.0013
#define MASS_CO 28.0101

/**
 * \returns The mass of the given peptide.
 */
float calc_peptide_mass(
    PEPTIDE_T* peptide);

#endif
