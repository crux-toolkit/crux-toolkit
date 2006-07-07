/**
 * \file mass.h 
 * $Revision: 1.4 $
 * \brief Provides constants and methods for calculating mass
 *****************************************************************************/
#ifndef _MASS_H
#define _MASS_H

#include "peptide.h"

/**
 * Mass of ammonia
 */
#define MASS_NH3 17.0306

/**
 * Mass of water
 */
#define MASS_H2O 18.0153

/**
 * Mass of hydrogen
 */
#define MASS_H 1.007

/**
 * Mass of oxygen
 */
#define MASS_O 16.0013

/**
 * Mass of carbon monoxide
 */
#define MASS_CO 28.0101

/**
 * \returns The mass of the given peptide.
 */
float calc_peptide_mass(PEPTIDE_T* peptide);

#endif
