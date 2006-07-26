#include "objects.h"
#include "carp.h"

/**
 * \file mass.h 
 * $Revision: 1.11 $
 * \brief Provides constants and methods for calculating mass
 *****************************************************************************/
#ifndef _MASS_H
#define _MASS_H

/**
 * Mass of ammonia
 */
#define MASS_NH3 17.0306

/**
 * Mass of water
 */
#define MASS_H2O_MONO 18.01056
#define MASS_H2O_AVERAGE 18.0153

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
 * \returns The mass of the given amino acid.
 */
float get_mass_amino_acid(
  char amino_acid, ///< the query amino acid -in
  MASS_TYPE_T mass_type ///< the isotopic mass type (AVERAGE, MONO) -in
  );

/**
 * \returns The average mass of the given amino acid.
 */
float get_mass_amino_acid_average(
  char amino_acid ///< the query amino acid -in
  );

/**
 * \returns The monoisotopic mass of the given amino acid.
 */
float get_mass_amino_acid_monoisotopic(
  char amino_acid ///< the query amino acid -in
  );

#endif
