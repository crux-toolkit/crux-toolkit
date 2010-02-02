#include "objects.h"
#include "carp.h"
//#include "modifications.h"

/**
 * \file mass.h 
 * $Revision: 1.18 $
 * \brief Provides constants and methods for calculating mass
 *****************************************************************************/
#ifndef _MASS_H
#define _MASS_H

/**
 * Mass of ammonia
 */
#define MASS_NH3_MONO  17.02655 ///< Mass of NH3 (monoisotopic)
#define MASS_NH3_AVERAGE 17.03056 ///< Mass of NH3 (average)


// mass of water
#define MASS_H2O_MONO 18.01056 ///< Mass of water (monoisotopic)
#define MASS_H2O_AVERAGE 18.0153 ///< Mass of water (average)

/**
 * Mass of hydrogen
 */
#define MASS_H_MONO 1.0078246 ///< Mass of hydrogen (monoisotopic)
#define MASS_H_AVERAGE 1.00794 ///< Mass of hydrogen (average)

// FIXME, change in spectrum, peptide to be able to pick mono, average
#define MASS_H 1.0078246 ///< mass of hydrogen
#define MASS_PROTON  1.00727646677 ///< mass of proton
/**
 * Mass of oxygen
 */
#define MASS_O 16.0013 ///< mass of oxygen

/**
 * Mass of carbon monoxide
 */
#define MASS_CO_MONO 27.9949 ///< Mass of  (monoisotopic)
#define MASS_CO_AVERAGE 28.0101 ///< Mass of  (average)

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \returns The mass of the given amino acid.
 */
FLOAT_T get_mass_amino_acid(
  char amino_acid, ///< the query amino acid -in
  MASS_TYPE_T mass_type ///< the isotopic mass type (AVERAGE, MONO) -in
  );

/**
 * \returns The mass of the given amino acid.
 */
FLOAT_T get_mass_mod_amino_acid(
  MODIFIED_AA_T amino_acid, ///< the query amino acid -in
  MASS_TYPE_T mass_type ///< the isotopic mass type (AVERAGE, MONO) -in
  );


/**
 * \returns The average mass of the given amino acid.
 */
FLOAT_T get_mass_amino_acid_average(
  char amino_acid ///< the query amino acid -in
  );

/**
 * \returns The average mass of the given amino acid.
 */
FLOAT_T get_mass_mod_amino_acid_average(
  MODIFIED_AA_T amino_acid ///< the query amino acid -in
  );

/**
 * \returns The monoisotopic mass of the given amino acid.
 */
FLOAT_T get_mass_amino_acid_monoisotopic(
  char amino_acid ///< the query amino acid -in
  );

/**
 * \returns The monoisotopic mass of the given amino acid.
 */
FLOAT_T get_mass_mod_amino_acid_monoisotopic(
  MODIFIED_AA_T amino_acid ///< the query amino acid -in
  );

/**
 * increase the amino acid mass for both mono and average
 */
void increase_amino_acid_mass(
  char amino_acid, ///< the query amino acid -in
  FLOAT_T update_mass ///< the mass amount to update for the amino acid -in
  );

#ifdef __cplusplus
}
#endif


#endif
