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
#define MASS_H2O_MONO 18.010564684 ///< Mass of water (monoisotopic)
#define MASS_H2O_AVERAGE 18.0153 ///< Mass of water (average)

/**
 * Mass of hydrogen
 */
#define MASS_H_MONO 1.0078246 ///< Mass of hydrogen (monoisotopic)
#define MASS_H_AVERAGE 1.00794 ///< Mass of hydrogen (average)

/* As for the constants, these are for supporting mono-isotopic and
 * average masses.  When I was collaborating with Pragya with
 * cross-linking code, she wanted more accurate calculations on the
 * mass of the peptides and the mass of the spectrum neutral mass.
 * MASS_PROTON was added in order to improve the accuracy of the
 * calculated spectrum mass and the precision of the mono-isotopic
 * masses of the amino acids were extended.  When we are going to
 * support accurate masses for the precursor ions, I think that
 * MASS_PROTON will be the correct one to use for the spectrum
 * mass. --Sean McIlwain, 8 November 2010 */

// FIXME, change in spectrum, peptide to be able to pick mono, average
#define MASS_H       1.00782503207 ///< mass of hydrogen
#define MASS_PROTON  1.00727646677 ///< mass of proton
#define MASS_NEUTRON 1.00866491600 ///< mass of neutron

/**
 * Mass of oxygen
 */
#define MASS_O 16.0013 ///< mass of oxygen

#define MASS_OH 17.00274

/**
 * Mass of carbon monoxide
 */
#define MASS_CO_MONO 27.9949 ///< Mass of  (monoisotopic)
#define MASS_CO_AVERAGE 28.0101 ///< Mass of  (average)


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
 * Finds the modification identifier associated with the given mass
 * shift.  Can be the identifier from a single modfification or from
 * multiple modficiations to the same residue.  The returned
 * identifier can be used to modify a MODIFIED_AA_T so that it has the
 * given mass shift. 
 */
MODIFIED_AA_T get_mod_identifier(FLOAT_T mass_shift);

/**
 * increase the amino acid mass for both mono and average
 */
void increase_amino_acid_mass(
  char amino_acid, ///< the query amino acid -in
  FLOAT_T update_mass ///< the mass amount to update for the amino acid -in
  );

/**
 * \brief Populates the array aa_mod_masses with the mass change of
 * all possible combinations of aa_mods.  Gets the list of aa_mods
 * from parameter.c.  For example, if mod1 and a mass change of 50 and
 * mod5 has a mass change of 10, then the entry at index 
 * (binary 00010001 =) 17 is (50 + 10=) 60.
 */
void initialize_aa_mod_combinations_array();


#endif
