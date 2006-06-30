/**
 * \file peptide.h 
 * $Revision: 1.2 $
 * \brief Object for representing one peptide.
 *****************************************************************************/
#ifndef PEPTIDE_H 
#define PEPTIDE_H

#include "utils.h"
#include <stdio.h>
#include "objects.h"
#include "protein.h"

/**
 * \typedef PEPTIDE_TYPE_T The peptide type, with regard to trypticity.
 */
typedef enum _peptide_type { TRYPTIC, PARTIALLY_TRYPTIC, NON_TRYPTIC} 
  PEPTIDE_TYPE_T;

/**
 * \returns An (empty) peptide object.
 */
PEPTIDE_T* allocate_peptide(void);

/* FIXME what about peptides in multiple proteins? */
/**
 * \returns A new peptide object, populated with the user specified parameters.
 */
PEPTIDE_T* new_peptide(
  PROTEIN_T* my_protein,       ///< The protein from whence the peptide came
  unsigned short int start, ///< The starting idx of the peptide 
  unsigned char length,     ///< The length of the peptide
  double peptide_mass
);
  

/* possible additional fields and methods
 * - charge (possibly "unknown"),
 * - mass-to-charge. 
 *   The function get_proteins 
 *   returns a list of proteins in which the peptide occurs 
 *   fxn: add_protein
 *   */


/**
 * Frees an allocated peptide object.
 */
void free_peptide (PEPTIDE_T* peptide);

/**
 * Prints a peptide object to file.
 */
void print_peptide(PEPTIDE_T* peptide, FILE* file);

/**
 * Copies peptide object src to dest.
 */
void copy_peptide(
  PEPTIDE_T* src,
  PEPTIDE_T* dest);

/**
 * Parses a peptide from file.
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_peptide_file(
  PEPTIDE_T* peptide,
  FILE* file);

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 * FIXME Chris, could you create the get and set methods for the object fields?
 */

/**
 * Additional get and set methods
 */

/**
 * Iterator
 */

/**
 * Instantiates a new residue_iterator from a peptide.
 * \returns a RESIDUE_ITERATOR_T object.
 */
RESIDUE_ITERATOR_T* new_residue_iterator(PEPTIDE_T* peptide);        

/**
 * Frees an allocated residue_iterator object.
 */
void free_residue_iterator(RESIDUE_ITERATOR_T* residue_iterator);

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional residues to iterate over, FALSE if not.
 */
BOOLEAN_T residue_iterator_has_next(RESIDUE_ITERATOR_T* residue_iterator);

/**
 * \returns The next residue (a character) in the peptide.
 */
char residue_iterator_next(RESIDUE_ITERATOR_T* residue_iterator);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
