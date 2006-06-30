/**
 * \file protein.h 
 * $Revision: 1.3 $
 * \brief Object for representing one protein sequence.
 *****************************************************************************/
#ifndef PROTEIN_H 
#define PROTEIN_H

#include <stdio.h>
#include "utils.h"
#include "objects.h"
#include "peptide.h"

/* CHRIS This is probably an object for which you can crib code for from an outside source. Even from in-house (like Charles).*/

/**
 * \returns An (empty) protein object.
 */
PROTEIN_T* allocate_protein(void);

/**
 * \returns An (empty) protein object.
 */
PROTEIN_T* allocate_protein(void);

/**
 * \returns A new protein object.
 */
PROTEIN_T* new_protein(
  char*         id, ///< The protein sequence id.
  char*   sequence, ///< The protein sequence.
  int       length, ///< The length of the protein sequence.
  char* annotation  ///< Optional protein annotation. 
  );         

/**
 * Frees an allocated protein object.
 */
void free_protein(PROTEIN_T* protein);

/**
 * Prints a protein object to file.
 */
void print_protein(PROTEIN_T* protein, FILE* file);

/**
 * Copies protein object src to dest.
 */
void copy_protein(
  PROTEIN_T* src,
  PROTEIN_T* dest);

/**
 * Parses a protein from an open (FASTA) file.
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T parse_protein_fasta_file(
  PROTEIN_T* protein,
  FILE* file);

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 * FIXME Chris, could you create the get and set methods for the object fields?
 */

/**
 * Additional get and set methods
 */

/*int get_protein_peptides(PROTEIN_T* protein);*/
/*int get_protein_peptides(PROTEIN_T* protein);*/


/**
 * Instantiates a new peptide_iterator from a peptide.
 * \returns a PEPTIDE_ITERATOR_T object.
 */
PEPTIDE_ITERATOR_T* new_peptide_iterator(PROTEIN_T* protein);        

/**
 * Frees an allocated peptide_iterator object.
 */
void free_peptide_iterator(PEPTIDE_ITERATOR_T* peptide_iterator);

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T peptide_iterator_has_next(PEPTIDE_ITERATOR_T* peptide_iterator);

/**
 * \returns The next peptide (a character) in the peptide.
 */
PEPTIDE_T* peptide_iterator_next(PEPTIDE_ITERATOR_T* peptide_iterator);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
