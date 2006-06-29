/**
 * \file protein.h 
 * $Revision: 1.1 $
 * \brief Object for representing one protein sequence.
 *****************************************************************************/
#ifndef PROTEIN_H 
#define PROTEIN_H

#include <stdio.h>
#include "utils.h"

/* CHRIS This is probably an object for which you can crib code for from an outside source. Even from in-house (like Charles).*/

/**
 * \typedef PROTEIN_T
 */
typedef struct protein PROTEIN_T;

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

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
