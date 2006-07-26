/**
 * \file protein.h 
 * $Revision: 1.18 $
 * \brief Object for representing one protein sequence.
 *****************************************************************************/
#ifndef PROTEIN_H 
#define PROTEIN_H

#include <stdio.h>
#include "utils.h"
#include "objects.h"
#include "peptide.h"
#include "peptide_src.h"
#include "carp.h"

/* CHRIS This is probably an object for which you can crib code for from an outside source. Even from in-house (like Charles).*/

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
void free_protein(
  PROTEIN_T* protein ///< object to free -in
  );

/**
 * Prints a protein object to file.
 */
void print_protein(
  PROTEIN_T* protein, ///< protein to print -in
  FILE* file ///< output stream -out
  );

/**
 * Copies protein object src to dest.
 * dest must be a heap allocated object 
 */
void copy_protein(
  PROTEIN_T* src,///< protein to copy -in
  PROTEIN_T* dest ///< protein to copy to -out
  );

/**
 * Parses a protein from an open (FASTA) file.
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T parse_protein_fasta_file(
  PROTEIN_T* protein, ///< protein object to fill in -out
  FILE* file ///< fasta file -in
  );

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 * FIXME Chris, could you create the get and set methods for the object fields?
 */

/**
 * Additional get and set methods
 */

/*PEPTIDE_T** get_protein_peptides(PROTEIN_T* protein, PEPTIDE_CONSTRAINT*
 * peptide_constraint);*/

/**
 *\returns the id of the protein
 * returns a heap allocated new copy of the id
 * user must free the return id
 */
char* get_protein_id(
  PROTEIN_T* protein ///< the query protein -in 
);

/**
 *\returns a pointer to the id of the protein
 */
char* get_protein_id_pointer(
  PROTEIN_T* protein ///< the query protein -in 
  );

/**
 * sets the id of the protein
 */
void set_protein_id(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  char* id ///< the sequence to add -in
);

/**
 *\returns the sequence of the protein
 * returns a heap allocated new copy of the sequence
 * user must free the return sequence 
 */
char* get_protein_sequence(
  PROTEIN_T* protein ///< the query protein -in 
);

/**
 *\returns a pointer to the sequence of the protein
 */
char* get_protein_sequence_pointer(
  PROTEIN_T* protein ///< the query protein -in 
  );

/**
 * sets the sequence of the protein
 */
void set_protein_sequence(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  char* sequence ///< the sequence to add -in
);

/**
 *\returns the length of the protein
 */
int get_protein_length(
  PROTEIN_T* protein ///< the query protein -in 
);

/**
 * sets the id of the protein
 */
void set_protein_length(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  int length ///< the length to add -in
);

/**
 *\returns the annotation of the protein
 * returns a heap allocated new copy of the annotation
 * user must free the return annotation
 */
char* get_protein_annotation(
  PROTEIN_T* protein ///< the query protein -in 
);

/**
 * sets the annotation of the protein
 */
void set_protein_annotation(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  char* annotation ///< the sequence to add -in
);

/**
 * Iterator
 * iterates over the peptides given a partent protein and constraints
 */

/**
 * Instantiates a new peptide_iterator from a peptide.
 * \returns a PROTEIN_PEPTIDE_ITERATOR_T object.
 */
PROTEIN_PEPTIDE_ITERATOR_T* new_protein_peptide_iterator(PROTEIN_T* protein, 
    PEPTIDE_CONSTRAINT_T* peptide_constraint);

/**
 * Frees an allocated peptide_iterator object.
 */
void free_protein_peptide_iterator(
    PROTEIN_PEPTIDE_ITERATOR_T* protein_peptide_iterator);

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T protein_peptide_iterator_has_next(
    PROTEIN_PEPTIDE_ITERATOR_T* protein_peptide_iterator);

/**
 * \returns The next peptide in the protein, in an unspecified order
 */
PEPTIDE_T* protein_peptide_iterator_next(
    PROTEIN_PEPTIDE_ITERATOR_T* protein_peptide_iterator);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
