/**
 * \file protein.h 
 * $Revision: 1.25 $
 * \brief Object for representing one protein sequence.
 *****************************************************************************/
#ifndef PROTEIN_H 
#define PROTEIN_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "utils.h"
#include "objects.h"
#include "peptide.h"
#include "peptide_src.h"
#include "carp.h"
#include "peptide_constraint.h"

/* CHRIS This is probably an object for which you can crib code for from an outside source. Even from in-house (like Charles).*/

/**
 * \returns An (empty) protein object.
 */
PROTEIN_T* allocate_protein(void);

/**
 * \returns A new protein object(heavy).
 */
PROTEIN_T* new_protein(
  char*         id, ///< The protein sequence id.
  char*   sequence, ///< The protein sequence.
  unsigned int length, ///< The length of the protein sequence.
  char* annotation,  ///< Optional protein annotation.  -in
  unsigned long int offset, ///< The file location in the source file in the database -in
  unsigned int protein_idx, ///< The index of the protein in it's database. -in
  DATABASE_T* database ///< the database of its origin
  );         

/**
 * \returns A new light protein object.
 */
PROTEIN_T* new_light_protein(
  unsigned long int offset, ///< The file location in the source file in the database -in
  unsigned int protein_idx ///< The index of the protein in it's database. -in
  );

/**
 * convert light protein to heavy, by parsing all the sequence from fasta file
 * \returns TRUE if successfully converts the protein to heavy 
 */
BOOLEAN_T protein_to_heavy(
  PROTEIN_T* protein ///< protein to convert to heavy -in 
  );
                         
/**
 * covert heavy protein back to light
 * \returns TRUE if successfully converts the protein to light
 */
BOOLEAN_T protein_to_light(
  PROTEIN_T* protein ///< protein to convert back to light -in 
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
 * Parses a protein from an memory mapped binary fasta file
 * the protein_idx field of the protein must be added before or after you parse the protein
 * \returns TRUE if success. FALSE is failure.
 * protein must be a heap allocated
 * 
 * Assume memmap pointer is set at beginning of protein
 * Assume protein binary format
 * <int: id length><char: id><int: annotation length><char: annotation><int: sequence length><char: sequence>
 *
 * modifies the *memmap pointer!
 */
BOOLEAN_T parse_protein_binary_memmap(
  PROTEIN_T* protein, ///< protein object to fill in -out
  void** memmap ///< a pointer to a pointer to the memory mapped binary fasta file -in
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
unsigned int get_protein_length(
  PROTEIN_T* protein ///< the query protein -in 
);

/**
 * sets the id of the protein
 */
void set_protein_length(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  unsigned int length ///< the length to add -in
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
 * sets the offset of the protein in the fasta file
 */
void set_protein_offset(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  unsigned long int offset ///< The file location in the source file in the database -in
  );

/**
 *\returns the offset the protein
 */
unsigned long int get_protein_offset(
  PROTEIN_T* protein ///< the query protein -in 
  );

/**
 * sets the protein_idx (if, idx=n, nth protein in the fasta file)
 */
void set_protein_protein_idx(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  unsigned int protein_idx ///< The index of the protein in it's database. -in
  );

/**
 *\returns the protein_idx field
 */
unsigned int get_protein_protein_idx(
  PROTEIN_T* protein ///< the query protein -in 
  );


/**
 * sets the is_light field (is the protein a light protein?)
 */
void set_protein_is_light(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  BOOLEAN_T is_light ///< is the protein a light protein? -in
  );

/**
 *\returns TRUE if the protein is light protein
 */
BOOLEAN_T get_protein_is_light(
  PROTEIN_T* protein ///< the query protein -in 
  );

/**
 * sets the database for protein
 */
void set_protein_database(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  DATABASE_T*  database ///< Which database is this protein part of -in
  );

/**
 *\returns Which database is this protein part of
 */
DATABASE_T* get_protein_database(
  PROTEIN_T* protein ///< the query protein -in 
  );

/**
 * prints a binary representation of the protein
 * 
 * FORMAT
 * <int: id length><char: id><int: annotation length><char: annotation><int: sequence length><char: sequence>
 *
 * make sure when rading the binary data, add one to the length so that it will read in the terminating char as well
 */
void serialize_protein(
  PROTEIN_T* protein, ///< protein to print as binary -in
  FILE* file ///< output stream -out
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

/**
 *\returns the protein that the iterator was created on
 */
PROTEIN_T* get_protein_peptide_iterator_portein(
  PROTEIN_PEPTIDE_ITERATOR_T* protein_peptide_iterator ///< working protein_peptide_iterator -in
  );
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#ifdef __cplusplus
}
#endif
#endif
