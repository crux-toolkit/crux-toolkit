/**
 * \file protein_index.h
 * $Revision: 1.3 $
 * \brief Object for creating a protein index
 *****************************************************************************/
#ifndef PROTEIN_INDEX_H
#define PROTEIN_INDEX_H

#include <stdio.h>
#include "utils.h"
#include "objects.h"
#include "peptide.h"
#include "Protein.h"
#include "carp.h"
#include "peptide_constraint.h"
#include "sorter.h"

/**
 * creates a protein index on to the output_file
 * \returns TRUE if successfully creates a protein index, else false
 */
BOOLEAN_T create_protein_index(
  char* fasta_file ///< input fasta file stream -in
  );

/**
 * \returns An (empty) protein_index object.
 */
PROTEIN_INDEX_T* allocate_protein_index(void);

/**
 * creates a protein_index that contains the offset and protein index of the protein
 * in the fasta file.
 *\returns a new protein_index object
 */
PROTEIN_INDEX_T* new_protein_index(
  unsigned long int offset, ///< The file location in the source file in the database
  unsigned int protein_idx ///< The index of the protein in it's database.
  );

/**
 *
 * free a protein index object
 */
void free_protein_index(
  PROTEIN_INDEX_T* protein_index  ///< the protein index to free
  );                        

/**
 * input is the fasta file name which the protein index
 * should have been created.
 * or if creating binary fasta file, is that already on disk?
 *
 *\returns TRUE if protein index or binary fasta file is on disk, else FALSE
 */
BOOLEAN_T protein_index_on_disk(
  char* fasta_file, ///< input fasta file -in
  BOOLEAN_T is_binary ///< are we looking for the binary fasta file? or preotein index
  );

/**
 * protein index iterator
 * the protein index iterator parses the protein index file
 * for each protein index one by one and returns a new protein_index object
 */

/**
 * 
 *\returns a new heap allocated protein index iterator
 */
PROTEIN_INDEX_ITERATOR_T* new_protein_index_iterator(
  char* fasta_file ///< input fasta file -in
);

/**
 * Frees the allocated protein index iterator
 */
void free_protein_index_iterator(
  PROTEIN_INDEX_ITERATOR_T* protein_index_iterator ///< the iterator to free -in
  );

/**
 *
 *\returns TRUE if there is another protein index to return, else FALSE
 */
BOOLEAN_T protein_index_iterator_has_next(
  PROTEIN_INDEX_ITERATOR_T* protein_index_iterator ///< the iterator of interest -in
);

/**
 *
 *\return the next protein index in the protein index file
 */
Protein* protein_index_iterator_next(
  PROTEIN_INDEX_ITERATOR_T* protein_index_iterator ///< the iterator of interest -in
);


/**
 * creates a binary fasta file on to the output_file
 * \returns TRUE if successfully creates a binary fasta file, else false
 */
BOOLEAN_T create_binary_fasta(
  char* fasta_file  ///< input fasta file -in
  );

/**
 * creates a binary fasta file on to the output_file in currenty directory
 * sets the output file name to the pointer passed in as argument
 * \returns TRUE if successfully creates a binary fasta file, else false
 */
BOOLEAN_T create_binary_fasta_in_cur(
  char* fasta_file_w_path, ///< input fasta file with full path -in
  char* fasta_filename, ///< input fasta a file, only filename -in
  char** output_file_name ///< get output filename -out
  );

/**
 * wrapper for create_binary_fasta_file so that two filenames are
 * passed instead of a filename and a filestream.  Eventually should
 * merge to one method
 */
BOOLEAN_T create_binary_fasta_here(
  char* fasta_filename,
  char* binary_filename
);

/**
 * Heap allocated char*, user must free
 *\returns the binary fasta name which was created from the given fasta file
 */
char* get_binary_fasta_name(
  char* fasta_file  ///< input fasta file -in                            
  );


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
