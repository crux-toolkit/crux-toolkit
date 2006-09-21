/**
 * \file protein_index.h
 * $Revision: 1.1 $
 * \brief Object for creating a protein index
 *****************************************************************************/
#ifndef PROTEIN_INDEX_H
#define PROTEIN_INDEX_H

#include <stdio.h>
#include "utils.h"
#include "objects.h"
#include "peptide.h"
#include "protein.h"
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
 *\returns TRUE if protein index is on disk, else FALSE
 */
BOOLEAN_T protein_index_on_disk(
  char* fasta
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
PROTEIN_T* protein_index_iterator_next(
  PROTEIN_INDEX_ITERATOR_T* protein_index_iterator ///< the iterator of interest -in
);



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
