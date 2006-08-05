/**
 * \file index.h 
 * $Revision: 1.2 $
 * \brief Object for representing an index of a index
 *****************************************************************************/
#ifndef INDEX_H 
#define INDEX_H

#include <stdio.h>
#include "utils.h"
#include "objects.h"
#include "peptide.h"
#include "protein.h"
#include "carp.h"
#include "peptide_constraint.h"

/**
 * \returns An (empty) index object.
 */
INDEX_T* allocate_index(void);

/**
 * \returns A new index object.
 */
INDEX_T* new_index(
  char* fasta_filename,  ///< The fasta file
  PEPTIDE_CONSTRAINT_T* constraint,  ///< Constraint which these peptides satisfy
  float mass_range,  ///< the range of mass that each index file should be partitioned into
  unsigned int max_size  ///< maximum limit of each index file
);         

/**
 * Frees an allocated index object.
 */
void free_index(INDEX_T* index);

/**
 * Scans the index on disk to 
 * populate fields in index.
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_index(
  INDEX_T* index ///< An allocated index
  );

/**
 * The main index method. Does all the heavy lifting, creating files
 * serializing peptides, etc.
 *
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T create_index(
  INDEX_T* index ///< An allocated index
  );

/**
 * Does this index exist on disk?
 *
 * \returns TRUE if it does. FALSE if it does not.
 */
BOOLEAN_T index_exists(
  INDEX_T* index ///< An allocated index
  );

/**
 * prints all the index files to a file
 * this file is used later by the index object to index the files with a given interval
 * \returns TRUE if it creates a list of infex files. FALSE if it fails.
 */
BOOLEAN_T create_index_files(
  INDEX_T* index, ///< An allocated index
  FILE* file ///< output stream to print
  );

/***********************************************
 * Iterators
 ***********************************************/

/**
 * Instantiates a new peptide_iterator from an index, which returns peptides
 * that obey peptide constraint. At first will only accept constraints
 * that will require reading in one file (e.g a 1m/z range). 
 */
INDEX_PEPTIDE_ITERATOR_T* new_index_peptide_iterator(
  INDEX_T* index, ///< The index -in
  PEPTIDE_CONSTRAINT_T* constraint ///< The constraint to satisfy -in
  );

/**
 * Frees an allocated index_peptide_iterator object.
 */
void free_index_peptide_iterator(
    INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator  ///< the iterator to free -in
    );

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T index_peptide_iterator_has_next(
    INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator ///< the iterator of interest -in
    );

/**
 * \returns The next peptide in the index.
 */
PEPTIDE_T* index_peptide_iterator_next(
    INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator ///< the iterator of interest -in
    );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
