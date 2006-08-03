/*****************************************************************************
 * \file index.c
 * $Revision: 1.2 $
 * \brief: Object for representing an index of a database
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"
#include "crux-utils.h"
#include "peptide.h"
#include "protein.h"
#include "index.h"
#include "carp.h"
#include "objects.h"

#define INDEX_MAX_FILES 10000


/* 
 * How does the create_index (the routine) work?
 *
 * - Creates a database object from a the fasta file member variable.
 *
 * - Depending on the size of the database, determines how many passes it
 *   will need to make over the database peptides (first implementation
 *   will only make one pass, later implementation can make multiple passes
 *   with multiple iterators)
 *
 * - From the directory and mass resolution member variables, 
 *   creates the list of filenames necessary for storing the peptides 
 *   
 * - Creates a database_peptide_iterator object from the database 
 *   and peptide constraint objects
 *
 * - Then starts iterating through peptides 
 *
 *    - From the peptide mass determines which file to place it in, and
 *      then calls the serialize_peptide method
 *
 *      - The serialize peptide method is an object method specifically
 *        designed for this purpose, and writes a peptide in a format from
 *        which it can be reconstructed
 *
 * LATER At some point we will index each of the peptide files too (in a TOC), 
 * so that we can retrieve them rapidly later. 
 */

/* 
 * How does generate_peptides (the executable) with --from-index work?
 *
 * - Given a fasta file, looks for the corresponding index on disk. Fail if it
 *   can't find the corresponding index files.
 *
 * - Instantiates an index object from the fasta file with new_index.
 *
 * - Attempts to parse_index the index.
 *
 * - Instantiates an index_peptide_iterator from the index, according to
 *   constraints passed on the command line to generate_peptides.
 *
 * - Then starts iterating through peptides, which are being loaded from
 *   disk, and outputs them as before
 * 
 */

/**
 * \struct index
 * \brief A index of a database
 */
struct index{
  DATABASE_T* database; ///< The database that has been indexed.
  char* directory; ///< The directory containing the indexed files
  char* filenames[INDEX_MAX_FILES]; ///< The files that contain the peptides
  PEPTIDE_CONSTRAINT_T* constraint; ///< Constraint which these peptides satisfy
  BOOLEAN_T is_on_disk; ///< Does this index exist on disk yet?
};    

/**
 * The main index method. Does all the heavy lifting, creating files
 * serializing peptides, etc. The index directory itself should have 
 * a standard suffix (e.g. cruxidx), so that a given fasta file will have
 * an obvious index location.
 *
 * Note: create an index .info file as to map masses to files and to store 
 * all information that was used to create this index.
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T create_index(
  INDEX_T* index ///< An allocated index
  );

/*
 * Private methods
 */

/*
 * Returns the index filename appropriate for this peptide
 */
char* get_peptide_file_name(
    INDEX_T* index,
    PEPTIDE_T* peptide
    );
/*
 * Iterators
 */

/**
 * \struct index_peptide_iterator
 * \brief An iterator to iterate over the peptides in a database
 */
struct index_peptide_iterator{
  INDEX_T* index; ///< The index object which we are iterating over
  PEPTIDE_CONSTRAINT_T* constraint; ///< The constraints which peptides should satisfy
  unsigned int peptide_idx; ///< The (non-object) index of the current peptide.
  char* current_index_filename; ///< The current file that we are reading from
  BOOLEAN_T has_next; ///< Is there another peptide?
};    

/// TODO new_index_peptide_iterator should to see if the index exists on disk

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
