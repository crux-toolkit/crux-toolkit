/**
 * \file index.h 
 * $Revision: 1.25 $
 * \brief Object for representing an index of a index
 *****************************************************************************/
#ifndef INDEX_H 
#define INDEX_H

#include <stdio.h>

#include "utils.h"
#include "objects.h"
#include "peptide.h"
#include "Protein.h"
#include "sorter.h"
#include "carp.h"
#include "peptide_constraint.h"
#include "database.h"

static const int MAX_INDEX_FILES = 1024;


/**
 * clean_up
 *
 */
void clean_up( int dummy );

/**
 * \returns An (empty) index object.
 */
INDEX_T* allocate_index(void);

/**
 * \returns A new index object.
 */
INDEX_T* new_index(
  const char* fasta_filename,  ///< The fasta file
  const char* output_dir,      ///< The name of the new index
  PEPTIDE_CONSTRAINT_T* constraint,///< Constraint which these peptides satisfy
  FLOAT_T mass_range  ///< the range of masses contained in each index file
);         

/**
 * wrapper function, create index object for search purpose
 * from disk. If no crux_index files have been created on disk, returns null
 * \returns A new index object ready for search.
 */
INDEX_T* new_index_from_disk(
  const char* fasta_filename  ///< The fasta file
  //  BOOLEAN_T is_unique ///< only unique peptides? -in
  );

/**
 * Merely increments the index ptr count
 */
INDEX_T* copy_index_ptr(
    INDEX_T* index
    );

/**
 * Get index ptr count
 */
int get_index_pointer_count(
    INDEX_T* index
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
  INDEX_T* index, ///< An allocated index -in/out
  BOOLEAN_T create_text_file ///< Should an ASCII text file be create? -in
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

/**
 * \brief Looks in given directory for a file ending in
 * "-binary-fasta" and returns a heap-allocated string of the full
 * name including the index directory.
 *
 * Exits with error if index_name does not exist, no file
 * *-binary-fasta exists, or more than one *binary-fasta file exists.
 * \returns A string with the name of the existing binary fasta file
 * for this index.
 */

char* get_index_binary_fasta_name(const char* index_name);
/**
 * foo.fasta --> foo_crux_index/foo_binary_fasta
 * \returns the binary fasta file name with crux directory name
 */
char* get_binary_fasta_name_in_crux_dir(
  char* fasta_filename  ///< fasta file name -in
  );

/*********************************************
 * set and get methods for the object fields
 *********************************************/

/**
 *\returns the directory of the index
 * returns a heap allocated new copy of the directory
 * user must free the return directory name
 */
char* get_index_directory(
  INDEX_T* index ///< The index -in
  );

/**
 * sets the directory of the index
 * index->directory must been initiailized
 */
void set_index_directory(
  INDEX_T* index, ///< The index -in
  const char* directory ///< the directory to add -in
  );

/**
 *\returns a pointer to the database
 */
DATABASE_T* get_index_database(
  INDEX_T* index ///< The index -in
  );

/**
 * sets the database of the index
 */
void set_index_database(
  INDEX_T* index, ///< The index -in
  DATABASE_T* database ///< The database that has been indexed. -in
  );

/**
 *\returns a pointer to the peptides constraint
 */
PEPTIDE_CONSTRAINT_T* get_index_constraint(
  INDEX_T* index ///< The index -in
  );

/**
 * sets the peptides constraint
 */
/*
void set_index_constraint(
  INDEX_T* index, ///< The index -in
  PEPTIDE_CONSTRAINT_T* constraint ///< Constraint which these peptides satisfy -in
  );
*/

/**
 * \brief Sets the peptide search constraint to be used by the
 * generate_peptides_iterator.  Makes a copy of the constraint pointer.
 * Deletes any existing search constraint. 
 */
void set_index_search_constraint(
  INDEX_T* index, ///< The index -in
  PEPTIDE_CONSTRAINT_T* constraint ///< Constraint for the next iterator
  );

/**
 *\returns TRUE if index files are on disk else FALSE
 */
BOOLEAN_T get_index_on_disk(
  INDEX_T* index ///< The index -in
  );

/**
 * sets the on disk field of index
 */
void set_index_on_disk(
  INDEX_T* index, ///< The index -in
  BOOLEAN_T on_disk ///< Does this index exist on disk yet? -in
  );

/**
 *\returns the range of mass that each index file should be partitioned into
 */
FLOAT_T get_index_mass_range(
  INDEX_T* index ///< The index -in
  );

/**
 * sets the mass_range field of index
 */
void set_index_mass_range(
  INDEX_T* index, ///< The index -in
  FLOAT_T mass_range  ///< the range of mass that each index file should be partitioned into -in
  );


/**
 *\returns TRUE if only allow unique peptides else FALSE
 */
BOOLEAN_T get_index_is_unique(
  INDEX_T* index ///< The index -in
  );

/**
 * sets the is_unique field
 */
void set_index_is_unique(
  INDEX_T* index, ///< The index -in
  BOOLEAN_T is_unique ///< do you allow duplicate peptides? -in
  );

int get_index_num_proteins( INDEX_T* index );

/***********************************************
 * Iterators index
 ***********************************************/

/**
 * Instantiates a new peptide_iterator from an index, which returns peptides
 * that obey peptide constraint. At first will only accept constraints
 * that will require reading in one file (e.g a 1m/z range). 
 */
INDEX_PEPTIDE_ITERATOR_T* new_index_peptide_iterator(
  INDEX_T* index ///< The index -in
  // BOOLEAN_T seq ///< output sequence -in
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

/**********************************************************************
 * wrapper, for generate_peptides_iterator, cast back to original type
 ***********************************************************************/
/**
 * Frees an allocated index_peptide_iterator object.
 */
void void_free_index_peptide_iterator(
    void* index_peptide_iterator  ///< the iterator to free -in
    );

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T void_index_peptide_iterator_has_next(
    void* index_peptide_iterator ///< the iterator of interest -in
    );

/**
 * \returns The next peptide in the index.
 */
PEPTIDE_T* void_index_peptide_iterator_next(
    void* index_peptide_iterator ///< the iterator of interest -in
    );


/***********************************************
 * index_filtered_peptide_iterator
 ***********************************************/

/**
 * Instantiates a new index_filtered_peptide_iterator from a index.
 * \returns a new heap allocated index_filtered_peptide_iterator object
 */
INDEX_FILTERED_PEPTIDE_ITERATOR_T* new_index_filtered_peptide_iterator(
  INDEX_T* index ///< The index object which we are iterating over -in
  );

/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
PEPTIDE_T* index_filtered_peptide_iterator_next(
  INDEX_FILTERED_PEPTIDE_ITERATOR_T* index_filtered_peptide_iterator ///< the index_filtered_peptide_iterator to initialize -in
  );

/**
 * The basic iterator functions.
 * check to see if the index_filtered_peptide_iterator has more peptides to return
 *\returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T index_filtered_peptide_iterator_has_next(
  INDEX_FILTERED_PEPTIDE_ITERATOR_T* index_filtered_peptide_iterator ///< the index_filtered_peptide_iterator to initialize -in
  );

/**
 * Frees an allocated index_filtered_peptide_iterator object.
 */
void free_index_filtered_peptide_iterator(
  INDEX_FILTERED_PEPTIDE_ITERATOR_T* index_filtered_peptide_iterator ///< the index_filtered_peptide_iterator to initialize -in
  );

/**********************************************************************
 * wrapper, for generate_peptides_iterator, cast back to original type
 ***********************************************************************/

/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
PEPTIDE_T* void_index_filtered_peptide_iterator_next(
  void* index_filtered_peptide_iterator ///< the index_filtered_peptide_iterator to initialize -in
  );

/**
 * The basic iterator functions.
 * check to see if the index_filtered_peptide_iterator has more peptides to return
 *\returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T void_index_filtered_peptide_iterator_has_next(
  void* index_filtered_peptide_iterator ///< the index_filtered_peptide_iterator to initialize -in
  );

/**
 * Frees an allocated index_filtered_peptide_iterator object.
 */
void void_free_index_filtered_peptide_iterator(
  void* index_filtered_peptide_iterator ///< the index_filtered_peptide_iterator to initialize -in
  );



/***********************************************
 * Iterators BIN
 ***********************************************/


/**
 * Instantiates a new bin_peptide_iterator from a gvien bin file handler.
 * \returns a new heap allocated bin_peptide_iterator object
 */
BIN_PEPTIDE_ITERATOR_T* new_bin_peptide_iterator(
  INDEX_T* index, ///< The index object which we are iterating over -in
  FILE* file, ///< the bin to parse peptides
  BOOLEAN_T use_array  ///< should I use array peptide_src or link list when parsing peptides -in
  );

/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
PEPTIDE_T* bin_peptide_iterator_next(
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator ///< the bin_peptide_iterator to get peptide -in
  );

/**
 * The basic iterator functions.
 * check to see if the bin_peptide_iterator has more peptides to return
 *\returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T bin_peptide_iterator_has_next(
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator ///< the bin_peptide_iterator to initialize -in
  );

/**
 * Frees an allocated bin_peptide_iterator object.
 */
void free_bin_peptide_iterator(
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator ///< the iterator to free -in
  );


/***********************************************
 * Iterators sorted BIN
 ***********************************************/

/**
 * Instantiates a new sorted_bin_peptide_iterator from a gvien bin file handler.
 * \returns a new heap allocated sorted_bin_peptide_iterator object
 */
BIN_SORTED_PEPTIDE_ITERATOR_T* new_bin_sorted_peptide_iterator(
  INDEX_T* index, ///< The index object which we are iterating over -in
  FILE* file,///< the working file handler to the bin -in
  unsigned int peptide_count ///< the total peptide count in the bin -in
  );

/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
PEPTIDE_T* bin_sorted_peptide_iterator_next(
  BIN_SORTED_PEPTIDE_ITERATOR_T* bin_sorted_peptide_iterator ///< the bin_peptide_iterator to get peptide -in
  );

/**
 * The basic iterator functions.
 * check to see if the bin_sorted_peptide_iterator has more peptides to return
 *\returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T bin_sorted_peptide_iterator_has_next(
  BIN_SORTED_PEPTIDE_ITERATOR_T* bin_sorted_peptide_iterator ///< the bin_peptide_iterator to initialize -in
  );

/**
 * Frees an allocated bin_peptide_iterator object.
 */
void free_bin_sorted_peptide_iterator(
  BIN_SORTED_PEPTIDE_ITERATOR_T* bin_sorted_peptide_iterator ///< the iterator to free -in
  );


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
