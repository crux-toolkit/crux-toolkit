/**
 * \file database.h 
 * $Revision: 1.27 $
 * \brief Object for representing a database of protein sequences.
 *****************************************************************************/
#ifndef DATABASE_H
#define DATABASE_H

#include <stdio.h>
#include "utils.h"
#include "objects.h"
#include "peptide.h"
#include "protein.h"
#include "carp.h"
#include "peptide_constraint.h"
#include "sorter.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \returns An (empty) database object.
 */
DATABASE_T* allocate_database(void);

/**
 * \returns A new database object.
 */
DATABASE_T* new_database(
  char*         filename, ///< The file from which to parse the database. either text fasta file or binary fasta file -in
  BOOLEAN_T is_memmap ///< are we using a memory mapped binary fasta file, thus proteins are all memory mapped -in
  );         

/**
 * Frees an allocated protein object.
 */
void free_database(
  DATABASE_T* database ///< An allocated database -in
  );

/**
 * Prints a database object to file.
 */
void print_database(
  DATABASE_T* database,  ///< database to print -in
  FILE* file    ///< output file stream -out             
  );

/**
 * Parses a database from the file in the filename member variable
 * reads in all proteins in the fasta file and creates a protein object
 * and adds them to the database protein array
 * total proteins in fasta file must not exceed MAX_PROTEIN constant
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_database(
  DATABASE_T* database ///< An allocated database to parse -in
  );

/**
 * \brief Changes a database from one that reads from a fasta file to
 * one that reads from a binary/memmory mapped protein file.
 *
 * If database already has binary source (i.e. is_memmap == TRUE), 
 * returns TRUE.  
 * Opens the fasta file pointed to by filename for reading.  Creates a
 * file with the name given.  Reads in each protein from the text file
 * and serializes it to the output file.  Closes both files.  Changes
 * filename to point to new output file and sets is_memmap to true.
 * Parses the database.
 * \returns TRUE if all processes succeed, else FALSE.
 */
BOOLEAN_T transform_database_text_to_memmap(
  DATABASE_T* database,
  char* binary_protein_filename
  );

/**
 * \returns FALSE if nth protein cannot be parsed or does not exist 
 */
/**
BOOLEAN_T get_database_protein_at_idx(
    DATABASE_T* database, ///< A parsed database object -in
    unsigned int protein_idx, ///< The index of the protein to retrieve -in
    PROTEIN_T** protein   ///< A pointer to a pointer to a PROTEIN object -out
    );
**/

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 * FIXME Chris, could you create the get and set methods for the object fields?
 */

/**
 * Additional get and set methods
 */

/**
 *\returns the filename of the database
 * returns a heap allocated new copy of the filename
 * user must free the return filename
 */
char* get_database_filename(
  DATABASE_T* database ///< the query database -in 
  );

/**
 *\returns the pointer to the filename of the database
 * user must not free or change the filename
 */
char* get_database_filename_pointer(
  DATABASE_T* database ///< the query database -in 
  );

/**
 * sets the filename of the database
 * protein->sequence must been initiailized
 */
void set_database_filename(
  DATABASE_T* database, ///< the database to set it's fields -out
  char* filename ///< the filename to add -in
  );

/**
 *\returns TRUE|FALSE whether the database has been parsed?
 */
BOOLEAN_T get_database_is_parsed(
  DATABASE_T* database ///< the query database -in 
  );

/**
 *\returns the total number of proteins of the database
 */
unsigned int get_database_num_proteins(
  DATABASE_T* database ///< the query database -in 
  );

/**
 *\returns the src FILE* of the database
 */
FILE* get_database_file(
  DATABASE_T* database ///< the query database -in 
  );

/**
 * sets the src FILE* of the database
 */
void set_database_file(
  DATABASE_T* database, ///< the database to set it's fields -out
  FILE* file ///< the src file to add -in
  );

/**
 *\returns the nth protein of the database
 */
PROTEIN_T* get_database_protein_at_idx(
  DATABASE_T* database, ///< the query database -in 
  unsigned int protein_idx ///< The index of the protein to retrieve -in
  );

/**
 * sets the use_light_protein of the database
 */
void set_database_use_light_protein(
  DATABASE_T* database, ///< the database to set it's fields -out
  BOOLEAN_T use ///< should I use the light/heavy functionality?
  );

/**
 *\returns TRUE|FALSE whether the database uses light/heavy
 */
BOOLEAN_T get_database_use_light_protein(
  DATABASE_T* database ///< the query database -in 
  );

/**
 *sets TRUE,FALSE whether the database uses memory mapped
 */
void set_database_memmap(
  DATABASE_T* database, ///< the query database -in 
  BOOLEAN_T is_memmap  ///< is the database memory mapped?
  );


/**
 * increase the pointer_count produced by this database.
 * \returns database pointer
 */
DATABASE_T* copy_database_ptr(
  DATABASE_T* database ///< the query database -in/out
  );

/***********************************************
 * Iterators
 ***********************************************/

/**
 * Instantiates a new database_protein_iterator from a database.
 * \returns a DATABASE_PROTEIN_ITERATOR_T object.
 */
DATABASE_PROTEIN_ITERATOR_T* new_database_protein_iterator(
    DATABASE_T* database ///< the database to create a protein iterator -in
    );        

/**
 * Frees an allocated database_protein_iterator object.
 */
void free_database_protein_iterator(
    DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator  ///< the iterator to free -in
    );

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional proteins to iterate over, FALSE if not.
 */
BOOLEAN_T database_protein_iterator_has_next(
    DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator ///< the iterator of interest -in
    );

/**
 * \returns The next protein in the database.
 */
PROTEIN_T* database_protein_iterator_next(
    DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator ///< the iterator of interest -in
    );

/**
 * \returns the protein to the corresponding protein_idx in the database.
 */
PROTEIN_T* database_protein_iterator_protein_idx(
    DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator, ///< the iterator of interest -in
    unsigned int protein_idx ///< protein_idx to which protein to return -in
    );


/***********************************************
 * database_peptide_Iterators
 ***********************************************/

/**
 * Instantiates a new database_peptide_iterator from a database.
 * \returns a DATABASE_PEPTIDE_ITERATOR_T object.
 */
DATABASE_PEPTIDE_ITERATOR_T* new_database_peptide_iterator(
  DATABASE_T* database, ///< the database of interest -in
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide_constraint to filter peptides -in
  );

/**
 * Frees an allocated database_peptide_iterator object.
 */
void free_database_peptide_iterator(
  DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator ///< the iterator to free -in
  );

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T database_peptide_iterator_has_next(
  DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator ///< the iterator of interest -in
  );

/**
 * \returns The next peptide in the database.
 */
PEPTIDE_T* database_peptide_iterator_next(
  DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator ///< the iterator of interest -in
  );

/**********************************************************************
 * wrapper, for generate_peptides_iterator, cast back to original type
 ***********************************************************************/

/**
 * Frees an allocated database_peptide_iterator object.
 */
void void_free_database_peptide_iterator(
  void* database_peptide_iterator ///< the iterator to free -in
  );

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T void_database_peptide_iterator_has_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  );

/**
 * \returns The next peptide in the database.
 */
PEPTIDE_T* void_database_peptide_iterator_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  );

/***********************************
 * database sorted peptide iterator
 ***********************************/

/**
 * Instantiates a new database_sorted_peptide_iterator from a database.
 * uses a sorted_peptide_iterator as it's engine
 * \returns a DATABASE_SORTED_PEPTIDE_ITERATOR_T object.
 */
DATABASE_SORTED_PEPTIDE_ITERATOR_T* new_database_sorted_peptide_iterator(
  DATABASE_T* database, ///< the database of interest -in
  PEPTIDE_CONSTRAINT_T* peptide_constraint, ///< the peptide_constraint to filter peptides -in
  SORT_TYPE_T sort_type, ///< the sort type for this iterator
  BOOLEAN_T unique ///< only return unique peptides? -in
  );

/**
 * Frees an allocated database_sorted_peptide_iterator object.
 */
void free_database_sorted_peptide_iterator(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_peptide_iterator ///< the iterator to free -in
  );

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T database_sorted_peptide_iterator_has_next(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_peptide_iterator ///< the iterator of interest -in
  );

/**
 * returns each peptide in sorted order
 * \returns The next peptide in the database.
 */
PEPTIDE_T* database_sorted_peptide_iterator_next(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_peptide_iterator ///< the iterator of interest -in
  );

/**********************************************************************
 * wrapper, for generate_peptides_iterator, cast back to original type
 ***********************************************************************/

/**
 * Frees an allocated database_sorted_peptide_iterator object.
 */
void void_free_database_sorted_peptide_iterator(
  void* database_peptide_iterator ///< the iterator to free -in
  );

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T void_database_sorted_peptide_iterator_has_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  );

/**
 * returns each peptide in sorted order
 * \returns The next peptide in the database.
 */
PEPTIDE_T* void_database_sorted_peptide_iterator_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  );

#ifdef __cplusplus
}
#endif


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
