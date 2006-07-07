/**
 * \file database.h 
 * $Revision: 1.1 $
 * \brief Object for representing a database of protein sequences.
 *****************************************************************************/
#ifndef DATABASE_H
#define DATABASE_H

#include <stdio.h>
#include "utils.h"
#include "peptide.h"
#include "protein.h"

/* CHRIS This is probably an object for which you can crib code for from an outside source. Even from in-house (like Charles).*/

/**
 * \typedef DATABASE_T
 * \brief A database of protein sequences.
 */
typedef struct database DATABASE_T;

/**
 * \typedef DATABASE_PROTEIN_ITERATOR_T
 * \brief An object to iterate over the proteins in a database 
 */
typedef struct database_protein_iterator DATABASE_PROTEIN_ITERATOR_T;

/**
 * \typedef DATABASE_PEPTIDE_ITERATOR_T
 * \brief An object to iterate over the peptides in a database 
 */
typedef struct database_peptide_iterator DATABASE_PEPTIDE_ITERATOR_T;

/**
 * \returns An (empty) database object.
 */
DATABASE_T* allocate_database(void);

/**
 * \returns A new database object.
 */
DATABASE_T* new_database(
  char*         filename ///< The file from which to parse the database.
  );         

/**
 * Frees an allocated protein object.
 */
void free_database(DATABASE_T* protein);

/**
 * Prints a database object to file.
 */
void print_database(DATABASE_T* database, FILE* file);

// START HERE

/**
 * Copies database object src to dest.
 */
void copy_database(
  DATABASE_T* src,
  DATABASE_T* dest);

/**
 * Parses a database from the file in the filename member variable
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_database(
  DATABASE_T* database ///< An allocated database
  );

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 * FIXME Chris, could you create the get and set methods for the object fields?
 */

/**
 * Additional get and set methods
 */

/*int get_protein_peptides(PROTEIN_T* protein);*/
/*int get_protein_peptides(PROTEIN_T* protein);*/

/***********************************************
 * Iterators
 ***********************************************/

/**
 * Instantiates a new database_protein_iterator from a database.
 * \returns a DATABASE_PROTEIN_ITERATOR_T object.
 */
DATABASE_PROTEIN_ITERATOR_T* new_database_protein_iterator(
    DATABASE_T* protein);        

/**
 * Frees an allocated database_protein_iterator object.
 */
void free_database_protein_iterator(
    DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator);

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional proteins to iterate over, FALSE if not.
 */
BOOLEAN_T database_protein_iterator_has_next(
    DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator);

/**
 * \returns The next protein in the database.
 */
PROTEIN_T* database_protein_iterator_next(
    DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator);

/**
 * Instantiates a new database_peptide_iterator from a database.
 * \returns a DATABASE_PEPTIDE_ITERATOR_T object.
 */
DATABASE_PEPTIDE_ITERATOR_T* new_database_peptide_iterator(
    DATABASE_T* protein,
    PEPTIDE_CONSTRAINT_T* peptide_constraint);

/**
 * Frees an allocated database_peptide_iterator object.
 */
void free_database_peptide_iterator(
    DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator);

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T database_peptide_iterator_has_next(
    DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator);

/**
 * \returns The next peptide in the database.
 */
PEPTIDE_T* database_peptide_iterator_next(
    DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator);


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
