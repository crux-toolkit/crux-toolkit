/*************************************************************************//**
 * \file DatabaseSortedPeptideIterator.cpp
 * \brief Object to iterator over the peptides within a database, in an
 * specified sorted order. (mass, length, lexical)
 ****************************************************************************/

#include "DatabaseSortedPeptideIterator.h"
#include "DatabasePeptideIterator.h"

/***********************************
 * database sorted peptide iterator
 ***********************************/

/**
 * Instantiates a new database_sorted_peptide_iterator from a database.
 * uses a sorted_peptide_iterator as it's engine
 * \returns a DATABASE_SORTED_PEPTIDE_ITERATOR_T object.
 */
DatabaseSortedPeptideIterator::DatabaseSortedPeptideIterator(
  Database* database, ///< the database of interest -in
  PeptideConstraint* peptide_constraint, 
    ///< the peptide_constraint to filter peptides -in
  SORT_TYPE_T sort_type, ///< the sort type for this iterator -in
  bool unique ///< only return unique peptides? -in
  )
{
  // initialize
  sorted_peptide_iterator_ = NULL;

  // create the database peptide iterator
  DatabasePeptideIterator* db_peptide_iterator =
    new DatabasePeptideIterator(database, peptide_constraint, 
                                  false);// do not store peptides

  // create a sorted peptide iterator from db peptide iterator
  SORTED_PEPTIDE_ITERATOR_T* sorted_peptide_iterator =
    new_sorted_peptide_iterator_database(
        db_peptide_iterator, sort_type, unique);

  // set sorted_peptide_iterator
  sorted_peptide_iterator_ = sorted_peptide_iterator;
  
  delete db_peptide_iterator; 
}

/**
 * The basic iterator functions.
 * \returns true if there are additional peptides, false if not.
 */
bool DatabaseSortedPeptideIterator::hasNext() {

  return sorted_peptide_iterator_has_next(sorted_peptide_iterator_);
}

/**
 * returns each peptide in sorted order
 * \returns The next peptide in the database.
 */
Peptide* DatabaseSortedPeptideIterator::next() {

  return sorted_peptide_iterator_next(sorted_peptide_iterator_);
}

/**
 * Frees an allocated database_sorted_peptide_iterator object.
 */
DatabaseSortedPeptideIterator::~DatabaseSortedPeptideIterator() {

  free_sorted_peptide_iterator(sorted_peptide_iterator_);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
