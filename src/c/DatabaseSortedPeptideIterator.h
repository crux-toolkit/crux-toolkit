/**
 * \file DatabaseSortedPeptideIterator.h 
 * $Revision: 1.27 $
 * \brief Object to iterator over the peptides within a database, in an
 * specified sorted order. (mass, length, lexical)
 *****************************************************************************/

#ifndef DATABASESORTEDPEPTIDEITERATOR_H
#define DATABASESORTEDPEPTIDEITERATOR_H

#include "objects.h"

/**
 * \class DatabaseSortedPeptideIterator
 * \brief Object to iterate over the peptides within a database, in an
 * specified sorted order.(mass, length, lexical)
 */
class DatabaseSortedPeptideIterator {
 protected:
  SORTED_PEPTIDE_ITERATOR_T* sorted_peptide_iterator_; ///< the peptide iterator that sorts the peptides
 public:

  /**
   * Instantiates a new database_sorted_peptide_iterator from a database.
   * uses a sorted_peptide_iterator as it's engine
   * \returns a DATABASE_SORTED_PEPTIDE_ITERATOR_T object.
   */
  DatabaseSortedPeptideIterator(
    Database* database, ///< the database of interest -in
    PeptideConstraint* peptide_constraint, ///< the peptide_constraint to filter peptides -in
    SORT_TYPE_T sort_type, ///< the sort type for this iterator
    bool unique ///< only return unique peptides? -in
    );

  /**
   * Frees an allocated database_sorted_peptide_iterator object.
   */
  virtual ~DatabaseSortedPeptideIterator();

  /**
   * The basic iterator functions.
   * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
   */
  bool hasNext();

  /**
   * returns each peptide in sorted order
   * \returns The next peptide in the database.
   */
  Peptide* next();

};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
