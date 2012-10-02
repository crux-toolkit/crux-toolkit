/**
 * \file sorter.h
 * $Revision: 1.4 $
 * \brief Object to sort objects
 ****************************************************************************/
#ifndef SORTER_H
#define SORTER_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#ifndef _MSC_VER
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <time.h>

#include "utils.h"
#include "crux-utils.h"
#include "Peptide.h"
#include "Protein.h"
#include "Index.h"
#include "carp.h"
#include "objects.h"
#include "PeptideConstraint.h"
#include "Database.h"


/***********************************
 * sorted peptide iterator
 ***********************************/


/**
 * Instantiates a new sorted_peptide_iterator from a bin_peptide_iterator
 * \returns a SORTED_PEPTIDE_ITERATOR_T object.
 */
SORTED_PEPTIDE_ITERATOR_T* new_sorted_peptide_iterator_bin(
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator, ///< the peptide iterator to extend -in
  SORT_TYPE_T sort_type, ///< the sort type for this iterator -in
  bool unique, ///< only return unique peptides? -in
  unsigned int peptide_count ///< the total peptide count in the bin -in
  );

/**
 * Frees an allocated database_sorted_peptide_iterator object.
 */
void free_sorted_peptide_iterator(
  SORTED_PEPTIDE_ITERATOR_T* database_peptide_iterator ///< the iterator to free -in
  );

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
bool sorted_peptide_iterator_has_next(
  SORTED_PEPTIDE_ITERATOR_T* peptide_iterator ///< the iterator of interest -in
  );

/**
 * returns each peptide in sorted order
 * \returns The next peptide in the database.
 */
Crux::Peptide* sorted_peptide_iterator_next(
  SORTED_PEPTIDE_ITERATOR_T* peptide_iterator ///< the iterator of interest -in
  );



#endif // SORTER_H
