/**
 * \file generate_peptides_iterator.h 
 * $Revision: 1.10 $
 * \brief object to return candidate peptides from database
 *****************************************************************************/
#ifndef GENERATE_PEPTIDES_ITERATOR_H 
#define GENERATE_PEPTIDES_ITERATOR_H 

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <dirent.h>
#include "utils.h"
#include "crux-utils.h"
#include "peptide.h"
#include "protein.h"
#include "index.h"
#include "carp.h"
#include "objects.h"
#include "peptide_constraint.h"
#include "database.h"


/**
 *\returns a empty generate_peptide_iterator object
 */
GENERATE_PEPTIDES_ITERATOR_T* allocate_generate_peptides_iterator(void);

/**
 *\returns a new generate_peptide_iterator object
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator(void);

/**
 *\returns a new generate_peptide_iterator object
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator_from_mass(
  float neutral_mass, ///< the neutral_mass that which the peptides will be searched -in
  INDEX_T* index,
  DATABASE_T* database
);

/**
 * Used for when need to resue genearte peptide iterator mutiple times
 * only changing by the mass window
 * MUST use it with index
 * MUST call set_generate_peptides_mutable before using this iterator
 *\returns a new generate_peptides_iterator object that can is mutable
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator_mutable(void);

/**
 * sets the mutable generate_peptides_iterator for the next iteration of creating peptides
 * user provides the mass window the iterator will operate
 */
void set_generate_peptides_mutable(
  GENERATE_PEPTIDES_ITERATOR_T* gen_peptides_iterator, ///< working mutable iterator
  float max_mass, ///< the max mass which the peptides will be searched -in
  float min_mass ///< the min mass that which the peptides will be searched -in
  );

/***********************************************************/

/**
 *\returns TRUE, if there is a next peptide, else FALSE
 */
BOOLEAN_T generate_peptides_iterator_has_next(
  GENERATE_PEPTIDES_ITERATOR_T* generate_peptide_iterator ///< working iterator
  );

/**
 *\returns the next peptide in the iterator
 */
PEPTIDE_T* generate_peptides_iterator_next(
  GENERATE_PEPTIDES_ITERATOR_T* generate_peptide_iterator ///< working iterator
  );

/**
 * Don't free the iterator until completed with the peptides generated
 * Frees an allocated generate_peptide_iterator object
 */
void free_generate_peptides_iterator(
  GENERATE_PEPTIDES_ITERATOR_T* generate_peptide_iterator ///< iterator to free
  );
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
