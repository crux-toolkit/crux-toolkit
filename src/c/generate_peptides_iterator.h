/**
 * \file generate_peptides_iterator.h 
 * $Revision: 1.5 $
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
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator_sp(
  float neutral_mass ///< the neutral_mass that which the peptides will be searched -in
);

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

#endif
