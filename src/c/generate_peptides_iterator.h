/*****************************************************************************
 * \file generate_peptide_iterator.h 
 * AUTHOR: Chris Park
 * CREATE DATE: 8 Nov 2007
 * $Revision: 1.1 $
 * DESCRIPTION: object to return candidate peptides
 *****************************************************************************/
#ifndef GENERATE_PEPTIDE_ITERATOR_H 
#define GENERATE_PEPTIDE_ITERATOR_H 

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
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
GENERATE_PEPTIDE_ITERATOR_T* allocate_generate_peptide_iterator(void);

/**
 *\returns a new generate_peptide_iterator object
 */
GENERATE_PEPTIDE_ITERATOR_T* new_generate_peptide_iterator(void);

/**
 *\returns TRUE, if there is a next peptide, else FALSE
 */
BOOLEAN_T generate_peptide_iterator_has_next(
  GENERATE_PEPTIDE_ITERATOR_T* generate_peptide_iterator ///< working iterator
  );

/**
 *\returns the next peptide in the iterator
 */
PEPTIDE_T* generate_peptide_iterator_next(
  GENERATE_PEPTIDE_ITERATOR_T* generate_peptide_iterator ///< working iterator
  );

/**
 * Frees an allocated generate_peptide_iterator object
 */
void free_generate_peptide_iterator(
  GENERATE_PEPTIDE_ITERATOR_T* generate_peptide_iterator ///< iterator to free
  );

