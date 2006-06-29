/**
 * \file peptide_iterator.h 
 * $Revision: 1.2 $
 * \brief Object for iterator over peptides.
 *****************************************************************************/
#ifndef PEPTIDE_ITERATOR_H
#define PEPTIDE_ITERATOR_H

#include <stdio.h>
#include "utils.h"
#include "peptide.h"
#include "protein.h"

/**
 * Instantiates a new peptide_iterator from a peptide.
 * \returns a PEPTIDE_ITERATOR_T object.
 */

/**
 * \typedef PEPTIDE_ITERATOR_T
 */
typedef struct peptide_iterator PEPTIDE_ITERATOR_T;


PEPTIDE_ITERATOR_T* new_peptide_iterator(PROTEIN_T* protein);        

/**
 * Frees an allocated peptide_iterator object.
 */
void free_peptide_iterator(PEPTIDE_ITERATOR_T* peptide_iterator);

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T peptide_iterator_has_next(PEPTIDE_ITERATOR_T* peptide_iterator);

/**
 * \returns The next peptide (a character) in the peptide.
 */
PEPTIDE_T* peptide_iterator_next(PEPTIDE_ITERATOR_T* peptide_iterator);
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
