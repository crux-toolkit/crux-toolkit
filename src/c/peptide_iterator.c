/*****************************************************************************
 * \file peptide.c
 * $Revision: 1.1 $
 * \brief: Object for iterating over peptides, either a single protein.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"
#include "peptide_iterator.h"

/**
 * \struct peptide_iterator
 */
struct peptide_iterator {
  PROTEIN_T* protein; ///< The spectrum whose peaks to iterate over. 
  unsigned short int cur_start; ///< Start in protein of the current peptide.
  unsigned char cur_length; ///< The length of the current peptide.
};

/* FIXME need to allow for iteration over peptides in proteins and in 
 * a larger database.*/
  
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

