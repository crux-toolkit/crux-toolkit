/*****************************************************************************
 * \file database.c
 * $Revision: 1.1 $
 * \brief: Object for representing a database of protein sequences.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"
#include "peptide.h"
#include "protein.h"
#include "database.h"

/**
 * \struct database
 * \brief A database of protein sequences
 */
struct database{
  char*          filename; ///< The original database filename.
  PROTEIN_T**    proteins; ///< The proteins in this database.
};    

/**
 * \struct database_protein_iterator
 * \brief Object to iterate over the proteins within a database
 */
struct database_protein_iterator {
  DATABASE_T* database;  ///< The database whose proteins to iterate over 
  int cur_protein;      ///< The index of the current protein
};

/**
 * \struct database_peptide_iterator
 * \brief Object to iterate over the peptides within a database, in an
 * unspecified order.
 */
struct database_peptide_iterator {
  DATABASE_PROTEIN_ITERATOR_T* 
    database_protein_iterator; ///<The protein iterator. 
  PROTEIN_PEPTIDE_ITERATOR_T* 
    cur_protein_peptide_iterator; ///< The peptide iterator for the current protein.
  PEPTIDE_CONSTRAINT_T* peptide_constraint; ///< The constraints for the kind of peptide to iterate over.
  };

// FIXME I think all of these fields are necessary? But maybe some can go
// in the PROTEIN_PEPTIDE_ITERATOR?

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

