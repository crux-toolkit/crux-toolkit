/**
 * \file generate_peptides_iterator.h 
 * $Revision: 1.14 $
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
#include "Protein.h"
#include "index.h"
#include "carp.h"
#include "objects.h"
#include "PeptideConstraint.h"
#include "peptide_modifications.h"
#include "Database.h"
#include "linked_list.h"
#include "modified_peptides_iterator.h"

// TODO (BF 10-Apr-08) should this be private?
/**
 *\returns a empty generate_peptide_iterator object
 */
GENERATE_PEPTIDES_ITERATOR_T* allocate_generate_peptides_iterator(void);

/**
 * \brief Create a peptide iterator based entirely on parameter values
 * (defaults and those given by user).
 *
 * The peptides are drawn from either a fasta file or an index based
 * on the value of the use-index parameter.  With default values, mass
 * range is wide so this effectively iterates over "all" peptides in
 * the protein input.  If no peptides in the protein source meet the
 * criteria, a peptide iterator is still returned, but when passed to
 * *_has_next() it will always return FALSE.
 *\returns A new generate_peptide_iterator object
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator(void);

/**
 * \brief Create a peptide iterator for peptides of a target mass.
 *
 * Peptides with mass between neutral_mass +/- "mass-window"
 * (a user-defined parameter).  Peptides are drawn from either the
 * given index or the given database, if the index is NULL. If no
 * peptides in the protein source meet the criteria, a peptide
 * iterator is still returned, but when passed to *_has_next() it will
 * always return FALSE.  
 *
 *\returns A new generate_peptide_iterator object
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator_from_mass(
  FLOAT_T neutral_mass, ///< The target mass (uncharged) for peptides
  INDEX_T* index,     ///< The index from which to draw peptides OR
  Database* database///< The database from which to draw peptides
);

/**
 * \brief Create a peptide iterator for peptides in a specific range
 * of masses.
 *
 * This is the version of new_* that is called by the others and does
 * all the work of setting the member variable fields.  Parameters
 * other than min and max mass are taken from parameter.c.  If no
 * peptides in the protein source meet the criteria, a peptide
 * iterator is still returned, but when passed to *_has_next() it will
 * always return FALSE.  
 *
 * \returns a new generate_peptides_iterator object, with fasta file input
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator_from_mass_range(
  double min_mass,     ///< The min mass of peptides to generate -in
  double max_mass,     ///< The maximum mas of peptide to generate -in
  INDEX_T* index,      ///< The index
  Database* database ///< The database
  );

/**
 * \brief Create a new peptide iterator to return modified peptides.
 *
 * Peptides are generated based on values in parameter.c with the
 * exception of the mass range.  Min and max mass are given by target
 * mass +/- mass-window.  Only those peptides that can be modified by
 * the peptide_mod are returned.  (A peptide_mod may contain no
 * AA_MODS, in which case all peptides are returned.)  Peptides are
 * taken either from an index or from a database (fasta file).  If no
 * peptides pass the criteria, a new iterator is still returned, but
 * when passed to has_next() it will always return FALSE.
 * \returns A newly allocated peptide iterator.
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator_mods(
  double mass,                ///< target mass of peptides
  PEPTIDE_MOD_T* peptide_mod, ///< the peptide mod to apply
  INDEX_T* index,             ///< index from which to draw peptides OR
  Database* dbase           ///< database from which to draw peptides
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
