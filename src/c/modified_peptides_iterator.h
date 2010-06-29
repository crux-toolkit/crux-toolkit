/**
 * \file modified_peptides_iterator.h
 * AUTHOR: Barbara Frewen
 * DATE: April 15, 2008
 * DESCRIPTION: Header file for peptide iterator that includes
 * modified peptides.
 * $Revision: 1.5 $
 */
#ifndef MODIFIED_PEPTIDES_ITERATOR_H
#define MODIFIED_PEPTIDES_ITERATOR_H

#include "utils.h"
#include "objects.h"
#include "linked_list.h"
#include "generate_peptides_iterator.h"

/**
 * \brief Create a new modified_PEPTIDES_iterator for a specific mass.
 *
 * The returned iterator is initialized with the first peptide queued
 * up and ready to return.  Also creates a generate_peptides_iterator
 * from which it gets the peptides to modify. Peptides are of mass +/-
 * mass-window taken from parameter.c.  All other peptide
 * specifications are taken from parameter.c.  If no peptides meet the
 * specifications, an iterator is still returned and when given to
 * has_next will always return FALSE.
 * 
 * \returns A newly allocated modified_peptides_iterator.
 */
MODIFIED_PEPTIDES_ITERATOR_T* new_modified_peptides_iterator_from_mass(
  double mass,         ///< Target mass of peptides
  PEPTIDE_MOD_T* pmod, ///< Peptide mod to apply
  BOOLEAN_T is_decoy,  ///< generate decoy peptides
  INDEX_T* index,      ///< Index from which to draw peptides OR
  DATABASE_T* dbase    ///< Database from which to draw peptides
  );

/**
 * \brief Create a new modified_PEPTIDES_iterator for a specific mass/charge.
 *
 * The returned iterator is initialized with the first peptide queued
 * up and ready to return.  Also creates a generate_peptides_iterator
 * from which it gets the peptides to modify. Peptides are of mass +/-
 * mass-window or mz +/- mass-window(mz) taken from parameter.c.  All other peptide
 * specifications are taken from parameter.c.  If no peptides meet the
 * specifications, an iterator is still returned and when given to
 * has_next will always return FALSE.
 * 
 * \returns A newly allocated modified_peptides_iterator.
 */
MODIFIED_PEPTIDES_ITERATOR_T* new_modified_peptides_iterator_from_mz(
  double mz,         ///< Target mz of peptides
  int charge,        ///< Charge of peptides
  PEPTIDE_MOD_T* pmod, ///< Peptide mod to apply
  BOOLEAN_T is_decoy,  ///< generate decoy peptides
  INDEX_T* index,      ///< Index from which to draw peptides OR
  DATABASE_T* dbase    ///< Database from which to draw peptides
  );

/**
 * \brief Create a new modified_PEPTIDES_iterator for all peptides in
 * the database or index.
 *
 * The returned iterator is initialized with the first peptide queued
 * up and ready to return.  Also creates a generate_peptides_iterator
 * from which it gets the peptides to modify. Peptides range from mass
 * 'min mass' + pmod->delta_mass to 'max mass' + pmod->delta_mass (min
 * and max taken from parameter.c).  All other peptide
 * specifications are taken from parameter.c.  If no peptides meet the
 * specifications, an iterator is still returned and when passed to
 * has_next() will always return FALSE.
 * 
 * \returns A newly allocated modified_peptides_iterator.
 */
MODIFIED_PEPTIDES_ITERATOR_T* new_modified_peptides_iterator(
  PEPTIDE_MOD_T* pmod, ///< Peptide mod to apply
  INDEX_T* index,      ///< Index from which to draw peptides OR
  DATABASE_T* dbase    ///< Database from which to draw peptides
  );

/**
 * \brief Check to see if the iterator has more peptides to return,
 * i.e. a call to has_next will return non-NULL.
 * \returns TRUE if the iterator has more peptides to return.
 */
BOOLEAN_T modified_peptides_iterator_has_next(
  MODIFIED_PEPTIDES_ITERATOR_T* modified_peptides_iterator);

/**
 * \brief Return the next peptide or NULL if no peptides remain.
 * \returns A modified peptide.
 */
PEPTIDE_T* modified_peptides_iterator_next(
  MODIFIED_PEPTIDES_ITERATOR_T* modified_peptides_iterator);

/**
 * \brief Free the memory used by this iterator.
 */
void free_modified_peptides_iterator(
  MODIFIED_PEPTIDES_ITERATOR_T* modified_peptides_iterator);

/**
 * \brief Void wrapper of modified_peptides_iterator_has_next to be
 * used by generate_peptides_iterator. 
 */
BOOLEAN_T void_modified_peptides_iterator_has_next(
  void* modified_peptides_iterator);

/**
 * \brief Void wrapper of modified_peptides_iterator_next to be
 * used by generate_peptides_iterator. 
 */
PEPTIDE_T* void_modified_peptides_iterator_next(
  void* modified_peptides_iterator);

/**
 * \brief Void wrapper of modified_peptides_iterator_free to be
 * used by generate_PEPTIDES_iterator. 
 */
void void_modified_peptides_iterator_free(
  void* modified_PEPTIDES_iterator);


#endif // MODIFIED_PEPTIDES_ITERATOR_H
