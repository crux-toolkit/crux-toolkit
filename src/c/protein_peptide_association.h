/*****************************************************************************
 * \file protein_peptide_association.h
 * $Revision: 1.3 $
 * \brief: Object for mapping a peptide to it's parent protein.
 ****************************************************************************/
#ifndef PROTEIN_PEPTIDE_ASSOCIATION_H
#define PROTEIN_PEPTIDE_ASSOCIATION_H

#include "utils.h"
#include "mass.h"
#include <stdio.h>
#include "objects.h"

/**
 * \returns An (empty) protein_peptide_association object.
 */
PROTEIN_PEPTIDE_ASSOCIATION_T* allocate_protein_peptide_association(void);

/**
 *\returns a PROTEIN_PEPTIDE_ASSOCIATION object, populated with user specified parameters
 */
PROTEIN_PEPTIDE_ASSOCIATION_T* new_protein_peptide_association(
  PEPTIDE_TYPE_T peptide_type, ///< the peptide type for the corresponding protein -in
  PROTEIN_T* parent_protein, ///< the parent of this preptide -in
  int start_idx ///< start index of the peptide in the protein sequence -in
  );


/**
 * Frees the entire allocated protein_peptide_association linklist object
 */
void free_protein_peptide_association(
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association ///< object to free -in 
  );

/**
 * Frees the an individual allocated protein_peptide_association object
 * assumes that new_association pointer is NULL or some other pointer exist for the rest of the linklist 
 */
void free_one_protein_peptide_association(
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association///< object to free -in 
  );

/**
 * Prints a peptide object to file.
 */
void print_protein_peptide_association(
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association, ///< object to print -in 
  FILE* file  ///< the out put stream -out
  );

/**
 * Copies the entire linklist of protein_peptide_association object src to dest.
 * dest must be a heap allocated protein_peptide_association
 */
void copy_protein_peptide_association(
  PROTEIN_PEPTIDE_ASSOCIATION_T* src, ///< source protein_peptide_association -in
  PROTEIN_PEPTIDE_ASSOCIATION_T* dest ///< destination protein_peptide_association -out
  );

/**
 * sets the peptide type
 * peptide type: TRYPTIC, PARTIALLY_TRYPTIC, NON_TRYPTIC
 */
void set_protein_peptide_association_peptide_type( 
  PROTEIN_PEPTIDE_ASSOCIATION_T* new_association, ///< the protein_peptide_association to set -out   
  PEPTIDE_TYPE_T peptide_type ///< the type of the peptide -in
  );

/**
 * \returns the peptide type with association to the parent protein
 * peptide type: TRYPTIC, PARTIALLY_TRYPTIC, NON_TRYPTIC
 */
PEPTIDE_TYPE_T get_protein_peptide_association_peptide_type( 
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association ///< the query protein_peptide_association -in   
  );

/**
 * sets the parent protein
 */
void set_protein_peptide_association_parent_protein(
  PROTEIN_PEPTIDE_ASSOCIATION_T* new_association, ///< the protein_peptide_association to set -out   
  PROTEIN_T* parent_protein ///< the parent of this preptide -in  
  );

/**
 * \returns a pointer to the parent protein
 */
PROTEIN_T* get_protein_peptide_association_parent_protein( 
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association ///< the query protein_peptide_association -in   
  );

/**
 * sets the start index of the peptide in the protein sequence
 */
void set_protein_peptide_association_start_idx(
  PROTEIN_PEPTIDE_ASSOCIATION_T* new_association, ///< the protein_peptide_association to set -out   
  int start_idx ///< start index of the peptide in the protein sequence -in
  );

/**
 * \returns the start index of the peptide in the protein sequence
 */
int get_protein_peptide_association_start_idx( 
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association ///< the query protein_peptide_association -in   
  );

/**
 * sets the next protein_peptide_association on the link list
 */
void set_protein_peptide_association_next_association(
  PROTEIN_PEPTIDE_ASSOCIATION_T* src_association, ///< the protein_peptide_association to set -out   
  PROTEIN_PEPTIDE_ASSOCIATION_T* new_association ///< the new protein_peptide_association to add -in   
  );

/**
 * \returns the next protein_peptide_association on the link list
 */
PROTEIN_PEPTIDE_ASSOCIATION_T* get_protein_peptide_association_next_association( 
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association ///< the query protein_peptide_association -in   
  );



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
