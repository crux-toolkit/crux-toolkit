/**
 * \file peptide_src.h
 * $Revision: 1.4 $
 * \brief Object for mapping a peptide to it's parent protein.
 */
#ifndef PEPTIDE_SRC_H
#define PEPTIDE_SRC_H

#include "utils.h"
#include "mass.h"
#include <stdio.h>
#include "objects.h"
#include "carp.h"

/**
 * \returns An (empty) peptide_src object.
 */
PEPTIDE_SRC_T* allocate_peptide_src(void);

/**
 *\returns a PROTEIN_PEPTIDE_ASSOCIATION object, populated with user specified parameters
 */
PEPTIDE_SRC_T* new_peptide_src(
  PEPTIDE_TYPE_T peptide_type, ///< the peptide type for the corresponding protein -in
  PROTEIN_T* parent_protein, ///< the parent of this preptide -in
  int start_idx ///< start index of the peptide in the protein sequence, first residue is 1 -in
  );


/**
 * Frees the entire allocated peptide_src linklist object
 */
void free_peptide_src(
  PEPTIDE_SRC_T* peptide_src ///< object to free -in 
  );

/**
 * Frees the an individual allocated peptide_src object
 * assumes that new_association pointer is NULL or some other pointer exist for the rest of the linklist 
 */
void free_one_peptide_src(
  PEPTIDE_SRC_T* peptide_src///< object to free -in 
  );

/**
 * Prints a peptide object to file.
 */
void print_peptide_src(
  PEPTIDE_SRC_T* peptide_src, ///< object to print -in 
  FILE* file  ///< the out put stream -out
  );

/**
 * Copies the entire linklist of peptide_src object src to dest.
 * dest must be a heap allocated peptide_src
 */
void copy_peptide_src(
  PEPTIDE_SRC_T* src, ///< source peptide_src -in
  PEPTIDE_SRC_T* dest ///< destination peptide_src -out
  );

/**
 * sets the peptide type
 * peptide type: TRYPTIC, PARTIALLY_TRYPTIC, NON_TRYPTIC
 */
void set_peptide_src_peptide_type( 
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  PEPTIDE_TYPE_T peptide_type ///< the type of the peptide -in
  );

/**
 * \returns the peptide type with association to the parent protein
 * peptide type: TRYPTIC, PARTIALLY_TRYPTIC, NON_TRYPTIC
 */
PEPTIDE_TYPE_T get_peptide_src_peptide_type( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  );

/**
 * sets the parent protein
 */
void set_peptide_src_parent_protein(
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  PROTEIN_T* parent_protein ///< the parent of this preptide -in  
  );

/**
 * \returns a pointer to the parent protein
 */
PROTEIN_T* get_peptide_src_parent_protein( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  );

/**
 * sets the start index of the peptide in the protein sequence
 */
void set_peptide_src_start_idx(
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  int start_idx ///< start index of the peptide in the protein sequence -in
  );

/**
 * \returns the start index of the peptide in the protein sequence
 */
int get_peptide_src_start_idx( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  );

/**
 * sets the next peptide_src on the link list
 * assumes that the src_association's next_association feild is NULL
 */
void set_peptide_src_next_association(
  PEPTIDE_SRC_T* src_association, ///< the peptide_src to set -out   
  PEPTIDE_SRC_T* new_association ///< the new peptide_src to add -in   
  );

/**
 * \returns the next peptide_src on the link list
 */
PEPTIDE_SRC_T* get_peptide_src_next_association( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  );

/**
 * \returns a pointer to the start of the peptide with in it's parent protein sequence
 */
char* get_peptide_src_sequence_pointer(
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
