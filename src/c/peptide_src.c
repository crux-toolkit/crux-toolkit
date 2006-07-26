/*****************************************************************************
 * \file peptide_src.c
 * $Revision: 1.3 $
 * \brief: Object for mapping a peptide to it's parent protein.
 ****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "carp.h"
#include "utils.h"
#include "mass.h"
#include "objects.h"
#include "peptide.h"
#include "protein.h"
#include "peptide_src.h"


/**
 * \struct peptide_src
 * \brief object for mapping a peptide to it's parent protein.
 */
struct peptide_src{
  PEPTIDE_TYPE_T peptide_type; ///< the peptide type for the corresponding protein
  PROTEIN_T* parent_protein; ///< the parent of this preptide
  int start_idx; ///< start index of the peptide in the protein sequence, first residue is 1 
  PEPTIDE_SRC_T* next_association; ///< a linklist of peptide_src   
};


/**
 * \returns An (empty) peptide_src object.
 */
PEPTIDE_SRC_T* allocate_peptide_src(void){
  PEPTIDE_SRC_T* peptide_src =
    mymalloc(sizeof(PEPTIDE_SRC_T));

  peptide_src->next_association = NULL;
  return peptide_src;
}

/**
 *\returns a PROTEIN_PEPTIDE_ASSOCIATION object, populated with user specified parameters
 */
PEPTIDE_SRC_T* new_peptide_src(
  PEPTIDE_TYPE_T peptide_type, ///< the peptide type for the corresponding protein -in
  PROTEIN_T* parent_protein, ///< the parent of this preptide -in
  int start_idx ///< start index of the peptide in the protein sequence -in
  )
{
  PEPTIDE_SRC_T* new_association = allocate_peptide_src();
  set_peptide_src_peptide_type(new_association, peptide_type);
  set_peptide_src_parent_protein(new_association, parent_protein);
  set_peptide_src_start_idx(new_association, start_idx);
  return new_association;
}

/**
 * Frees the entire allocated peptide_src linklist object
 */
void free_peptide_src(
  PEPTIDE_SRC_T* peptide_src  ///< object to free -in 
  )
{
  if(peptide_src->next_association != NULL){
    free_peptide_src(peptide_src->next_association);
  }    
  free(peptide_src);
}

/**
 * Frees the an individual allocated peptide_src object
 * assumes that new_association pointer is NULL or some other pointer exist for the rest of the linklist 
 */
void free_one_peptide_src(
  PEPTIDE_SRC_T* peptide_src  ///< object to free -in 
  )
{
  free(peptide_src);
}

//FIXME might need to change how this is printed
/**
 * Prints a peptide object to file.
 */
void print_peptide_src(
  PEPTIDE_SRC_T* peptide_src, ///< object to print -in 
  FILE* file  ///< the out put stream -out
  )
{
  PEPTIDE_TYPE_T peptide_type = get_peptide_src_peptide_type(peptide_src);
  char* sequence = get_protein_sequence(peptide_src->parent_protein);
  fprintf(file, "parent protein:%s\n", sequence);

  fprintf(file, "peptide start: %d\n", peptide_src->start_idx);

  if(peptide_type == TRYPTIC){
    fprintf(file, "peptide type:%s\n", "TRYPTIC");
  }
  else if(peptide_type == PARTIALLY_TRYPTIC){
    fprintf(file, "peptide type:%s\n", "PARTIALLY_TRYPTIC");
  }
  else if(peptide_type == NOT_TRYPTIC){
    fprintf(file, "peptide type:%s\n", "NOT_TRYPTIC");
  }
  else if(peptide_type == ANY_TRYPTIC){
    fprintf(file, "peptide type:%s\n", "ANY_TRYPTIC");
  }
  free(sequence);
}

/**
 * Copies the entire linklist of peptide_src object src to dest.
 * dest must be a heap allocated peptide_src
 */
void copy_peptide_src(
  PEPTIDE_SRC_T* src, ///< source peptide_src -in
  PEPTIDE_SRC_T* dest ///< destination peptide_src -out
  )
{
  PEPTIDE_SRC_T* next_association;
  set_peptide_src_peptide_type(dest, src->peptide_type);
  set_peptide_src_parent_protein(dest, src->parent_protein);
  set_peptide_src_start_idx(dest, src->start_idx);
  //check if end of the linklist
  if(get_peptide_src_next_association(src) != NULL){
    next_association = allocate_peptide_src();
    dest->next_association = next_association;
    copy_peptide_src(src->next_association, next_association);
  }
}

/**
 * sets the peptide type
 * peptide type: TRYPTIC, PARTIALLY_TRYPTIC, NON_TRYPTIC
 */
void set_peptide_src_peptide_type( 
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  PEPTIDE_TYPE_T peptide_type ///< the type of the peptide -in
  )
{
  new_association->peptide_type = peptide_type;
}

/**
 * \returns the peptide type with association to the parent protein
 * peptide type: TRYPTIC, PARTIALLY_TRYPTIC, NON_TRYPTIC
 */
PEPTIDE_TYPE_T get_peptide_src_peptide_type( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  )
{
  return peptide_src->peptide_type;
}


/**
 * sets the parent protein
 */
void set_peptide_src_parent_protein(
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  PROTEIN_T* parent_protein ///< the parent of this preptide -in  
  )
{
  new_association->parent_protein = parent_protein;

}

/**
 * \returns a pointer to the parent protein
 */
PROTEIN_T* get_peptide_src_parent_protein( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  )
{
  return peptide_src->parent_protein;
}

/**
 * sets the start index of the peptide in the protein sequence
 */
void set_peptide_src_start_idx(
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  int start_idx ///< start index of the peptide in the protein sequence -in
  )
{
  new_association->start_idx = start_idx;
}

/**
 * \returns the start index of the peptide in the protein sequence
 */
int get_peptide_src_start_idx( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  )
{
  return peptide_src->start_idx;
}

/**
 * sets the next peptide_src on the link list
 * assumes that the src_association's next_association feild is NULL
 */
void set_peptide_src_next_association(
  PEPTIDE_SRC_T* src_association, ///< the peptide_src to set -out   
  PEPTIDE_SRC_T* new_association ///< the new peptide_src to add -in   
  )
{
  src_association->next_association = new_association;
}

/**
 * \returns the next peptide_src on the link list
 */
PEPTIDE_SRC_T* get_peptide_src_next_association( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  )
{
  return peptide_src->next_association;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
