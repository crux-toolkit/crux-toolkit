/*****************************************************************************
 * \file protein_peptide_association.c
 * $Revision: 1.7 $
 * \brief: Object for mapping a peptide to it's parent protein.
 ****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "mass.h"
#include "objects.h"
#include "peptide.h"
#include "protein.h"
#include "protein_peptide_association.h"


/**
 * \struct protein_peptide_association
 * \brief object for mapping a peptide to it's parent protein.
 */
struct protein_peptide_association{
  PEPTIDE_TYPE_T peptide_type; ///< the peptide type for the corresponding protein
  PROTEIN_T* parent_protein; ///< the parent of this preptide
  int start_idx; ///< start index of the peptide in the protein sequence, first residue is 1 
  PROTEIN_PEPTIDE_ASSOCIATION_T* next_association; ///< a linklist of protein_peptide_association   
};


/**
 * \returns An (empty) protein_peptide_association object.
 */
PROTEIN_PEPTIDE_ASSOCIATION_T* allocate_protein_peptide_association(void){
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association =
    mymalloc(sizeof(PROTEIN_PEPTIDE_ASSOCIATION_T));

  protein_peptide_association->next_association = NULL;
  return protein_peptide_association;
}

/**
 *\returns a PROTEIN_PEPTIDE_ASSOCIATION object, populated with user specified parameters
 */
PROTEIN_PEPTIDE_ASSOCIATION_T* new_protein_peptide_association(
  PEPTIDE_TYPE_T peptide_type, ///< the peptide type for the corresponding protein -in
  PROTEIN_T* parent_protein, ///< the parent of this preptide -in
  int start_idx ///< start index of the peptide in the protein sequence -in
  )
{
  PROTEIN_PEPTIDE_ASSOCIATION_T* new_association = allocate_protein_peptide_association();
  set_protein_peptide_association_peptide_type(new_association, peptide_type);
  set_protein_peptide_association_parent_protein(new_association, parent_protein);
  set_protein_peptide_association_start_idx(new_association, start_idx);
  return new_association;
}

/**
 * Frees the entire allocated protein_peptide_association linklist object
 */
void free_protein_peptide_association(
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association  ///< object to free -in 
  )
{
  if(protein_peptide_association->next_association != NULL){
    free_protein_peptide_association(protein_peptide_association->next_association);
  }    
  free(protein_peptide_association);
}

/**
 * Frees the an individual allocated protein_peptide_association object
 * assumes that new_association pointer is NULL or some other pointer exist for the rest of the linklist 
 */
void free_one_protein_peptide_association(
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association  ///< object to free -in 
  )
{
  free(protein_peptide_association);
}

//FIXME might need to change how this is printed
/**
 * Prints a peptide object to file.
 */
void print_protein_peptide_association(
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association, ///< object to print -in 
  FILE* file  ///< the out put stream -out
  )
{
  PEPTIDE_TYPE_T peptide_type = get_protein_peptide_association_peptide_type(protein_peptide_association);
  char* sequence = get_protein_sequence(protein_peptide_association->parent_protein);
  fprintf(file, "parent protein:%s\n", sequence);
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
 * Copies the entire linklist of protein_peptide_association object src to dest.
 * dest must be a heap allocated protein_peptide_association
 */
void copy_protein_peptide_association(
  PROTEIN_PEPTIDE_ASSOCIATION_T* src, ///< source protein_peptide_association -in
  PROTEIN_PEPTIDE_ASSOCIATION_T* dest ///< destination protein_peptide_association -out
  )
{
  PROTEIN_PEPTIDE_ASSOCIATION_T* next_association;
  set_protein_peptide_association_peptide_type(dest, src->peptide_type);
  set_protein_peptide_association_parent_protein(dest, src->parent_protein);
  set_protein_peptide_association_start_idx(dest, src->start_idx);
  //check if end of the linklist
  if(get_protein_peptide_association_next_association(src) != NULL){
    next_association = allocate_protein_peptide_association();
    dest->next_association = next_association;
    copy_protein_peptide_association(src->next_association, next_association);
  }
}

/**
 * sets the peptide type
 * peptide type: TRYPTIC, PARTIALLY_TRYPTIC, NON_TRYPTIC
 */
void set_protein_peptide_association_peptide_type( 
  PROTEIN_PEPTIDE_ASSOCIATION_T* new_association, ///< the protein_peptide_association to set -out   
  PEPTIDE_TYPE_T peptide_type ///< the type of the peptide -in
  )
{
  new_association->peptide_type = peptide_type;
}

/**
 * \returns the peptide type with association to the parent protein
 * peptide type: TRYPTIC, PARTIALLY_TRYPTIC, NON_TRYPTIC
 */
PEPTIDE_TYPE_T get_protein_peptide_association_peptide_type( 
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association ///< the query protein_peptide_association -in   
  )
{
  return protein_peptide_association->peptide_type;
}


/**
 * sets the parent protein
 */
void set_protein_peptide_association_parent_protein(
  PROTEIN_PEPTIDE_ASSOCIATION_T* new_association, ///< the protein_peptide_association to set -out   
  PROTEIN_T* parent_protein ///< the parent of this preptide -in  
  )
{
  new_association->parent_protein = parent_protein;

}

/**
 * \returns a pointer to the parent protein
 */
PROTEIN_T* get_protein_peptide_association_parent_protein( 
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association ///< the query protein_peptide_association -in   
  )
{
  return protein_peptide_association->parent_protein;
}

/**
 * sets the start index of the peptide in the protein sequence
 */
void set_protein_peptide_association_start_idx(
  PROTEIN_PEPTIDE_ASSOCIATION_T* new_association, ///< the protein_peptide_association to set -out   
  int start_idx ///< start index of the peptide in the protein sequence -in
  )
{
  new_association->start_idx = start_idx;
}

/**
 * \returns the start index of the peptide in the protein sequence
 */
int get_protein_peptide_association_start_idx( 
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association ///< the query protein_peptide_association -in   
  )
{
  return protein_peptide_association->start_idx;
}

/**
 * sets the next protein_peptide_association on the link list
 * assumes that the src_association's next_association feild is NULL
 */
void set_protein_peptide_association_next_association(
  PROTEIN_PEPTIDE_ASSOCIATION_T* src_association, ///< the protein_peptide_association to set -out   
  PROTEIN_PEPTIDE_ASSOCIATION_T* new_association ///< the new protein_peptide_association to add -in   
  )
{
  src_association->next_association = new_association;
}

/**
 * \returns the next protein_peptide_association on the link list
 */
PROTEIN_PEPTIDE_ASSOCIATION_T* get_protein_peptide_association_next_association( 
  PROTEIN_PEPTIDE_ASSOCIATION_T* protein_peptide_association ///< the query protein_peptide_association -in   
  )
{
  return protein_peptide_association->next_association;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
