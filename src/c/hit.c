/*************************************************************************//**
 * \file hit.c
 * AUTHOR: Aaron Klammer
 * CREATE DATE: 2008 March 11
 * DESCRIPTION: Object for collecting the evidence for a particular protein.
 * REVISION: $Revision: 1.3 $
 ****************************************************************************/
#include "hit.h"

/**
 *\struct hit 
 *\brief Object for collecting the evidence for a particular protein.
 */
struct hit{
  PROTEIN_T* protein; ///< the protein for which we collecting evidence
  double score; ///< The "score" for this protein hit. Should probably be
                ///< replaced with an array at some point
  int pointer_count; ///< the number of pointers to this hit object 
  // Possible fields 
  // BOOLEAN_T null_protein; ///< Is this hit a null protein hit?
};

/**
 * \returns a new memory allocated hit
 */
HIT_T* new_hit(void){
  HIT_T* hit = (HIT_T*)mycalloc(1, sizeof(HIT_T));
  
  ++hit->pointer_count;

  return hit;
}

/**
 * free the memory allocated hit
 * spectrum is not freed by hit
 */
void free_hit(
  HIT_T* hit ///< the hit to free -in
  )
{
  --hit->pointer_count;
  
  // only free hit when pointer count reaches
  if(hit->pointer_count == 0){

    /* if (hit->peptide != NULL){
      free_peptide(hit->peptide);
    } */

    free(hit);  
  }
}

/**
 * \returns Increments the hit score by score
 */
void hit_increment_score(
    HIT_T* hit,
    double score){
  hit->score = hit->score + score;
}

/**
 * \returns Gets the hit protein;
 */
PROTEIN_T* get_hit_protein(
    HIT_T* hit){
  return hit->protein;
}

/**
 * \returns Sets the hit protein
 */
void set_hit_protein(
    HIT_T* hit,
    PROTEIN_T* protein){
  hit->protein = protein;
}

/**
 * print the information of the hit
 */
void print_hit(
  FILE* file,  ///< output stream -out
  HIT_T* hit ///< the hit to print -in  
  )
{
  fprintf(file, "%s\t%.6f\n", get_protein_id(hit->protein), hit->score);
}

/**
 * \returns the hit score in the hit object
 */
double get_hit_score(
  HIT_T* hit          ///< the hit to work -in  
  )
{
  return hit->score;
}

/**
 *Increments the pointer count to the hit object
 */
void increment_hit_pointer_count(
  HIT_T* hit ///< the hit to work -in  
  )
{
  ++hit->pointer_count;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

