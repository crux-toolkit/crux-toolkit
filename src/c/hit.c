/*************************************************************************//**
 * \file hit.c
 * AUTHOR: Aaron Klammer
 * CREATE DATE: 2008 March 11
 * DESCRIPTION: Object for collecting the evidence for a particular protein.
 * REVISION: $Revision: 1.4 $
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
  int num_scores;
  // Possible fields 
  // BOOLEAN_T null_protein; ///< Is this hit a null protein hit?
};

/**
 * \returns a new memory allocated hit
 */
HIT_T* new_hit(void){
  HIT_T* hit = (HIT_T*)mycalloc(1, sizeof(HIT_T));
  
  ++hit->pointer_count;
  hit->num_scores = 0;

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
  carp(CARP_DETAILED_INFO, "Incrementing score %d", hit->num_scores);
  hit->num_scores += 1;
  carp(CARP_DETAILED_INFO, "Incrementing score %d", hit->num_scores);
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
void hit_recalibrate_score(
  HIT_T* hit
  )
{
  double ln_k = - hit->score; 
  int n = hit->num_scores; 

  double new_score = 0.0;

  int i;
  for (i=0; i<n-1; i++){
    double inner = 1.0;
    int j;
    for (j=1; j<i+1; j++){
      inner *= -ln_k;
      inner /= j;
    }
    carp(CARP_DETAILED_INFO, "Inner for %i = %.9f", i, inner);
    new_score += inner;
  }

  carp(CARP_DETAILED_INFO, "New_score = %.9f", new_score);
  new_score = -log(new_score);
  new_score += ln_k;
  hit->score = new_score;
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

