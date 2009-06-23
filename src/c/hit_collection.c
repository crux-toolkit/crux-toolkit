/*********************************************************************//**
 * FILE: \file hit_collection.c
 * AUTHOR: Aaron Klammer
 * DESCRIPTION: \brief A collection of hits.
 * CREATE DATE: 2008 March 11
 * REVISION: $Revision: 1.9 $
 ****************************************************************************/
#include "hit_collection.h"

/**
 * \struct hit_collection
 * \brief An object that contains a set of hit objects.
 */
 struct hit_collection {
  HIT_T* hits[_MAX_NUMBER_HITS]; ///< Array of hit objects
  HASH_T* id_to_hit;
  int hit_total; ///< Total number of hits
 };

/**
 *\struct hit_iterator
 *\brief An object that iterates over the hit objects in the
 * specified hit_collection for the specified score type (SP, XCORR)
 */
struct hit_iterator {
  HIT_COLLECTION_T* hit_collection; ///< the hit collection to iterate over
  int hit_idx;   ///< current hit to return
  int hit_total; ///< total number of hits
};

/**
 * \returns An (empty) hit_collection object.
 */
HIT_COLLECTION_T* allocate_hit_collection(void)
{
  HIT_COLLECTION_T* hit_collection =
    (HIT_COLLECTION_T*)mycalloc(1, sizeof(HIT_COLLECTION_T));
    
  // carp(CARP_DETAILED_DEBUG, "Allocated hit collection");
  
  hit_collection->id_to_hit = new_hash(1000000);  // TODO
  return hit_collection;
}

/**
 * /brief Free the memory allocated for a hit collection
 */
void free_hit_collection(
  HIT_COLLECTION_T* hit_collection ///< the hit collection to free -out
  )
{
  // carp(CARP_DETAILED_DEBUG, "Freeing hit collection"); 
  // decrement the pointer count in each hit object
  // MEMLEAK
  while(hit_collection->hit_total > 0){
    --hit_collection->hit_total;
    free_hit(hit_collection->hits[hit_collection->hit_total]);
    hit_collection->hits[hit_collection->hit_total] = NULL;
  }
  
  free(hit_collection);
}

/**
 * \brief Creates a new hit collection from a scored match collection.
 * \detail This is the main protein assembly routine.  Allocates memory for
 * the hit collection. 
 * \returns A new hit_collection object.
 */
HIT_COLLECTION_T* new_hit_collection_from_match_collection(
 MATCH_COLLECTION_T* match_collection ///< the match collection -in
 // PROTEIN_SCORER_TYPE_T protein_scorer_type  ///< the type of protein score -in
 )
{
  HIT_COLLECTION_T* hit_collection = allocate_hit_collection();
  
  // assert pvalue only

  // three data structures, one for storing peptide2max_score, 
  // peptide2number proteins, and protein2score
  // an array of all peptides in results
  HASH_T* peptide_to_max_match = new_hash(1000000); // TODO 

  SCORER_TYPE_T scorer_type = LOGP_BONF_WEIBULL_XCORR;

  MATCH_ITERATOR_T* match_iterator = 
    new_match_iterator(match_collection, scorer_type, FALSE);


  // iterate through the matches
  MATCH_T* match = NULL;
  while(match_iterator_has_next(match_iterator)){

    match = match_iterator_next(match_iterator);
    FLOAT_T new_score = get_match_score(match, scorer_type);

    char* peptide_sequence = get_match_sequence(match);
    FLOAT_T max_score = 0.0;

    MATCH_T* max_match = 
      (MATCH_T*) get_hash_value(peptide_to_max_match, peptide_sequence);
    if (max_match != NULL){
      max_score = get_match_score(max_match, scorer_type);
    }

    // note the peptide max score in peptide2max_score hash
    if (max_score < new_score){
      add_or_update_hash(
          peptide_to_max_match, peptide_sequence, (void*)match);
    }
  }
  free_match_iterator(match_iterator);


  // examine each peptide
  HASH_ITERATOR_T* hash_iterator = new_hash_iterator(peptide_to_max_match);
  while(hash_iterator_has_next(hash_iterator)){
    char* key = hash_iterator_next(hash_iterator);
    MATCH_T* match = (MATCH_T*) get_hash_value(peptide_to_max_match, key);
    PEPTIDE_T* peptide = get_match_peptide(match);

    // count the number of proteins for this peptide
    PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator 
      = new_peptide_src_iterator(peptide);
    int protein_count = 0;
    while(peptide_src_iterator_has_next(peptide_src_iterator)){
      peptide_src_iterator_next(peptide_src_iterator);
      protein_count++;
    }
    free_peptide_src_iterator(peptide_src_iterator);

    // take nth root of this peptide's max score. this has the effect of 
    // distribution the peptide score evenly across all its proteins
    FLOAT_T nth_root_score = get_match_score(match, scorer_type) / protein_count;

    // note the score in each of the protein hits in the hit collection
    peptide_src_iterator = new_peptide_src_iterator(peptide);
    while(peptide_src_iterator_has_next(peptide_src_iterator)){
      PEPTIDE_SRC_T* peptide_src 
        = peptide_src_iterator_next(peptide_src_iterator);
      PROTEIN_T* protein = get_peptide_src_parent_protein(peptide_src);
      hit_collection_add_protein_score(hit_collection, protein, nth_root_score);
    }
    free_peptide_src_iterator(peptide_src_iterator);
  }
  free_hash_iterator(hash_iterator);
  
  return hit_collection;
}

/**
 * Adds the current score to the appropriate hit for the protein id
 * If protein id does not exist, add it to the collection
 */
BOOLEAN_T hit_collection_add_protein_score(
    HIT_COLLECTION_T* hit_collection,
    PROTEIN_T* protein,
    FLOAT_T score
  ){
  char* protein_id = get_protein_id(protein);
  HIT_T* hit = NULL;
  if (hit_collection_has_protein_id(hit_collection, protein_id)){
    hit = (HIT_T*) get_hash_value(hit_collection->id_to_hit, protein_id);
  } else {
    hit = new_hit();
    set_hit_protein(hit, protein);
    hit_collection_add_hit(hit_collection, hit);
  }
  hit_increment_score(hit, score);
  return TRUE;
}

/**
 * \brief Print the hits to file in sqt format.
 * \returns TRUE, if successfully prints to the output
 */
BOOLEAN_T print_hit_collection(
  FILE* output,                  ///< the output file -out
  HIT_COLLECTION_T* hit_collection
  ){

  HIT_ITERATOR_T* hit_iterator = new_hit_iterator(hit_collection);
  while (hit_iterator_has_next(hit_iterator)){
    HIT_T* hit = hit_iterator_next(hit_iterator);
    print_hit(output, hit);
  }  

  return TRUE;
}

/**
 * Does this hit collection have the hit associated with this protein id?
 */
BOOLEAN_T hit_collection_has_protein_id(
  HIT_COLLECTION_T* hit_collection,
  char* protein_id){
  if (get_hash_value(hit_collection->id_to_hit, protein_id) == NULL){
    return FALSE;
  } else {
    return TRUE;
  }
}


/**
 * \brief Add the hit to the hit collection.
 */
BOOLEAN_T hit_collection_add_hit(
    HIT_COLLECTION_T* hit_collection,
    HIT_T* hit
    ){
  add_or_update_hash(
      hit_collection->id_to_hit, 
      get_protein_id(get_hit_protein(hit)), 
      (void *)get_hit_protein(hit));
  hit_collection->hits[hit_collection->hit_total++] = hit;
  return TRUE;
}

/**
 * hit_iterator routines!
 */

/**
 *\returns a new memory allocated hit iterator
 */
HIT_ITERATOR_T* new_hit_iterator(
  HIT_COLLECTION_T* hit_collection ///< the hit collection to iterate -out
  ){
  if (hit_collection == NULL){
    carp(CARP_FATAL, "Null hit collection passed to hit iterator");
  }
  
  // allocate a new hit iterator
  HIT_ITERATOR_T* hit_iterator = 
    (HIT_ITERATOR_T*)mycalloc(1, sizeof(HIT_ITERATOR_T));
  
  // set items
  hit_iterator->hit_collection = hit_collection;
  hit_iterator->hit_idx = 0;
  hit_iterator->hit_total = hit_collection->hit_total;

  return hit_iterator;
}

/**
 * Does the hit_iterator have another hit object to return?
 * \returns TRUE, if hit iterator has a next hit, else FALSE
 */
BOOLEAN_T hit_iterator_has_next(
  HIT_ITERATOR_T* hit_iterator ///< the working  hit iterator -in
  )
{
  return (hit_iterator->hit_idx < hit_iterator->hit_total);
}

/**
 * \returns the next the hit struct
 */
HIT_T* hit_iterator_next(
  HIT_ITERATOR_T* hit_iterator ///< the working hit iterator -in
  )
{
  return hit_iterator->hit_collection->hits[hit_iterator->hit_idx++];
}

/**
 * free the memory allocated iterator
 */
void free_hit_iterator(
  HIT_ITERATOR_T* hit_iterator ///< the hit iterator to free
  )
{
  if (hit_iterator != NULL){
    free(hit_iterator);
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

