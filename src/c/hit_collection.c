/*********************************************************************//**
 * FILE: \file hit_collection.c
 * AUTHOR: Aaron Klammer
 * DESCRIPTION: \brief A collection of hits.
 * CREATE DATE: 2008 March 11
 * REVISION: $Revision: 1.2 $
 ****************************************************************************/
#include "hit_collection.h"

/**
 * \struct hit_collection
 * \brief An object that contains a set of hit objects.
 */
struct hit_collection {
  HIT_T* hits[_MAX_NUMBER_HITS]; ///< Array of hit objects
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
 HIT_COLLECTION_T* match_collection ///< the match collection -in
 )
{
  HIT_COLLECTION_T* hit_collection = allocate_hit_collection();
  carp(CARP_DETAILED_DEBUG, "Finished creating hit collection");
  return hit_collection;
}

/**
 * \brief Print the hits to file in sqt format.
 * \returns TRUE, if successfully prints to the output
 */
BOOLEAN_T print_hit_collection(
  FILE* output,                  ///< the output file -out
  HIT_COLLECTION_T* hit_collection
  ){

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
    die("Null hit collection passed to hit iterator");
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
 *\returns the next the hit struct
 */
HIT_T* hit_iterator_next(
  HIT_ITERATOR_T* hit_iterator ///< the working hit iterator -in
  )
{
  return hit_iterator->hit_collection->hit[hit_iterator->hit_idx++];
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

