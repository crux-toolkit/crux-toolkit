/**
 * \file modified_peptides_iterator.c
 * AUTHOR: Barbara Frewen
 * DATE: April 15, 2008
 * DESCRIPTION: An iterator that can be used by
 * generate_peptides_iterator to include modified peptides.
 * $Revision: 1.1.2.2 $
 */
#include "modified_peptides_iterator.h"

/* Private data structures */
/**
 * \struct modified_peptides_iterator_t
 * \brief An object for keeping track of which peptide to return next.
 */
struct modified_peptides_iterator_t{
  GENERATE_PEPTIDES_ITERATOR_T* peptide_generator; 
  ///< for getting unmodified peptides
  PEPTIDE_MOD_T* peptide_mod;///< the modification to apply to peptides
  PEPTIDE_T* next_peptide;///< peptide queued to return next
  LINKED_LIST_T* temp_peptide_list;///< storage for modified peptides
};

/* Private functions */
MODIFIED_PEPTIDES_ITERATOR_T* allocate_modified_peptides_iterator(){
  MODIFIED_PEPTIDES_ITERATOR_T* new_iterator = 
    (MODIFIED_PEPTIDES_ITERATOR_T*)mycalloc(1, 
                                      sizeof(MODIFIED_PEPTIDES_ITERATOR_T));

  return new_iterator;
}

/**
 * \brief Queue next_peptide to be returned via
 * modified_peptides_iterator_next. 
 *
 * Private function that can be called on a newly created iterator or
 * in mid-iteration. Takes unmodified peptides from a
 * generate_peptides_iterator and applies the modification to each.
 * 
 * Filters out peptides that can't be modified. Stores multiple
 * versions of same peptide that can have the mod applied in more than
 * one way (ie on different residues) in temp_peptide_list.  Deletes
 * nodes from list as used. Sets next_peptide to NULL when there are
 * no more to be returned.
 * 
 */
void queue_next_peptide(
 MODIFIED_PEPTIDES_ITERATOR_T* iterator ///< working iterator
 ){
  if( iterator == NULL ){
    return;
  }

  // first, try getting next from the temp list
  //if( iterator->temp_peptide_list != NULL ){
  if( ! is_empty_linked_list( iterator->temp_peptide_list) ){
    //printf("Getting next peptide from temp list\n");
    /*
    iterator->next_peptide = 
      (PEPTIDE_T*)get_data_linked_list( iterator->temp_peptide_list );
    iterator->temp_peptide_list = 
      get_next_free_this_linked_list( iterator->temp_peptide_list );
    */
    iterator->next_peptide =pop_front_linked_list(iterator->temp_peptide_list);
    return;
  }

  // second, try getting next from iterator
  if( ! generate_peptides_iterator_has_next( iterator->peptide_generator ) ){
    //printf("No more peptides in the generator\n");
    iterator->next_peptide = NULL;
    return;   // no more peptides for this iterator
  }

  // else, get the unmodified peptide
  iterator->next_peptide = 
    generate_peptides_iterator_next( iterator->peptide_generator);

  printf("Next peptide from pep_gen is %s\n", get_peptide_sequence(iterator->next_peptide));

  // apply modifications, discarding peptides that can't be modified

  // keep looking until a peptide can be modified or we run out of peptides
  while( iterator->next_peptide != NULL &&
         ! is_peptide_modifiable(iterator->next_peptide, 
                                 iterator->peptide_mod) ){ 
    printf("Skipping peptide from the generator\n");
    free_peptide( iterator->next_peptide );
    iterator->next_peptide = 
      generate_peptides_iterator_next( iterator->peptide_generator);
  }

  if( iterator->next_peptide == NULL ){ 
    // none of the remaining peptides were modifiable
    //printf("Skipped all remaining peptides in generator\n");
    return;
  }

  printf("Modifying peptide\n");
  modify_peptide( iterator->next_peptide, 
                  iterator->peptide_mod, 
                  &(iterator->temp_peptide_list) );

  if( iterator->temp_peptide_list == NULL ){
    printf("Modifier didn't return any peptides\n");
  }

  // now set next_peptide to the first in the list and move list forward
  free_peptide(iterator->next_peptide); //we have a copy in the temp_list
  /*
  iterator->next_peptide = 
    get_data_linked_list( iterator->temp_peptide_list );
  iterator->temp_peptide_list = 
    get_next_free_this_linked_list( iterator->temp_peptide_list );
  */
  iterator->next_peptide = pop_front_linked_list(iterator->temp_peptide_list);
  if( iterator->next_peptide == NULL ){
    printf("Iterator's next peptide was lost\n");
  }
  //printf("Next peptide is %s\n", get_peptide_sequence(iterator->next_peptide));


}

/* Public functions */

/**
 * \brief Create a new modified_peptides_iterator.
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
MODIFIED_PEPTIDES_ITERATOR_T* new_modified_peptides_iterator(
  double mass,         ///< Target mass of peptides
  PEPTIDE_MOD_T* pmod, ///< Peptide mod to apply
  INDEX_T* index,      ///< Index from which to draw peptides OR
  DATABASE_T* dbase    ///< Database from which to draw peptides
  ){
  MODIFIED_PEPTIDES_ITERATOR_T* new_iterator = 
    allocate_modified_peptides_iterator();

  // init the peptide list
  new_iterator->temp_peptide_list = new_empty_list();

  // set peptide_mod field
  new_iterator->peptide_mod = pmod;

  // get the mass difference
  double delta_mass = peptide_mod_get_mass_change(pmod);

  // create peptide_generator
  new_iterator->peptide_generator = 
    new_generate_peptides_iterator_from_mass(mass + delta_mass, index, dbase);

  // queue first peptide
  printf("Queuing first peptide\n");
  queue_next_peptide( new_iterator );

  return new_iterator;
}

/**
 * \brief Check to see if the iterator has more peptides to return,
 * i.e. a call to has_next will return non-NULL.
 * \returns TRUE if the iterator has more peptides to return.
 */
BOOLEAN_T modified_peptides_iterator_has_next(
  MODIFIED_PEPTIDES_ITERATOR_T* iterator){
  if( iterator == NULL ){
    return FALSE;
  }

  if( iterator->next_peptide == NULL ){
    return FALSE;
  }// else
  return TRUE;
}

/**
 * \brief Return the next peptide or NULL if no peptides remain.
 * \returns A modified peptide.
 */
PEPTIDE_T* modified_peptides_iterator_next(
  MODIFIED_PEPTIDES_ITERATOR_T* iterator){
  if( iterator == NULL ){
    return NULL;
  }

  PEPTIDE_T* returnme = iterator->next_peptide;
  queue_next_peptide(iterator);

  return returnme;
}

/**
 * \brief Free the memory used by this iterator.
 */
void free_modified_peptides_iterator(
  MODIFIED_PEPTIDES_ITERATOR_T* iterator){
  if( iterator ){
    delete_linked_list( iterator->temp_peptide_list );
    free_peptide( iterator->next_peptide );
    free_peptide_mod( iterator->peptide_mod );
    free_generate_peptides_iterator( iterator->peptide_generator);
    free( iterator );
  }
}

/**
 * \brief Void wrapper of modified_peptides_iterator_has_next to be
 * used by generate_peptides_iterator. 
 */
BOOLEAN_T void_modified_peptides_iterator_has_next(
  void* modified_peptides_iterator){
  return modified_peptides_iterator_has_next(
           (MODIFIED_PEPTIDES_ITERATOR_T*)modified_peptides_iterator);
}

/**
 * \brief Void wrapper of modified_peptides_iterator_next to be
 * used by generate_peptides_iterator. 
 */
PEPTIDE_T* void_modified_peptides_iterator_next(
  void* modified_peptides_iterator){
  return modified_peptides_iterator_next(
           (MODIFIED_PEPTIDES_ITERATOR_T*)modified_peptides_iterator);
}

/**
 * \brief Void wrapper of modified_peptides_iterator_free to be
 * used by generate_peptides_iterator. 
 */
void void_modified_peptides_iterator_free(
  void* modified_peptides_iterator){
  return free_modified_peptides_iterator(
           (MODIFIED_PEPTIDES_ITERATOR_T*)modified_peptides_iterator);
}
