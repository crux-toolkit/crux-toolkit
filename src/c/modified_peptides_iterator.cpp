/**
 * \file modified_peptides_iterator.cpp
 * AUTHOR: Barbara Frewen
 * DATE: April 15, 2008
 * \brief An iterator that can be used by
 * generate_peptides_iterator to include modified peptides.
 */
#include "modified_peptides_iterator.h"
#include "SpectrumZState.h"


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
  int max_aas_modified; 
  ///< return peptides with no more than this many aas modified
  BOOLEAN_T is_decoy; ///< generate target or decoy peptides
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
  if( iterator == NULL ){ return; }

  // first, try getting next from the temp list
  if( ! is_empty_linked_list( iterator->temp_peptide_list) ){
    carp(CARP_DETAILED_DEBUG,"Queue is getting next peptide from temp list");
    iterator->next_peptide = (PEPTIDE_T*)pop_front_linked_list(iterator->temp_peptide_list);

    return;
  }

  // second, try getting next from iterator
  if( ! generate_peptides_iterator_has_next( iterator->peptide_generator ) ){
    carp(CARP_DETAILED_DEBUG, "Queue has no more peptides in the generator");
    iterator->next_peptide = NULL;
    return;   // no more peptides for this iterator
  }

  // else, get the next unmodified peptide
  PEPTIDE_T* unmod_peptide = 
    generate_peptides_iterator_next( iterator->peptide_generator);
  
  IF_CARP_DETAILED_DEBUG(
    char* debugseq = get_peptide_sequence(unmod_peptide);
    carp(CARP_DETAILED_DEBUG, "Next peptide in pep_gen is %s", debugseq);
    free(debugseq);
  )

  // turn the peptide into a decoy, if required
  if( iterator->is_decoy ){
    transform_peptide_to_decoy(unmod_peptide);
  }

  // apply modifications, discard peptides that can't be modified

  // keep looking until a peptide can be modified or we run out of peptides
  carp(CARP_DETAILED_DEBUG, "Queue is looking for modifyable peptide");
  while( unmod_peptide != NULL &&
         ! is_peptide_modifiable(unmod_peptide, 
                                 iterator->peptide_mod) ){ 
    // carp(CARP_DETAILED_DEBUG, "Skipping peptide %s from the generator",
    //   get_peptide_sequence(unmod_peptide)); //memleak
    free_peptide( unmod_peptide );
    unmod_peptide = 
      generate_peptides_iterator_next( iterator->peptide_generator );
  }

  if( unmod_peptide == NULL ){ 
    // none of the remaining peptides were modifiable
    carp(CARP_DETAILED_DEBUG, "Skipped all remaining peptides in generator");
    iterator->next_peptide = NULL;
    return;
  }
  
  IF_CARP_DETAILED_DEBUG(
    char* umodseq = get_peptide_sequence(unmod_peptide);
    carp(CARP_DETAILED_DEBUG, "Iterator is modifying peptide %s",
         umodseq);
    free(umodseq);
  )
  modify_peptide(unmod_peptide, 
                 iterator->peptide_mod, 
                 iterator->temp_peptide_list,
                 iterator->max_aas_modified );
  // this put a copy in the list, get rid of the original
  free_peptide(unmod_peptide);

  if( is_empty_linked_list(iterator->temp_peptide_list) ){
    carp(CARP_DETAILED_DEBUG, "Modifier didn't return any peptides");
    iterator->next_peptide = NULL;
    return;
  }

  // now set next_peptide to the first in the list and move list forward
  iterator->next_peptide = (PEPTIDE_T*)pop_front_linked_list(iterator->temp_peptide_list);
  if( iterator->next_peptide == NULL ){
    printf("Iterator's next peptide was lost\n");
  }
  
  IF_CARP_DETAILED_DEBUG(
    char* seq = 
    get_peptide_modified_sequence_with_masses(iterator->next_peptide, FALSE);
    carp(CARP_DETAILED_DEBUG, "Queue set next peptide as %s", seq);
    free(seq);
  )
    /*
  // now check that it does not exceed the max number of modified aas
  int max_aas_moded = get_int_parameter("max-aas-modified");
  if( count_peptide_modified_aas(iterator->next_peptide) > max_aas_moded ){
    free_peptide(iterator->next_peptide);
    queue_next_peptide(iterator);
  }
    */

}

/* Public functions */
/**
 * \brief Create a new modified_peptides_iterator for all peptides in
 * database or index.
 *
 * The returned iterator is initialized with the first peptide queued
 * up and ready to return.  Also creates a generate_peptides_iterator
 * from which it gets the peptides to modify. Peptides are from mass
 * 'min mass' + pmod->delta_mass to 'max mass' + pmod->delta_mass with
 * min and max mass taken from parameter.c.  All other peptide
 * specifications are taken from parameter.c.  If no peptides meet the
 * specifications, an iterator is still returned and when given to
 * has_next will always return FALSE.
 * 
 * \returns A newly allocated modified_peptides_iterator.
 */
MODIFIED_PEPTIDES_ITERATOR_T* new_modified_peptides_iterator(
  PEPTIDE_MOD_T* pmod, ///< Peptide mod to apply
  Index* index,      ///< Index from which to draw peptides OR
  Database* dbase    ///< Database from which to draw peptides
){
  if( index == NULL && dbase == NULL ){
    carp(CARP_FATAL, 
         "Cannot create modified peptides iterator from NULL protein source");
  }

  carp(CARP_DETAILED_DEBUG, 
       "Creating modified peptides iterator with %d aamods",
       peptide_mod_get_num_aa_mods(pmod));
  MODIFIED_PEPTIDES_ITERATOR_T* new_iterator = 
    allocate_modified_peptides_iterator();

  // init the maximum number of aas that can be modified
  new_iterator->max_aas_modified = get_int_parameter("max-aas-modified");

  // init the peptide list and peptide_mod
  new_iterator->temp_peptide_list = new_empty_list();
  new_iterator->peptide_mod = pmod;

  // get min and max masses
  double min_mass = get_double_parameter("min-mass");
  double max_mass = get_double_parameter("max-mass");

  carp(CARP_DETAILED_DEBUG, "min mass is %.2f, max is %.2f", 
       min_mass, max_mass);

  // TODO: this could be removed since we are searching all masses
  //       would remove warning
  // get the mass difference for the mod
  double delta_mass = peptide_mod_get_mass_change(pmod);
  min_mass += delta_mass;
  max_mass += delta_mass;

  new_iterator->peptide_generator = 
    new_generate_peptides_iterator_from_mass_range(min_mass,
                                                   max_mass,
                                                   index, dbase);

  // queue first peptide
  carp(CARP_DETAILED_DEBUG, "Queueing first peptide.");
  queue_next_peptide( new_iterator );

  return new_iterator;
}

/**
 * \brief Create a new modified_peptides_iterator for a specified mass.
 *
 * The returned iterator is initialized with the first peptide queued
 * up and ready to return.  Also creates a generate_peptides_iterator
 * from which it gets the peptides to modify. Peptides are of mass +/-
 * mass-window taken from parameter.c.  All other peptide
 * specifications are taken from parameter.c.  If no peptides meet the
 * specifications, an iterator is still returned and when given to
 * has_next will always return FALSE.  Returns either target or decoy
 * peptides depending on the is_decoy argument.
 * 
 * \returns A newly allocated modified_peptides_iterator.
 */
MODIFIED_PEPTIDES_ITERATOR_T* new_modified_peptides_iterator_from_mass(
  double mass,         ///< Target mass of peptides BEFORE modification
  PEPTIDE_MOD_T* pmod, ///< Peptide mod to apply
  BOOLEAN_T is_decoy,  ///< generate decoy peptides
  Index* index,      ///< Index from which to draw peptides OR
  Database* dbase    ///< Database from which to draw peptides
  ){
  MODIFIED_PEPTIDES_ITERATOR_T* new_iterator = 
    allocate_modified_peptides_iterator();

  // init max aas modified
  new_iterator->max_aas_modified = get_int_parameter("max-aas-modified");

  // init the peptide list
  new_iterator->temp_peptide_list = new_empty_list();

  // set peptide_mod field
  new_iterator->peptide_mod = pmod;

  // set the decoy field
  new_iterator->is_decoy = is_decoy;

  // get the mass difference
  double delta_mass = peptide_mod_get_mass_change(pmod);

  // create peptide_generator
  new_iterator->peptide_generator = 
    new_generate_peptides_iterator_from_mass(mass - delta_mass, index, dbase);

  // queue first peptide
  queue_next_peptide( new_iterator );

  return new_iterator;
}


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
  Index* index,      ///< Index from which to draw peptides OR
  Database* dbase    ///< Database from which to draw peptides
  ) {

  WINDOW_TYPE_T precursor_window_type = 
    get_window_type_parameter("precursor-window-type");
  double window = get_double_parameter("precursor-window");
  double min_mass = 0;
  double max_mass = 0;
  
  MODIFIED_PEPTIDES_ITERATOR_T* new_iterator = 
    allocate_modified_peptides_iterator();

  // init max aas modified
  new_iterator->max_aas_modified = get_int_parameter("max-aas-modified");
  
  // init the peptide list
  new_iterator->temp_peptide_list = new_empty_list();
  
  // set peptide_mod field
  new_iterator->peptide_mod = pmod;

  // set is_decoy field
  new_iterator->is_decoy = is_decoy;
  
  // get the mass difference
  double delta_mass = peptide_mod_get_mass_change(pmod);

  if (precursor_window_type == WINDOW_MASS) {
    //TODO: should we change MASS_H to MASS_PROTON and regenerate the smokes?
    double mass = (mz - MASS_H) * charge - delta_mass; 
    min_mass = mass - window;
    max_mass = mass + window;
  } else if (precursor_window_type == WINDOW_MZ) {
      double min_mz = mz - window;
      double max_mz = mz + window;
      min_mass = (min_mz - MASS_PROTON) * (double)charge - delta_mass;
      max_mass = (max_mz - MASS_PROTON) * (double)charge - delta_mass;
  } else if (precursor_window_type == WINDOW_PPM) {
    double mass = (mz - MASS_PROTON) * (double)charge - delta_mass;
    min_mass = mass / (1.0 + window * 1e-6);
    max_mass = mass / (1.0 - window * 1e-6);
  } else {
    carp(CARP_FATAL,"Invalid window type");
  }
  
  carp(CARP_DETAILED_DEBUG, "Generating peptides from %.5f to %.5f for "
       "precursor m/z %.2f, mass %.5f, window size %.2f", min_mass, max_mass, 
       mz, ((mz - MASS_PROTON) * (double)charge - delta_mass), window);
  // create peptide_generator
    new_iterator->peptide_generator = 
      new_generate_peptides_iterator_from_mass_range(min_mass, max_mass, index, dbase);
  
  // queue first peptide
  queue_next_peptide( new_iterator );

  return new_iterator;
}

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
MODIFIED_PEPTIDES_ITERATOR_T* new_modified_peptides_iterator_from_zstate(
  double mz, ///< Spectrum precursor mz.
  SpectrumZState& zstate, ///< zstate to search.
  PEPTIDE_MOD_T* pmod, ///< Peptide mod to apply
  BOOLEAN_T is_decoy,  ///< generate decoy peptides
  Index* index,      ///< Index from which to draw peptides OR
  Database* dbase    ///< Database from which to draw peptides
  ) {

  WINDOW_TYPE_T precursor_window_type = 
    get_window_type_parameter("precursor-window-type");
  double window = get_double_parameter("precursor-window");
  double min_mass = 0;
  double max_mass = 0;
  
  MODIFIED_PEPTIDES_ITERATOR_T* new_iterator = 
    allocate_modified_peptides_iterator();

  // init max aas modified
  new_iterator->max_aas_modified = get_int_parameter("max-aas-modified");
  
  // init the peptide list
  new_iterator->temp_peptide_list = new_empty_list();
  
  // set peptide_mod field
  new_iterator->peptide_mod = pmod;

  // set is_decoy field
  new_iterator->is_decoy = is_decoy;
  
  // get the mass difference
  double delta_mass = peptide_mod_get_mass_change(pmod);

  if (precursor_window_type == WINDOW_MASS) {
    double mass = zstate.getNeutralMass() - delta_mass;
    min_mass = mass - window;
    max_mass = mass + window;
  } else if (precursor_window_type == WINDOW_MZ) {
    double min_mz = mz - window;
    double max_mz = mz + window;
    min_mass = (min_mz - MASS_PROTON) * (double)zstate.getCharge() - delta_mass;
    max_mass = (max_mz - MASS_PROTON) * (double)zstate.getCharge() - delta_mass;
  } else if (precursor_window_type == WINDOW_PPM) {
    double mass = zstate.getNeutralMass() - delta_mass;
    min_mass = mass / (1.0 + window * 1e-6);
    max_mass = mass / (1.0 - window * 1e-6);
    carp(CARP_DEBUG,"mass:%f charge:%i min_mass:%f max_mass:%f",mass, zstate.getCharge(), min_mass, max_mass);
  } else {
    carp(CARP_FATAL,"Invalid window type");
  }
  
  // create peptide_generator
    new_iterator->peptide_generator = 
      new_generate_peptides_iterator_from_mass_range(min_mass, max_mass, index, dbase);
  
  // queue first peptide
  queue_next_peptide( new_iterator );

  return new_iterator;
}


/**
 * \brief Check to see if the iterator has more peptides to return,
 * i.e. a call to next() will return non-NULL.
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
 * Filters out peptides that have too many aas modified.
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
    // maybe don't??
    //free_peptide_mod( iterator->peptide_mod );
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
