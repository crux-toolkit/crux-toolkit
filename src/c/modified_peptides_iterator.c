// write header etc

/**
 * \brief Queue next_peptide to be returned via
 * generate_peptides_iterator_get_next. 
 *
 * Private function that can be called on a newly created iterator or
 * in mid-iteration. 
 * Filters peptides for modifciations and applies modifications,
 * storing multiple versions of same peptide with mods on different
 * residues in temp_peptide_list.  Deletes nodes from list as used.
 * 
 */
void queue_next_peptide(
 GENERATE_PEPTIDES_ITERATOR_T* gp_iterator ///< working iterator
 ){

  // first, try getting next from the temp list
  if( gp_iterator->temp_peptide_list != NULL ){
    gp_iterator->next_peptide = 
      (PEPTIDE_T*)get_data_linked_list( gp_iterator->temp_peptide_list );
    gp_iterator->temp_peptide_list = 
      get_next_free_this_linked_list( gp_iterator->temp_peptide_list );
  }

  // second, try getting next from iterator
  if( ! gp_iterator->has_next(gp_iterator->iterator) ){
    gp_iterator->next_peptide = NULL;
    return;
  }

  gp_iterator->next_peptide = gp_iterator->next( gp_iterator->iterator );

  // keep looking until a peptide can be modified
  while( ! is_peptide_modifiable(gp_iterator->next_peptide, 
                                 gp_iterator->peptide_mod) ){ 
    gp_iterator->next_peptide = gp_iterator->next( gp_iterator->iterator );
  }

  if( gp_iterator->next_peptide == NULL ){ 
    // none of the remaining peptides were modifiable
    return;
  }

  modify_peptide( gp_iterator->next_peptide, 
                  gp_iterator->peptide_mod, 
                  &(gp_iterator->temp_peptide_list) );

  // now set next_peptide to the first in the list and move list forward
  gp_iterator->next_peptide = 
    get_data_linked_list( gp_iterator->temp_peptide_list );
  gp_iterator->temp_peptide_list = 
    get_next_free_this_linked_list( gp_iterator->temp_peptide_list );

}
