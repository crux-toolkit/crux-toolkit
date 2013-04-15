/**
 * \file modifications.cpp
 * \brief Datatypes and methods for peptide modifications
 *
 * Two data structures define modifications.  The AA_MOD_T is the most
 * basic type.  It is the information provided by the user: mass
 * change caused by this mod, amino acids which may be modified in
 * this way, and the maximum number of this type of modification which
 * may occur on one peptide.  AA_MODs are defined in modifications.h.
 * A collection of AA_MODs that may occur
 * on some peptide are represented as a PEPTIDE_MOD_T.  This stores
 * a list of AA_MODS and the net mass change experienced by the
 * peptide.  PEPTIDE_MODs are defined here.
 * AA_MODs are instantiated once after parsing the parameter file.  All
 * possible PEPTIDE_MODs are calcualted once and reused for each
 * spectrum search.  One PEPTIDE_MOD corresponds to one mass window
 * that must be searched.
 * 
 * $Revision: 1.6 $
 */

#include "peptide_modifications.h"
#include "Peptide.h"

using namespace Crux;

/* Private data types, typedefed in objects.h */
/**
 * \struct _peptide_mod
 *  A collection of aa modifications that can occur on a single
 *  peptide.  Includes the list of aa_mods and the net mass
 *  change of all mods.  Each aa_mod is included in the list as many
 *  times at it appears in the collection. There is no information about
 *  which residues in a peptide are modified.
 */
struct _peptide_mod{
  double mass_change;     ///< the net mass change for the peptide
  int num_mods;           ///< the number of items in the list_of_mods
  int aa_mod_counts[MAX_AA_MODS]; ///< the number of each kind of aa mod
                                  ///as listed in parameter.c
};

/* Private functions */
int apply_mod_to_list(
  LINKED_LIST_T* mod_seqs, ///< a pointer to a list of seqs
  AA_MOD_T* mod_to_apply,  ///< the specific mod to apply
  int num_copies,
  int max_aas_modified
);

int apply_mod_to_seq(
  MODIFIED_AA_T* seq, ///< the seq to modify
  AA_MOD_T* mod,      ///< the mod to apply
  int skip_n,         ///< skip over n aas modified by this mod
  LINKED_LIST_T* return_list); ///< the newly modified versions of


/* Definitions of public methods */

/**
 * \brief Allocate a PEPTIDE_MOD and set all fields to default values
 * (i.e. no modifications).
 * \returns A heap allocated PEPTIDE_MOD with default values.
 */
PEPTIDE_MOD_T* new_peptide_mod(){
  PEPTIDE_MOD_T* new_mod = (PEPTIDE_MOD_T*)mymalloc(sizeof(PEPTIDE_MOD_T));
  new_mod->mass_change = 0;
  new_mod->num_mods = 0;
  int i=0;
  for(i=0; i<MAX_AA_MODS;i++){
    new_mod->aa_mod_counts[i] = 0;
  }
  // add any fixed terminal mods by setting their count to 1
  int index = get_fixed_mod_index(N_TERM); 
  AA_MOD_T** list = NULL;
  get_all_aa_mod_list(&list);
  if( index > -1 ){
    new_mod->aa_mod_counts[index] = 1;
    new_mod->num_mods++;
    new_mod->mass_change += aa_mod_get_mass_change(list[index]);
  }
  index = get_fixed_mod_index(C_TERM); 
  if( index > -1 ){
    new_mod->aa_mod_counts[index] = 1;
    new_mod->num_mods++;
    new_mod->mass_change += aa_mod_get_mass_change(list[index]);
  }

  return new_mod;
}

/**
 * \brief Allocate a new peptide mod and copy contents of given mod
 * into it.
 * \returns A pointer to a new peptide mod which is a copy of the
 * given one.
 */
PEPTIDE_MOD_T* copy_peptide_mod(PEPTIDE_MOD_T* original){
  PEPTIDE_MOD_T* copy = new_peptide_mod();
  copy->mass_change = original->mass_change;
  copy->num_mods = original->num_mods;
  int i=0;
  for(i=0; i<MAX_AA_MODS;i++){
    copy->aa_mod_counts[i] = original->aa_mod_counts[i];
  }

  return copy;
}

/**
 * \brief Free the memory for a PEPTIDE_MOD including the aa_list.
 */
void free_peptide_mod(PEPTIDE_MOD_T* mod){
  if(mod){
    free(mod);
  }
}

/* see below */
int generate_peptide_mod_list_TESTER(
 PEPTIDE_MOD_T*** peptide_mod_list,
 AA_MOD_T** aa_mod_list,
 int num_aa_mods);

/**
 * \brief Generate a list of all PEPTIDE_MODs that can be considered
 * given the list of AA_MODs and POSITION_MODs provided by the user
 * via the parameter file.
 *
 * This only needs to be called once for a search.  The list of
 * PEPTIDE_MODs can be reused for each spectrum.  Allows at most one
 * n-term and one c-term mod and at most mod->max_per_peptide
 * modifications of one type per peptide.  Sorts the list by number of
 * aa mods and returns only the number of peptide mods that have
 * between 0 and max-mods aa mods.
 *
 * The actual work is done in generate_peptide_mod_list_TESTER() so
 * that it can be tested without going through parameter.c.
 *
 * \returns The number of peptide_mods listed in the
 * peptide_mod_list argument.
 */
int generate_peptide_mod_list(
 PEPTIDE_MOD_T*** peptide_mod_list ///< return here the list of peptide_mods
 ){
  AA_MOD_T** aa_mod_list = NULL;
  int num_aa_mods = get_all_aa_mod_list(&aa_mod_list);

  // continued below
  int count = generate_peptide_mod_list_TESTER( peptide_mod_list, 
                                                aa_mod_list, num_aa_mods);
  return count; 
}

int generate_peptide_mod_list_TESTER(
 PEPTIDE_MOD_T*** peptide_mod_list,
 AA_MOD_T** aa_mod_list,
 int num_aa_mods)
{

  // initialize list of peptide mods with one unmodified peptide
  LINKED_LIST_T* final_list = new_list( new_peptide_mod() );
  int final_counter = 1;

  // for each non-fixed aa_mod
  int num_fixed_mods = get_num_fixed_mods();
  for(int mod_list_idx = num_fixed_mods; mod_list_idx < num_aa_mods; 
      mod_list_idx++){

    AA_MOD_T* cur_aa_mod = aa_mod_list[mod_list_idx];
    int cur_mod_max = aa_mod_get_max_per_peptide(cur_aa_mod);
    char cur_mod_id = aa_mod_get_symbol(cur_aa_mod);
    carp(CARP_DETAILED_DEBUG, "cur max %d, id %c", cur_mod_max, cur_mod_id);

    // initialize temp list 
    LINKED_LIST_T* temp_list = new_empty_list();
    int temp_counter = 0;

    int copies = 1;
    for(copies = 1; copies < cur_mod_max + 1; copies++){

      // for each pep_mod in final list
      LIST_POINTER_T* final_ptr = get_first_linked_list(final_list);
      while( final_ptr != NULL ){
        // get the current peptide mod
        PEPTIDE_MOD_T* cur_pep_mod = 
          (PEPTIDE_MOD_T*)get_data_linked_list(final_ptr);
        carp(CARP_DETAILED_DEBUG, "cur pep_mod has %d mods", 
             peptide_mod_get_num_aa_mods(cur_pep_mod)); 

        // copy the peptide_mod and add to it i copies of cur_mod
        PEPTIDE_MOD_T* mod_copy = copy_peptide_mod(cur_pep_mod);
        //peptide_mod_add_aa_mod(mod_cpy, cur_aa_mod, copies);
        peptide_mod_add_aa_mod(mod_copy, mod_list_idx, copies);
        carp(CARP_DETAILED_DEBUG, "adding %d %c's to temp", 
             copies, cur_mod_id);

        // add to temp list
        //temp_list_end = add_linked_list(temp_list_end, mod_cpy);
        push_back_linked_list(temp_list, mod_copy);

        // if this is the first added, set the head of temp list
        //temp_list = (temp_list==NULL) ? temp_list_end : temp_list;
        temp_counter++;
        final_ptr = get_next_linked_list(final_ptr);
      }// last in final list

    }// max copies

    // append temp list to final list
    carp(CARP_DETAILED_DEBUG, "adding temp list (%d) to final (%d)",
         temp_counter, final_counter);
    combine_lists(final_list, temp_list);
     free(temp_list);
    final_counter += temp_counter;
  }// last aa_mod

  // allocate an array of PEPTIDE_MOD_T* to return
  PEPTIDE_MOD_T** final_array = (PEPTIDE_MOD_T**)mycalloc(final_counter, 
                                                   sizeof(PEPTIDE_MOD_T*));
  // fill in the array and delete the list
  LIST_POINTER_T* final_list_ptr = get_first_linked_list(final_list);
  int pep_idx = 0;
  while(final_list_ptr != NULL){
    assert(pep_idx < final_counter);

    final_array[pep_idx] =(PEPTIDE_MOD_T*)get_data_linked_list(final_list_ptr);
    final_list_ptr = get_next_linked_list(final_list_ptr);
    pep_idx++;
  }

  delete_linked_list(final_list);

  // sort list
  qsort(final_array, final_counter, sizeof(PEPTIDE_MOD_T*),
        (QSORT_COMPARE_METHOD)compare_peptide_mod_num_aa_mods);

  // find the index of the first peptide mod with too many aa mods
  int max = get_int_parameter("max-mods");
  // should use binary search...
  int mod_idx = 0;
  for(mod_idx=0; mod_idx < final_counter; mod_idx++){
    if(final_array[mod_idx]->num_mods > max){
      break;
    }
  }

  carp(CARP_INFO, 
       "Created %d peptide mods, keeping %d with %d or fewer aa mods", 
       final_counter, mod_idx, max);
  // set return value
  *peptide_mod_list = final_array;
  //  return final_counter;
  return mod_idx;

  /*
    get the number of items in final list

    init final_list to have one unmodified pep
    init temp_list to empty

    for each mod in aa_mod_list
      for i 1 to mod_max
        for each pep_mod in final_list
          copy the pep_mod and add i aa_mods
          add to the temp_list
        end of final_list
      end i==max
      add the temp_list to the final_list
    last mod

    sort by num-mods, fewest to most
   */
}

/**
 * \brief Given the sequence, a modification, and the position for the sequence
 * modification, determine whether this modification is valid.
 * \returns true if the modification can be applied.
 */
bool is_sequence_position_modifiable(
  char* sequence, 
  int position, 
  AA_MOD_T* mod) {

  if(is_aa_modifiable(char_aa_to_modified(sequence[position]), mod)) {

    if (aa_mod_get_prevents_cleavage(mod)) {
      //make sure that position is not the beginning or end.
      if (position == 0 || sequence[position+1] == '\0') {
        return false;
      } else {
        return true;
      }
    } else {
      //modification doesn't prevent cleavage, return true.
      return true;
    }
  } else {
    //aa is not modifiable by this modificationa
    return false;
  }
}

/**
 * \brief Given the sequence, a modification, and the position for the sequence
 * modification, determine whether this modification is valid.
 * \returns true if the modification can be applied.
 */
bool is_sequence_position_modifiable(
  MODIFIED_AA_T* sequence,
  int position,
  AA_MOD_T* mod) {
  
  if( is_aa_modifiable( sequence[position], mod )){

    if (aa_mod_get_prevents_cleavage(mod)) {
      //make sure that position is not the beginning or end.
      if (position == 0 || sequence[position+1] == MOD_SEQ_NULL) {
        return false;
      } else {
        return true;
      }
    } else {
      //modification doesn't prevent cleavage, return true.
      return true;
    }
  } else {
    return false;
  }
}

/**
 * \brief Check a peptide sequence to see if the aa_mods in
 * peptide_mod can be applied. 
 *
 * Assumes that an amino acid can be modified by more than one aa_mod,
 * but not more than once by a sligle aa_mod.
 * \returns true if the sequence can be modified, else false
 */
bool is_peptide_modifiable
 (Peptide* peptide,          ///< The peptide to apply mods to
  PEPTIDE_MOD_T* peptide_mod){ ///< The mods to apply

  // a NULL peptide cannot be modified
  if( peptide == NULL ){
    return false;
  }
  // peptide mods with no aa mods can be applied to any peptide
  if( peptide_mod == NULL || 
      peptide_mod_get_num_aa_mods( peptide_mod ) == 0 ){
    return true;
  }

  char* sequence = peptide->getSequence();

  // for each aa_mod (skip those not in peptide_mod)
  
  AA_MOD_T** all_mods = NULL;
  //  int num_aa_mods = get_aa_mod_list(&all_mods);
  int num_aa_mods = get_all_aa_mod_list(&all_mods);
  assert( num_aa_mods < MAX_AA_MODS );

  bool success = true;
  for(int amod_idx = 0; amod_idx < num_aa_mods; amod_idx++){

    // FIRST: check to see if it is included in this pmod
    if( peptide_mod->aa_mod_counts[amod_idx] == 0 ){
      continue;
    }

    // SECOND: check modifyability based on type
    AA_MOD_T* cur_aa_mod = all_mods[amod_idx];
    // for position-based mods, check distance from protein end
    int max_distance = aa_mod_get_max_distance(cur_aa_mod);
    int locations_count = 0;
    int cur_seq_idx = 0;  
    switch( aa_mod_get_position(cur_aa_mod) ){
    case C_TERM: 
      if( max_distance < peptide->getCDistance()){ 
        success = false;
      }
      break;

    case N_TERM:
      if( max_distance < peptide->getNDistance()){ 
        success = false;
      }
      break;

      // count leagal locations for this aa mod, compare with counts
    case ANY_POSITION:
      // look for an aa in the seq where this mod can be placed
      while(sequence[cur_seq_idx] != '\0' ){
        if (is_sequence_position_modifiable(sequence, cur_seq_idx, cur_aa_mod)) {
            locations_count++;
        }// else keep looking
        cur_seq_idx++;
      }// end of sequence
      
      if( locations_count < peptide_mod->aa_mod_counts[amod_idx] ){
        //return false; 
        success = false;
      }
      break;

    }// end of switch

    if( success == false ){
      break; // done look at any more aa_mods
    }
  }// next in aa_mod list

  free(sequence);
  return success;
}

// move this to peptide.c
/*
void add_peptide_mod_seq(Peptide* peptide, MODIFIED_AA_T* cur_mod_seq){
  if( peptide == NULL || cur_mod_seq == NULL ){
    carp(CARP_ERROR, "Cannot add NULL modified sequence to null peptide");
  }

  // test out that the mod seq is valid
  int i = 0;
  while( cur_mod_seq[i] != MOD_SEQ_NULL){
    printf("%d(%c) ", cur_mod_seq[i], (char)cur_mod_seq[i] + 'A');
    i++;
  }
  printf("\n");
}
  */


/**
 * \brief Take a peptide and a peptide_mod and return via the third
 * arguement a list of modified peptides.
 *
 * The peptide_mod should be guaranteed to be applicable to
 * the peptide at least once.  However, there may not be any modified
 * forms that pass the max_aas_modified filter. For the peptide_mod to
 * be successfully 
 * applied, every aa_mod in its list must be applied to the sequence
 * as many times as its value in the aa_mods_counts array.  A single
 * amino acid can be modified multiple times by different aa_mods, but
 * only once for a given aa_mod (as defined by modifiable()).  The
 * modified peptides will be added to the end of the modified_peptides
 * list.  Any existing elements in the list will not be removed.
 * 
 * \returns The number of modified peptides in the list pointed to by
 * modified_peptides. 
 */
int modify_peptide(
  Peptide* peptide,             ///< the peptide to modify
  PEPTIDE_MOD_T* peptide_mod,     ///< the set of aa_mods to apply
  LINKED_LIST_T* modified_peptides,///< the returned modified peptides
  int max_aas_modified            ///< filter out peptides > m_a_m
){ 

  if( peptide == NULL ){
    carp(CARP_ERROR, "Cannot modify NULL peptide or use NULL peptide mod");
    return 0;
  }
  if( modified_peptides == NULL ){
    carp(CARP_ERROR, "Cannot return modified peptides to NULL list.");
    return 0;
  }
  // in case of no modifications
  if( peptide_mod == NULL || 
      peptide_mod_get_num_aa_mods(peptide_mod) == 0 ){
    carp(CARP_DETAILED_DEBUG, 
         "Modifying peptide with no aa_mods, return peptide copy");
    Peptide* peptide_copy = new Peptide(peptide); 
    push_back_linked_list(modified_peptides, peptide_copy);
    return 1;
  }

  // get the peptide sequence and convert to MODIFIED_AA_T*
  char* sequence = peptide->getSequence();
  MODIFIED_AA_T* pre_modified_seq = NULL;
  convert_to_mod_aa_seq(sequence, &pre_modified_seq);

  carp(CARP_DETAILED_DEBUG, "Modifying peptide %s", sequence);

  // get the aa_mod info
  int* aa_mod_counts = peptide_mod->aa_mod_counts;
  AA_MOD_T** aa_mod_list = NULL;
  int num_aa_mods = get_all_aa_mod_list(&aa_mod_list);

  // initialize a list of modified seqs with the unmod seq
  LINKED_LIST_T* modified_seqs = new_list(pre_modified_seq);
  int aa_mod_idx = 0;
  int total_count = 0;

  // try applying each aa_mod to the sequence
  for(aa_mod_idx = 0; aa_mod_idx < num_aa_mods; aa_mod_idx++){

    //printf("aaidx is %d and mod count is %d\n", aa_mod_idx, aa_mod_counts[aa_mod_idx]);
    int mod_count = aa_mod_counts[aa_mod_idx];
    if( mod_count == 0 ){ continue; } // do not apply this aa mod

    //fprintf(stderr, "applying to list, total count is %d\n", total_count);
    total_count = apply_mod_to_list(modified_seqs, 
                                    aa_mod_list[aa_mod_idx],
                                    mod_count,
                                    max_aas_modified);

    //printf("after applying count is %d\n", total_count);
    // the count should be > 0, but check for error case
    if( total_count == 0 || is_empty_linked_list(modified_seqs) ){
      carp(CARP_DETAILED_DEBUG,
           "Peptide modification could not be applied to sequence %s",
           sequence);
      return total_count;  // free things...
    }
  } // next aa_mod

  carp(CARP_DETAILED_DEBUG, "Sequence %s has %i modified forms.", sequence, total_count);
  free(sequence);

  // create a peptide for each sequence and add it to the list   
  // filter out those with more than max_aas_modified
  total_count = 0;
  while( ! is_empty_linked_list( modified_seqs ) ){

    MODIFIED_AA_T* cur_mod_seq = 
      (MODIFIED_AA_T*)pop_front_linked_list(modified_seqs);
    if( count_modified_aas(cur_mod_seq) > max_aas_modified ){
      continue;
    }

    Peptide* cur_peptide = new Peptide(peptide);
    cur_peptide->setMod(cur_mod_seq, peptide_mod);

    push_back_linked_list(modified_peptides, cur_peptide );
    total_count++;
  }
  carp(CARP_DETAILED_DEBUG, "There were %i modified seqs created", total_count);
  free(modified_seqs);
  return total_count;
}


/**
 * \brief Beginning with a list of modified seqs, apply another
 * modification to each in as many ways as possible.  Update the list
 * to contain the new seqs.
 * 
 * To be used for adding multiple modifications to a peptide
 * sequence.  Adding the first aa mod creates one list of seqs.  The
 * second aa mod is then applied to the list.  Any seqs that can have
 * both modifications applied is returned and those that cannot are
 * removed from the list.
 *
 * \returns The number of modified seqs in the list.
 */
int apply_mod_to_list(
  LINKED_LIST_T* apply_mod_to_these, ///< a pointer to a list of seqs
  AA_MOD_T* mod_to_apply,  ///< the specific mod to apply
  int num_to_apply,         ///< how many of this mod to apply
  int max_aas_modified ///< don't return seqs with more than this many aas moded
){

  /* Null cases */
  if( is_empty_linked_list(apply_mod_to_these) ){ // nothing to apply mod to
    return 0;
  }

  if( mod_to_apply == NULL ){  // no modification to add, return as is
  // get the current count
    LIST_POINTER_T* cur_seq = get_first_linked_list(apply_mod_to_these);
    int num_seqs = 0;
    while( cur_seq != NULL ){
      num_seqs++;
    }
    return num_seqs;
  }

  //printf("iterating over seqs in list, applying %d of mod %c\n", 
  //       num_to_apply, aa_mod_get_symbol(mod_to_apply));

  // initialize a temp list of seqs to return
  LINKED_LIST_T* completed_seqs = new_empty_list();

  int return_seq_count = 0;
  int times_idx = 0;

  // apply the mod num_to_apply times
  for(times_idx = 0; times_idx < num_to_apply; times_idx++){
    // reset count
    return_seq_count = 0;

    // for each seq
    while( ! is_empty_linked_list(apply_mod_to_these) ){
      // get modified seq
      MODIFIED_AA_T* cur_seq = 
        (MODIFIED_AA_T*)pop_front_linked_list(apply_mod_to_these);
      
      /*
      char* print_me = modified_aa_string_to_string(cur_seq);
      printf("Original: %s\n", print_me);
      free(print_me);
      */
      return_seq_count += apply_mod_to_seq(cur_seq,
                                           mod_to_apply,
                                           times_idx, // skip n modified aas
                                           completed_seqs);

      //printf("Now there are %d seqs to return\n", return_seq_count);
      free(cur_seq);
    }// apply to next seq until list is empty

    // make the collected results the input to the next application
    combine_lists(apply_mod_to_these, completed_seqs); 
    free(completed_seqs); // just deletes head, not full list
    completed_seqs = new_empty_list();

  }// apply next time

  free( completed_seqs ); // just the head of the list

  // Check all seqs in list and remove those with too many aas modified
  // push a dummy onto the front for ease of deleting
  // in a list nodeA->nodeB, we need nodeA to delete nodeB
  // call A and B prev and cur
  push_front_linked_list(apply_mod_to_these, NULL);
  LIST_POINTER_T* prev_node = get_first_linked_list(apply_mod_to_these);
  LIST_POINTER_T* cur_node = NULL;
  MODIFIED_AA_T* cur_seq = NULL;

  while( has_next_linked_list(prev_node) ){
    cur_node = get_next_linked_list(prev_node);
    cur_seq = (MODIFIED_AA_T*)get_data_linked_list(cur_node);

    if( count_modified_aas(cur_seq) > max_aas_modified){
      // delete cur, set cur to node after prev
      cur_node = delete_next_list_node(prev_node);
      free(cur_seq);  // list delete doesn't remove data

    }else
      prev_node = cur_node;
      cur_node = get_next_linked_list(prev_node);
      cur_seq = (MODIFIED_AA_T*)get_data_linked_list(cur_node);
  }// next seq
  // pop the dummy off the front
  pop_front_linked_list(apply_mod_to_these);

  return return_seq_count;
}

/**
 * \brief Beginning with a peptide sequence (of modified aas) and an
 * aa_mod, return via the 'return_list' parameter a list of modified
 * versions of the peptide sequence.
 *
 * The initial seq is left untouched.  If the seq cannot be modified
 * by the aa_mod, return_list is unchanged and 0 is returned.  If the
 * mod is NULL the return list is given one item, an unchanged copy of
 * the seq.
 *
 * \returns The number of ways in which the sequence could be
 * modified, i.e. the number of items added to return_list.
 */
// TODO (BF 17-Sep-08): could add seq length to make consistent with
// other functions with MODIFIED_AA_T array parameters
int apply_mod_to_seq(
  MODIFIED_AA_T* seq,          ///< the seq to modify
  AA_MOD_T* mod,               ///< the mod to apply
  int skip_n,                  ///< pass over n aas modified by this mod
  LINKED_LIST_T* return_list){ ///< the newly modified versions of the seq

  // Boundry cases
  if( seq == NULL ){
    return 0;
  }
  if( return_list == NULL ){
    carp(CARP_ERROR, "Cannot add modified sequences to a NULL list.");
    return 0;
  }

  // find seq length (see TODO above)
  int seq_len = 0;
  while( seq[seq_len] != MOD_SEQ_NULL ){
    seq_len++;
  }

  // no modification to add, return seq as is
  if( mod == NULL ){
    push_back_linked_list( return_list, copy_mod_aa_seq(seq, seq_len) );
    return 1;
  }

  // special case: positional mods
  if( N_TERM == aa_mod_get_position(mod) ){
    //MODIFIED_AA_T* seq_copy = copy_mod_aa_seq(seq);
    MODIFIED_AA_T* seq_copy = copy_mod_aa_seq(seq, seq_len);
    modify_aa( &seq_copy[0], mod );
    push_back_linked_list(return_list, seq_copy);
    return 1;  // one added to list
  }
  if( C_TERM == aa_mod_get_position(mod)){
    //MODIFIED_AA_T* seq_copy = copy_mod_aa_seq(seq);
    MODIFIED_AA_T* seq_copy = copy_mod_aa_seq(seq, seq_len);
    /*
    // find last index
    int seq_idx = 0;
    while( seq[seq_idx] != MOD_SEQ_NULL ){
      seq_idx++;
    }
    */
    //modify_aa( &seq_copy[seq_idx-1], mod );
    modify_aa( &seq_copy[seq_len-1], mod );
    push_back_linked_list(return_list, seq_copy);
    return 1;  // one added to list
  }

  // first loop, skip over n modified aas in the seq
  int seq_idx = 0;
  while( seq[seq_idx] != MOD_SEQ_NULL  && skip_n != 0){
    if( is_aa_modified( seq[seq_idx], mod ) ){
      skip_n--;
    }
    seq_idx++;
  }// nex mod_aa; stop at end of seq or having skipped n

  // second loop, for each aa: check for modifiability, copy and add

  // check each aa
  // copy seq when modifiable and add to list
  int count = 0;
  while( seq[seq_idx] != MOD_SEQ_NULL ){
    if (is_sequence_position_modifiable(seq, seq_idx, mod)) {
      //      MODIFIED_AA_T* seq_copy = copy_mod_aa_seq(seq);
      MODIFIED_AA_T* seq_copy = copy_mod_aa_seq(seq, seq_len);
      modify_aa( &seq_copy[seq_idx], mod );

      push_back_linked_list(return_list, seq_copy);
      count++;
    }
    
    seq_idx++;
  }

  return count;  
}

/**
 * print all fields in peptide mod. For debugging
 */
void print_p_mod(PEPTIDE_MOD_T* mod){
  printf("PMOD: mass %.2f, num %d, aa mods\n  ", 
         mod->mass_change, mod->num_mods);

  int i=0;
  AA_MOD_T** all_mods = NULL;
  //int num_aa_mods = get_aa_mod_list( &all_mods );
  int num_aa_mods = get_all_aa_mod_list( &all_mods );
  assert( num_aa_mods < MAX_AA_MODS );

  for(i = 0; i < num_aa_mods; i++){
    print_a_mod( all_mods[i] );
    printf("  ");
  }
  printf("\n");
}

/* Setters and Getters */
/**
 * \brief Add a new aa_mod to the peptide mod.  Updates mass_change,
 * num_mods and list_of_aa_mods.  Does not enforce the copy number of
 * an aa_mod to be less than max_per_peptide.
 * \returns void
 */
void peptide_mod_add_aa_mod(
  PEPTIDE_MOD_T* pep_mod, ///< The peptide mod being added to
  int aa_mod_idx,      ///< The index into the global list of aa mods
  int copies ){           ///< How many of the aa_mod to add

  AA_MOD_T** all_mods = NULL;
  //  int num_aa_mods = get_aa_mod_list( &all_mods );
  int num_aa_mods = get_all_aa_mod_list( &all_mods );
  assert( num_aa_mods < MAX_AA_MODS );

  //printf("mass change for amod of index %d is %.2f\n", aa_mod_idx, aa_mod_get_mass_change(all_mods[aa_mod_idx]));
  pep_mod->mass_change += aa_mod_get_mass_change(all_mods[aa_mod_idx])
                            * copies;
  //printf("pep mod mass change is now %.2f\n", pep_mod->mass_change );
  pep_mod->num_mods += copies;

  pep_mod->aa_mod_counts[aa_mod_idx] += copies;
}

/**
 * \brief Get the value of the net mass change for this peptide_mod.
 * \returns The mass change for the peptide mod.
 */
double peptide_mod_get_mass_change(PEPTIDE_MOD_T* mod){
  return mod->mass_change;
}

/**
 * \brief Get the number of aa_mods in this peptide_mod.
 * \returns The number of aa_mods in this peptide_mod.
 */
int peptide_mod_get_num_aa_mods(PEPTIDE_MOD_T* mod){
  return mod->num_mods;
}


/**
 * \brief Compares the number of aa mods in two peptide mods for
 * sorting.
 * \returns Negative int, 0, or positive int if the number of aa_mods
 * in pmod 1 is less than, equal to or greater than the number of
 * aa_mods in pmod2, respectively.
 */
int compare_peptide_mod_num_aa_mods(const void* pmod1, 
                                    const void* pmod2){
  return (*(PEPTIDE_MOD_T**)pmod1)->num_mods 
            - (*(PEPTIDE_MOD_T**)pmod2)->num_mods;
}
































  /* 
    initialize a counter to 0
  last aa_mod

  for each aa in peptide sequence
    for each aa_mod
      counter[aa_mod] += aa_mod.aa_list[aa]
    last aa_mod 
  last aa

  bool success = true
  for each aa_mod in peptide_mod
    success && (counter[aa_mod] >= peptide_mod.occurs[aa_mod])
  last aa_mod

  return success
  */

