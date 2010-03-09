/*************************************************************************//**
 * \file generate_peptides_iterator.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: Nov 8 2007
 * DESCRIPTION: object to return candidate peptides with a given restriction
 * REVISION: 
 ****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include "carp.h"
#include "peptide.h"
#include "peptide_src.h"
#include "protein.h"
#include "database.h"
#include "parameter.h"
#include "index.h"
#include "generate_peptides_iterator.h"

/**
 * \struct generate_peptides_iterator_t
 * \brief An object that navigates the options and selects the correct 
 * peptide iterator to use 
 */
struct generate_peptides_iterator_t{
  void* iterator;     ///< the index or database iterator we are wrapping 
  BOOLEAN_T (*has_next)(void*); ///< the function pointer to *_has_next
  PEPTIDE_T* (*next)(void*);    ///< the function pointer to *_next
  void (*free)(void*);          ///< the function pointer to *_free
  INDEX_T* index;               ///< the index object needed
  DATABASE_T* database;         ///< the database object needed
  PEPTIDE_CONSTRAINT_T* constraint; ///< peptide constraint
};
//do index,database, and constraint need to be members since they are
//not used in has_next, get_next, or free?

/**
 * \returns a empty generate_peptides_iterator object
 */
GENERATE_PEPTIDES_ITERATOR_T* allocate_generate_peptides_iterator(){
  // allocate an empty iterator
  GENERATE_PEPTIDES_ITERATOR_T* iterator = (GENERATE_PEPTIDES_ITERATOR_T*)
    mycalloc(1, sizeof(GENERATE_PEPTIDES_ITERATOR_T));

  // TODO (BF 10-Apr-08): set member fields to NULL
  return iterator;
}

/**
 * \brief Create a peptide iterator based entirely on parameter values
 * (defaults and those given by user).
 *
 * The peptides are drawn from either a fasta file or an index based
 * on the value of the use-index parameter.  With default values, mass
 * range is wide so this effectively iterates over "all" peptides in
 * the protein input.  If no peptides in the protein source meet the
 * criteria, a peptide iterator is still returned, but when passed to
 * *_has_next() it will always return FALSE.
 *\returns A new generate_peptide_iterator object
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator(void){
  // get parameters from parameter.c
  double min_mass = get_double_parameter("min-mass");
  double max_mass = get_double_parameter("max-mass");
  //  BOOLEAN_T use_index = get_boolean_parameter("use-index");

  //  BOOLEAN_T is_unique = get_boolean_parameter("unique-peptides");
  const char*  protein_input_name = get_string_parameter_pointer("protein input");
  BOOLEAN_T use_index = is_directory(protein_input_name);

  INDEX_T* index = NULL;
  DATABASE_T* database = NULL;
  // TODO (BF 27-Feb-08): set use_index according to the input, true if dir
  if (use_index == TRUE){
    //index = new_index_from_disk(protein_input_name, is_unique);
    index = new_index_from_disk(protein_input_name);
  } else {
    // FALSE indicates that we are not using a binary fasta file
    database = new_database(protein_input_name, FALSE);
  }

  return new_generate_peptides_iterator_from_mass_range(min_mass, max_mass, 
      index, database);
}

/**
 * \brief Create a peptide iterator for peptides of a target mass.
 *
 * Peptides with mass between neutral_mass +/- "mass-window"
 * (a user-defined parameter).  Peptides are drawn from either the
 * given index or the given database, if the index is NULL. If no
 * peptides in the protein source meet the criteria, a peptide
 * iterator is still returned, but when passed to *_has_next() it will
 * always return FALSE.  
 *
 *\returns A new generate_peptide_iterator object
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator_from_mass(
  FLOAT_T neutral_mass, ///< The target mass (uncharged) for peptides
  INDEX_T* index,     ///< The index from which to draw peptides OR
  DATABASE_T* database///< The database from which to draw peptides
  )
{
  // get parameters
  double mass_window = get_double_parameter("precursor-window");
  double min_mass = neutral_mass - mass_window;
  double max_mass = neutral_mass + mass_window;

  carp(CARP_DETAILED_DEBUG,"Generating peptides in %.2f ~ %.2f", 
       min_mass, max_mass); 

  return new_generate_peptides_iterator_from_mass_range(min_mass, max_mass, 
      index, database);
}

/**
 * \brief Create a peptide iterator for peptides in a specific range
 * of masses.
 *
 * This is the version of new_* that is called by the others and does
 * all the work of setting the member variable fields.  Parameters
 * other than min and max mass are taken from parameter.c.  If no
 * peptides in the protein source meet the criteria, a peptide
 * iterator is still returned, but when passed to *_has_next() it will
 * always return FALSE.  
 *
 * \returns a new generate_peptides_iterator object
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator_from_mass_range(
  double min_mass,     ///< The min mass of peptides to generate -in
  double max_mass,     ///< The maximum mas of peptide to generate -in
  INDEX_T* index,      ///< The index
  DATABASE_T* database ///< The database
  )
{
  // get parameters
  int min_length = get_int_parameter("min-length");
  int max_length = get_int_parameter("max-length");
  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  DIGEST_T digestion = get_digest_type_parameter("digestion");
  MASS_TYPE_T mass_type = get_mass_type_parameter("isotopic-mass");
  BOOLEAN_T missed_cleavages = get_boolean_parameter("missed-cleavages");
  SORT_TYPE_T sort_type = get_sort_type_parameter("sort");
  BOOLEAN_T is_unique = get_boolean_parameter("unique-peptides");

    
  // allocate an empty iterator
  GENERATE_PEPTIDES_ITERATOR_T* gen_peptide_iterator 
    = allocate_generate_peptides_iterator();
  
  // peptide constraint
  PEPTIDE_CONSTRAINT_T* constraint 
    //= new_peptide_constraint(peptide_type, min_mass, max_mass, 
    = new_peptide_constraint(enzyme, digestion, min_mass, max_mass, 
        min_length, max_length, missed_cleavages, mass_type);
  
  // assign to iterator
  gen_peptide_iterator->constraint = copy_peptide_constraint_ptr(constraint); 

  // Check that index OR database exists
  if (database == NULL && index == NULL ){
    carp(CARP_FATAL, "Cannot genrate peptides when index and database are both NULL.");
  }

  /***********************
   * use index file
   **********************/
  if(index != NULL){
    
    if((sort_type != SORT_MASS && sort_type != SORT_NONE)){
      carp(CARP_FATAL, "Cannot sort other than by mass when using index.");
    }
   
    // already tested for 
    if(index == NULL){
      carp(CARP_FATAL, "Failed to create peptides from index");
    } 

    // use array implementation of peptide_src
    set_peptide_src_implementation(FALSE);

    // create index and set to generate_peptides_iterator
    //    set_index_constraint(index, constraint); 
    set_index_search_constraint(index, constraint); 
    gen_peptide_iterator->index = copy_index_ptr(index);
    
    // only resrict peptide by mass and length, default iterator
    //if(peptide_type == ANY_TRYPTIC){ //BF: == ALL? 
    if(digestion == NON_SPECIFIC_DIGEST){
      // create index peptide interator & set generate_peptides_iterator
      INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator
        = new_index_peptide_iterator(index);        
      gen_peptide_iterator->iterator = index_peptide_iterator;
      gen_peptide_iterator->has_next = &void_index_peptide_iterator_has_next;
      gen_peptide_iterator->next = &void_index_peptide_iterator_next;
      gen_peptide_iterator->free = &void_free_index_peptide_iterator;
    }
    // if need to select among peptides by peptide_type and etc.
    // if need to select among peptides by enzyme
    else{

      INDEX_FILTERED_PEPTIDE_ITERATOR_T* index_filtered_peptide_iterator 
        = new_index_filtered_peptide_iterator(index);
      gen_peptide_iterator->iterator = index_filtered_peptide_iterator;
      gen_peptide_iterator->has_next 
        = &void_index_filtered_peptide_iterator_has_next;
      gen_peptide_iterator->next = &void_index_filtered_peptide_iterator_next;
      gen_peptide_iterator->free = &void_free_index_filtered_peptide_iterator;
    }
  }

  /*********************************************
   * Read in from fasta file, don't use index
   ************************************************/
  else{

    // def used for each iterator

    // set for all peptide src use link list implementation
    // this routine sets the static global in peptide.c
    set_peptide_src_implementation(TRUE);

    // create a new database & set generate_peptides_iterator
    gen_peptide_iterator->database = copy_database_ptr(database);
    
    // no sort
    //if(!is_unique && sort_type == NONE){ 
    if( sort_type == SORT_NONE ){ 
      carp(CARP_DETAILED_DEBUG, "Creating database peptide iterator");
      // create peptide iterator  & set generate_peptides_iterator
      DATABASE_PEPTIDE_ITERATOR_T* iterator 
        = new_database_peptide_iterator(database, constraint); 
      gen_peptide_iterator->iterator = iterator;
      gen_peptide_iterator->has_next = &void_database_peptide_iterator_has_next;
      gen_peptide_iterator->next = &void_database_peptide_iterator_next;
      gen_peptide_iterator->free = &void_free_database_peptide_iterator;
      
    }      
    // sort or check for unique
    else{   // should only be used for generate-peptides
      carp(CARP_DETAILED_DEBUG, "Creating sorted database peptide iterator");
      // only sort, by default will be sorted by mass
      DATABASE_SORTED_PEPTIDE_ITERATOR_T* sorted_iterator = NULL;
      if(sort_type == SORT_NONE){
        // create peptide iterator
        sorted_iterator = new_database_sorted_peptide_iterator(
            database, constraint, 
            SORT_MASS, TRUE);       
      }
      // create peptide iterator
      else{
        sorted_iterator = new_database_sorted_peptide_iterator(
            database, constraint, 
            sort_type, is_unique);
      }
      
      // set generate_peptides_iterator 
      gen_peptide_iterator->iterator = sorted_iterator;
      gen_peptide_iterator->has_next 
        = &void_database_sorted_peptide_iterator_has_next;
      gen_peptide_iterator->next = &void_database_sorted_peptide_iterator_next;
      gen_peptide_iterator->free = &void_free_database_sorted_peptide_iterator;
    }
  }
  free_peptide_constraint(constraint);
  return gen_peptide_iterator;
}

/**
 * \brief Create a new peptide iterator to return modified peptides.
 *
 * Peptides are generated based on values in parameter.c with the
 * exception of the mass range.  Min and max mass are given by target
 * mass +/- mass-window.  Only those peptides that can be modified by
 * the peptide_mod are returned.  (A peptide_mod may contain no
 * AA_MODS, in which case all peptides are returned.)  Peptides are
 * taken either from an index or from a database (fasta file).  If no
 * peptides pass the criteria, a new iterator is still returned, but
 * when passed to has_next() it will always return FALSE.
 *
 * \returns A newly allocated peptide iterator.
 */
/*
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator_mods(
  double mass,                ///< target mass of peptides
  PEPTIDE_MOD_T* pmod,        ///< the peptide mod to apply
  INDEX_T* index,             ///< index from which to draw peptides OR
  DATABASE_T* dbase){         ///< database from which to draw peptides
  
  // allocate a new iterator
  GENERATE_PEPTIDES_ITERATOR_T* new_iterator =
    allocate_generate_peptides_iterator();

  // create a mod_pept-iter and point to it (mass, pmod, index, dbase)
  MODIFIED_PEPTIDES_ITERATOR_T* mod_peptide_iterator = 
    new_modified_peptides_iterator_from_mass( mass, pmod, index, dbase );
  new_iterator->iterator = mod_peptide_iterator;

  // set has_next, get_next, free
  new_iterator->has_next = &void_modified_peptides_iterator_has_next;
  new_iterator->next     = &void_modified_peptides_iterator_next;
  new_iterator->free     = &void_modified_peptides_iterator_free;

  return new_iterator;
}
*/

/****************************************************************************/

/**
 * \brief Establish if an iterator has more peptides to return.
 *
 * Typically, the iterator will have next_peptide queued up and ready
 * to return.  But do not trust that the new_*_iterator function
 * initialized it correctly.  If next_peptdide is NULL, first check
 * that has_next is also FALSE.  Initialize next_peptide if necessary.
 * \returns TRUE, if there is a next peptide, else FALSE
 */
BOOLEAN_T generate_peptides_iterator_has_next(
 GENERATE_PEPTIDES_ITERATOR_T* generate_peptides_iterator ///< working iterator
  )
{
  return generate_peptides_iterator->has_next(generate_peptides_iterator->iterator);
}

/**
 *\returns the next peptide in the iterator
 */
PEPTIDE_T* generate_peptides_iterator_next(
 GENERATE_PEPTIDES_ITERATOR_T* generate_peptides_iterator ///< working iterator
  )
{
  return generate_peptides_iterator->next(generate_peptides_iterator->iterator);
}

/**
 * Don't free the iterator until completed with the peptides generated
 * Frees an allocated generate_peptides_iterator object
 */
void free_generate_peptides_iterator(
  GENERATE_PEPTIDES_ITERATOR_T* generate_peptides_iterator ///<iterator to free
  )
{
  // free the nested iterator
  generate_peptides_iterator->free(generate_peptides_iterator->iterator);
  free_peptide_constraint(generate_peptides_iterator->constraint);
  free(generate_peptides_iterator);
}
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
