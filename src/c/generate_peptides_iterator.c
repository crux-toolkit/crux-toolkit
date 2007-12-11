/*****************************************************************************
 * \file generate_peptides_iterator
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
  void* iterator;     ///< the object iterator that will be selected
  BOOLEAN_T (*has_next)(void*);     ///< the function pointer to *_has_next
  PEPTIDE_T* (*next)(void*);         ///< the function pointer to *_next
  void (*free)(void*);         ///< the function pointer to *_free
  INDEX_T* index;  ///< the index object needed
  DATABASE_T* database; ///< the database object needed
  PEPTIDE_CONSTRAINT_T* constraint; ///< peptide constraint
};


/**
 * \returns a empty generate_peptides_iterator object
 */
GENERATE_PEPTIDES_ITERATOR_T* allocate_generate_peptides_iterator(){
  // allocate an empty iterator
  GENERATE_PEPTIDES_ITERATOR_T* iterator = (GENERATE_PEPTIDES_ITERATOR_T*)
    mycalloc(1, sizeof(GENERATE_PEPTIDES_ITERATOR_T));

  return iterator;
}

/**
 * \returns a new generate_peptides_iterator object, with fasta file input
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator_from_mass_range(
  double min_mass,  ///< the min mass of peptides to generate -in
  double max_mass,  ///< the maximum mas of peptide to generate -in
  INDEX_T* index, ///< the index
  DATABASE_T* database ///< the database
  )
{
  // get parameters
  int min_length = get_int_parameter("min-length");
  int max_length = get_int_parameter("max-length");
  //char* isotopic_mass = get_string_parameter_pointer("isotopic-mass");
  //char* sort = get_string_parameter_pointer("sort");  // sort order
  BOOLEAN_T use_index_boolean = get_boolean_parameter("use-index");
  PEPTIDE_TYPE_T peptide_type = get_peptide_type_parameter("cleavages");

  // TODO put in parameter retrieval routines
  MASS_TYPE_T mass_type = get_mass_type_parameter("isotopic-mass");
  BOOLEAN_T missed_cleavages = get_boolean_parameter("missed-cleavages");
  //SORT_TYPE_T sort_type = NONE;
  //string_to_sort_type(sort, &sort_type);
  SORT_TYPE_T sort_type = get_sort_type_parameter("sort");
  BOOLEAN_T is_unique = get_boolean_parameter("unique-peptides");

  // check if maximum length is with in range <= 255
  /*  This is already checked
  if(max_length > 255){
    carp(CARP_FATAL, "maximum length:%d over limit 255.", max_length);
    exit(1);
  }
  */
  /*
  // determine isotopic mass option
  if(strcmp(isotopic_mass, "average")==0){
    mass_type = AVERAGE;
  }
  else if(strcmp(isotopic_mass, "mono")==0){
    mass_type = MONO;
  }
  else{
    carp(CARP_ERROR, "Incorrect argument %s", isotopic_mass);
  }
  */
  
  // determine sort type option
  /*
  if(strcmp(sort, "mass")==0){
    sort_type = MASS;
  }
  else if(strcmp(sort, "length")==0){
    sort_type = LENGTH;
  }
  else if(strcmp(sort, "lexical")==0){
    sort_type = LEXICAL;
  }
  else if(strcmp(sort, "none")==0){
    sort_type = NONE;
  }
  else{
    carp(CARP_ERROR, "Incorrect argument %s, using default value", sort);
  }
  */
    
  // allocate an empty iterator
  GENERATE_PEPTIDES_ITERATOR_T* gen_peptide_iterator 
    = allocate_generate_peptides_iterator();
  
  // peptide constraint
  PEPTIDE_CONSTRAINT_T* constraint 
    = new_peptide_constraint(peptide_type, min_mass, max_mass, 
        min_length, max_length, missed_cleavages, mass_type);
  
  // assign to iterator
  gen_peptide_iterator->constraint = copy_peptide_constraint_ptr(constraint); 

  /***********************
   * use index file
   **********************/
  if(use_index_boolean){
    
    carp(CARP_INFO, "Using index for peptide generation");

    if((sort_type != MASS && sort_type != NONE)){
      carp(CARP_FATAL, "Cannot sort other than by mass when using index.");
      exit(1);
    }
   
    if(index == NULL){
      carp(CARP_FATAL, "Failed to create peptides from index");
      free(gen_peptide_iterator);
      exit(1);
    }

    // use array implementation of peptide_src
    set_peptide_src_implementation(FALSE);

    // create index and set to generate_peptides_iterator
    set_index_constraint(index, constraint); 
    gen_peptide_iterator->index = copy_index_ptr(index);
    
    // only resrict peptide by mass and length, default iterator
    if(peptide_type == ANY_TRYPTIC){
      // create index peptide interator & set generate_peptides_iterator
      INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator
        = new_index_peptide_iterator(index);        
      gen_peptide_iterator->iterator = index_peptide_iterator;
      gen_peptide_iterator->has_next = &void_index_peptide_iterator_has_next;
      gen_peptide_iterator->next = &void_index_peptide_iterator_next;
      gen_peptide_iterator->free = &void_free_index_peptide_iterator;
    }
    // if need to select among peptides by peptide_type and etc.
    else{
      carp(CARP_INFO, "using filtered index peptide generation");
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
   
    carp(CARP_INFO, "Using fasta file for peptide generation");

    // set for all peptide src use link list implementation
    // this routine sets the static global in peptide.c
    set_peptide_src_implementation(TRUE);

    // create a new database & set generate_peptides_iterator
    gen_peptide_iterator->database = copy_database_ptr(database);
    
    // no sort, redundant
    if(!is_unique && sort_type == NONE){ 
      // create peptide iterator  & set generate_peptides_iterator
      DATABASE_PEPTIDE_ITERATOR_T* iterator 
        = new_database_peptide_iterator(database, constraint); 
      gen_peptide_iterator->iterator = iterator;
      gen_peptide_iterator->has_next = &void_database_peptide_iterator_has_next;
      gen_peptide_iterator->next = &void_database_peptide_iterator_next;
      gen_peptide_iterator->free = &void_free_database_peptide_iterator;
      
    }      
    // sort or check for unique
    else{
      // only sort, by default will be sorted by mass
      DATABASE_SORTED_PEPTIDE_ITERATOR_T* sorted_iterator = NULL;
      if(sort_type == NONE){
        // create peptide iterator
        sorted_iterator = new_database_sorted_peptide_iterator(
            database, constraint, 
            MASS, TRUE);       
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
 *\returns a new generate_peptide_iterator object with custom min, max mass for SP
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator_from_mass(
  float neutral_mass, ///< the neutral_mass that which the peptides will be searched -in
  INDEX_T* index,
  DATABASE_T* database
  )
{
  // get parameters
  double mass_window = get_double_parameter("mass-window");
  double min_mass = neutral_mass - mass_window;
  double max_mass = neutral_mass + mass_window;

  carp(CARP_DEBUG,"searching peptide in %.2f ~ %.2f", min_mass, max_mass); 

  return new_generate_peptides_iterator_from_mass_range(min_mass, max_mass, 
      index, database);
}

/**
 * \returns a new generate_peptides_iterator object from all parameters 
 * from parameters.c hash structure
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator(void){
  // get parameters from parameter.c
  double min_mass = get_double_parameter("min-mass");
  double max_mass = get_double_parameter("max-mass");
  BOOLEAN_T use_index = get_boolean_parameter("use-index");

  BOOLEAN_T is_unique = get_boolean_parameter("unique-peptides");
  //  char*  fasta_file = get_string_parameter_pointer("fasta-file");
  char*  fasta_file = get_string_parameter_pointer("protein input");

  INDEX_T* index = NULL;
  DATABASE_T* database = NULL;
  if (use_index == TRUE){
    index = new_index_from_disk(fasta_file, is_unique);
  } else {
    // FALSE indicates that we are not using a binary fasta file
    database = new_database(fasta_file, FALSE);
  }

  return new_generate_peptides_iterator_from_mass_range(min_mass, max_mass, 
      index, database);
}

/**********************************************************************************/

/**
 *\returns TRUE, if there is a next peptide, else FALSE
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
  GENERATE_PEPTIDES_ITERATOR_T* generate_peptides_iterator ///< iterator to free
  )
{
  // free the nested iterator
  generate_peptides_iterator->free(generate_peptides_iterator->iterator);
  free(generate_peptides_iterator);
}
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
