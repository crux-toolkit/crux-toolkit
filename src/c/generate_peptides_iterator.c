  
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
 *\struct generate_peptides_iterator
 *\brief An object that navigates the options and selects the correct peptide iterator to use 
 */
struct generate_peptides_iterator_t{
  void* iterator;     ///< the object iterator that will be selected
  BOOLEAN_T (*has_next)(void*);     ///< the function pointer to *_has_next
  PEPTIDE_T* (*next)(void*);         ///< the function pointer to *_next
  void (*free)(void*);         ///< the function pointer to *_free
  void (*free_peptide)(PEPTIDE_T*); ///< the function pointer to *_free_peptide
  INDEX_T* index;  ///< the index object needed
  DATABASE_T* database; ///< the database object needed
  PEPTIDE_CONSTRAINT_T* constraint; ///< peptide constraint
};


/**
 *\returns a empty generate_peptides_iterator object
 */
GENERATE_PEPTIDES_ITERATOR_T* allocate_generate_peptides_iterator(){
  //allocate an empty iterator
  GENERATE_PEPTIDES_ITERATOR_T* iterator = 
    (GENERATE_PEPTIDES_ITERATOR_T*)mycalloc(1, sizeof(GENERATE_PEPTIDES_ITERATOR_T));

  return iterator;
}

/**
 *\returns a new generate_peptides_iterator object
 */
GENERATE_PEPTIDES_ITERATOR_T* new_generate_peptides_iterator(){
  //get perameters
  double min_mass = get_double_parameter("min-mass", 200);
  double max_mass = get_double_parameter("max-mass", 2400);
  int min_length = get_int_parameter("min-length", 6);
  int max_length = get_int_parameter("max-length", 50);
  char* cleavages = get_string_parameter_pointer("cleavages");
  char* isotopic_mass = get_string_parameter_pointer("isotopic-mass");
  char* redundancy = get_string_parameter_pointer("redundancy");
  char* use_index = get_string_parameter_pointer("use-index");
  char* sort = get_string_parameter_pointer("sort");   // mass, length, lexical, none  

  BOOLEAN_T use_index_boolean = FALSE;
  MASS_TYPE_T mass_type = AVERAGE;
  PEPTIDE_TYPE_T peptide_type = TRYPTIC;
  BOOLEAN_T missed_cleavages = get_boolean_parameter("missed-cleavages", FALSE);
  char* in_file = get_string_parameter_pointer("fasta-file");
  BOOLEAN_T is_unique = FALSE;
  SORT_TYPE_T sort_type = NONE;

  //def used for each iterator
  DATABASE_PEPTIDE_ITERATOR_T* iterator = NULL;
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* sorted_iterator = NULL;
  DATABASE_T* database = NULL;
  INDEX_T* index = NULL;
  INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator = NULL;
  INDEX_FILTERED_PEPTIDE_ITERATOR_T* index_filtered_peptide_iterator = NULL;
  
  //parse all the necessary parameters
  //FIXME may add additional types such as non-trypticc or partially-tryptic
  if(strcmp(cleavages, "all")==0){
    peptide_type = ANY_TRYPTIC;
  }
  else if(strcmp(cleavages, "tryptic")==0){
    peptide_type = TRYPTIC;
  }
  else if(strcmp(cleavages, "partial")==0){
    peptide_type = PARTIALLY_TRYPTIC;
  }
  else{
    carp(CARP_ERROR, "incorrect argument %s, using default value", cleavages);
  }
  
  //check if maximum length is with in range <= 255
  if(max_length > 255){
    carp(CARP_FATAL, "maximum length:%d over limit 255.", max_length);
    exit(1);
  }
  
  //determine isotopic mass option
  if(strcmp(isotopic_mass, "average")==0){
    mass_type = AVERAGE;
  }
  else if(strcmp(isotopic_mass, "mono")==0){
    mass_type = MONO;
  }
  else{
    carp(CARP_ERROR, "incorrect argument %s, using default value", isotopic_mass);
  }
   
  //determine redundancy option
  if(strcmp(redundancy, "redundant")==0){
    is_unique = FALSE;
  }
  else if(strcmp(redundancy, "unique")==0){
    is_unique = TRUE;
  }
  else{
    carp(CARP_ERROR, "incorrect argument %s, using default value", redundancy);
  }
  
  //determine sort type option
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
    carp(CARP_ERROR, "incorrect argument %s, using default value", sort);
  }
    
  //determine use index command
  if(strcmp(use_index, "F")==0){
    use_index_boolean = FALSE;
  }
  else if(strcmp(use_index, "T")==0){
    use_index_boolean = TRUE;
  }
  else{
    carp(CARP_ERROR, "incorrect argument %s, using default value", use_index);
  }

  //check if input file exist
  if(access(in_file, F_OK)){
    carp(CARP_FATAL, "The file \"%s\" does not exist (or is not readable, or is empty).", in_file);
    exit(1);
  }
 
  //allocate an empty iterator
  GENERATE_PEPTIDES_ITERATOR_T* gen_peptide_iterator = allocate_generate_peptides_iterator();
  
  //peptide constraint
  PEPTIDE_CONSTRAINT_T* constraint = 
    new_peptide_constraint(peptide_type, min_mass, max_mass, min_length, max_length, missed_cleavages, mass_type);
  
  //asign to iterator
  gen_peptide_iterator->constraint = constraint;


  /***********************
   * use index file
   **********************/
  if(use_index_boolean){
    if((sort_type != MASS && sort_type != NONE)){
      carp(CARP_ERROR, " when using index, cannot sort other than by mass");
      carp(CARP_ERROR, "failed to perform search");
      exit(1);
    }
    
    //create index and set to generate_peptides_iterator
    index = new_search_index(in_file, constraint, is_unique);
    gen_peptide_iterator->index = index;
    
    if(index != NULL){
      //only resrict peptide by mass and length, default iterator
      if(peptide_type == ANY_TRYPTIC){
        //create index peptide interator & set generate_peptides_iterator
        index_peptide_iterator = new_index_peptide_iterator(index);        
        gen_peptide_iterator->iterator = index_peptide_iterator;
        gen_peptide_iterator->has_next = &void_index_peptide_iterator_has_next;
        gen_peptide_iterator->next = &void_index_peptide_iterator_next;
        gen_peptide_iterator->free = &void_free_index_peptide_iterator;
      }
      //if need to select among peptides by peptide_type and etc.
      else{
        //create index_filtered_peptide_iterator  & set generate_peptides_iterator
        index_filtered_peptide_iterator = new_index_filtered_peptide_iterator(index);
        gen_peptide_iterator->iterator = index_filtered_peptide_iterator;
        gen_peptide_iterator->has_next = &void_index_filtered_peptide_iterator_has_next;
        gen_peptide_iterator->next = &void_index_filtered_peptide_iterator_next;
        gen_peptide_iterator->free = &void_free_index_filtered_peptide_iterator;
      }
      //set the correct free peptide method
      gen_peptide_iterator->free_peptide = &free_peptide_for_array;
    }
  }
  /*********************************************
   *read in from fasta file, don't use index file
   ************************************************/
  else{
    //create a new database & set generate_peptides_iterator
    database = new_database(in_file, FALSE);         //needs to change this....by given option
    gen_peptide_iterator->database = database;
    
    //no sort, redundant
    if(!is_unique && sort_type == NONE){ 
      
      //create peptide iterator  & set generate_peptides_iterator
      iterator = new_database_peptide_iterator(database, constraint);
      gen_peptide_iterator->iterator = iterator;
      gen_peptide_iterator->has_next = &void_database_peptide_iterator_has_next;
      gen_peptide_iterator->next = &void_database_peptide_iterator_next;
      gen_peptide_iterator->free = &void_free_database_peptide_iterator;
      
    }      
    //sort or check for unique
    else{
      //only sort, by default will be sorted by mass
      if(sort_type == NONE){
        //create peptide iterator
        sorted_iterator = 
          new_database_sorted_peptide_iterator(database, constraint, MASS, TRUE);       
      }
      //create peptide iterator
      else{
        sorted_iterator = 
          new_database_sorted_peptide_iterator(database, constraint, sort_type, is_unique);
      }
      
      // set generate_peptides_iterator 
      gen_peptide_iterator->iterator = sorted_iterator;
      gen_peptide_iterator->has_next = &void_database_sorted_peptide_iterator_has_next;
      gen_peptide_iterator->next = &void_database_sorted_peptide_iterator_next;
      gen_peptide_iterator->free = &void_free_database_sorted_peptide_iterator;
      
    }
    //set the correct free peptide method
    gen_peptide_iterator->free_peptide = &free_peptide;
  }
  return gen_peptide_iterator;
}

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
  //free the nested iterator
  generate_peptides_iterator->free(generate_peptides_iterator->iterator);

  //free database or index, if exist
  if(generate_peptides_iterator->database != NULL){
    free_database(generate_peptides_iterator->database);
  }
  if(generate_peptides_iterator->index != NULL){
    free_index(generate_peptides_iterator->index);
  }
  
  //free peptide constraint
  free_peptide_constraint(generate_peptides_iterator->constraint);

  free(generate_peptides_iterator);
}

/**
 * Always free peptides created by the generate_peptides_iterator through this method
 * Frees the allocated peptide with the correct free-method according to it's type
 */
void free_peptide_produced_by_iterator(
  GENERATE_PEPTIDES_ITERATOR_T* generate_peptides_iterator, ///< the iterator which the peptide was produced -in
  PEPTIDE_T* peptide ///< the peptide to free -in
  )
{
  //free peptide with the correct free peptide method
  generate_peptides_iterator->free_peptide(peptide);
}
