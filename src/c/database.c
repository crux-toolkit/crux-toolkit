/*****************************************************************************
 * \file database.c
 * $Revision: 1.24 $
 * \brief: Object for representing a database of protein sequences.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"
#include "crux-utils.h"
#include "peptide.h"
#include "protein.h"
#include "database.h"
#include "carp.h"
#include "objects.h"
#include "peptide_constraint.h"

#define MAX_PROTEINS 2500000 ///< The maximum number of proteins in a database.

/**
 * \struct database
 * \brief A database of protein sequences
 */
struct database{
  char*        filename;      ///< Original database filename.
  FILE*        file;          ///< Open filehandle for this database.
                                ///  A database has only one
                                ///  associated file.
  int num_proteins;             ///< Number of proteins in this database.
  BOOLEAN_T is_parsed;          ///< Has this database been parsed yet.
  PROTEIN_T* proteins[MAX_PROTEINS];   ///< Proteins in this database 
  unsigned long int size; ///< The size of the database in bytes (convenience)
};    

/**
 * \struct database_protein_iterator
 * \brief Object to iterate over the proteins within a database
 */
struct database_protein_iterator {
  DATABASE_T* database;  ///< The database whose proteins to iterate over 
  int cur_protein;      ///< The index of the current protein
};

/**
 * \struct database_peptide_iterator
 * \brief Object to iterate over the peptides within a database, in an
 * unspecified order.
 */
struct database_peptide_iterator {
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator; ///<The protein iterator. 
  PROTEIN_PEPTIDE_ITERATOR_T* 
    cur_protein_peptide_iterator; ///< The peptide iterator for the current protein.
  PEPTIDE_CONSTRAINT_T* peptide_constraint; ///< The constraints for the kind of peptide to iterate over.
  };

/**
 * \struct database_sorted_peptide_iterator
 * \brief Object to iterate over the peptides within a database, in an
 * specified sorted order.(mass, length, lexical)
 */
struct database_sorted_peptide_iterator {
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator; ///<The protein iterator. 
  PROTEIN_PEPTIDE_ITERATOR_T* 
    cur_protein_peptide_iterator; ///< The peptide iterator for the current protein.
  PEPTIDE_CONSTRAINT_T* peptide_constraint; ///< The constraints for the kind of peptide to iterate over.
  SORT_TYPE_T sort_type; ///< The sort type for this iterator (MASS, LENGTH, LEXICAL);
  PEPTIDE_WRAPPER_T* peptide_wrapper; ///< a linklist of peptide wrappers
};

/**
 * \struct peptide_wrapper
 * \brief Object to wrap the peptide to build a linklist
 */
struct peptide_wrapper{
  PEPTIDE_WRAPPER_T* next_wrapper; ///< the next peptide wrapper
  PEPTIDE_T* peptide;   ///< the core, the peptide
};    

/**
 * \returns An (empty) database object.
 */
DATABASE_T* allocate_database(void){
  DATABASE_T* database = (DATABASE_T*)mycalloc(1,sizeof(DATABASE_T));
  database->is_parsed = FALSE;
  return database;
}

/**
 * \returns A new database object.
 */
DATABASE_T* new_database(
  char*         filename ///< The file from which to parse the database. -in
  )
{
  DATABASE_T* database = allocate_database();
  set_database_filename(database, filename);
  return database;
}  

/**
 * Frees an allocated protein object.
 */
void free_database(
  DATABASE_T* database ///< An allocated database -in
  )
{
  int protein_idx = 0;
  
  free(database->filename);
  
  //only free proteins if been parsed and file has been opened
  if(database->is_parsed){
    fclose(database->file);
    //free each protein in the array
    for(; protein_idx < database->num_proteins; ++protein_idx){
      free_protein(database->proteins[protein_idx]);
    }
  }
  free(database);
}

/**
 * Prints a database object to file.
 */
void print_database(
  DATABASE_T* database,  ///< database to print -in
  FILE* file    ///< output file stream -out             
  )
{
 
  fprintf(file, "filename:%s\n", database->filename);
  fprintf(file, "is_parsed:");
  
  //has the database been parsed?
  if(database->is_parsed){
    fprintf(file, "TRUE\n");
    DATABASE_PROTEIN_ITERATOR_T* iterator = new_database_protein_iterator(database);
 
    while(database_protein_iterator_has_next(iterator)){
      print_protein(database_protein_iterator_next(iterator), stdout);
    }
    free_database_protein_iterator(iterator);
  }
  else{
    fprintf(file, "FALSE\n");
  }
}

/**
 * Parses a database from the file in the filename member variable
 * reads in all proteins in the fasta file and creates a protein object
 * and adds them to the database protein array
 * total proteins in fasta file must not exceed MAX_PROTEIN constant
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_database(
  DATABASE_T* database ///< An allocated database -in
  )
{
  unsigned long working_index;
  FILE* file = NULL;
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  PROTEIN_T* new_protein;
  int protein_idx = 0;

  //check if already parsed
  if(database->is_parsed){
    return TRUE;
  }
  
  //open file and 
  file = fopen(database->filename, "r");

  working_index = ftell(file);
  //check each line until reach '>' line
  while((line_length =  getline(&new_line, &buf_length, file)) != -1){
    if(new_line[0] == '>'){
      if(database->num_proteins == MAX_PROTEINS+1){
        fclose(file);
        free(new_line);
        carp(CARP_ERROR, "exceeds protein index array size");
        return FALSE;
      }
      //the new protein to be parsed
      new_protein = allocate_protein();

      //rewind to the begining of the protein to include ">" line
      fseek(file, working_index, SEEK_SET);
      
      //failed to parse the protein from fasta file
      //protein offset is set in the parse_protein_fasta_file method
      if(!parse_protein_fasta_file(new_protein ,file)){
        fclose(file);
        free_protein(new_protein);
        for(; protein_idx < database->num_proteins; ++protein_idx){
          free_protein(database->proteins[protein_idx]);
        }
        database->num_proteins = 0;
        carp(CARP_ERROR, "failed to parse fasta file");
        return FALSE;
      }

      database->proteins[database->num_proteins] = new_protein;

      /* insert code need for light/heavy functionality for proteins
      set_protein_offset(new_protein, working_index);
      */
      ++database->num_proteins;
      //set protein index
      set_protein_protein_idx(new_protein, database->num_proteins);
    }
    working_index = ftell(file);
  }
  //yes the database is paresed now..!!
  database->is_parsed = TRUE;
  free(new_line);
  database->file = file;
  return TRUE;
}


//FIXME needs to be implemented at some stage..if needed...
/**
 * \returns FALSE if database has not yet been parsed or if the nth protein
 * cannot be parsed.
 */
/*
BOOLEAN_T get_database_protein_at_idx(
    DATABASE_T* database, ///< A parsed database object -in
    int protein_idx,      ///< The index of the protein to retrieve -in
    PROTEIN_T** protein   ///< A pointer to a pointer to a PROTEIN object -out
    );
*/

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 */

/**
 *\returns the filename of the database
 * returns a heap allocated new copy of the filename
 * user must free the return filename
 */
char* get_database_filename(
  DATABASE_T* database ///< the query database -in 
  )
{
  return my_copy_string(database->filename);
}

/**
 * sets the filename of the database
 * protein->sequence must been initiailized
 */
void set_database_filename(
  DATABASE_T* database, ///< the database to set it's fields -out
  char* filename ///< the filename to add -in
  )
{
  free(database->filename);
  database->filename = my_copy_string(filename);
}

/**
 *\returns TRUE|FALSE whether the database has been parsed?
 */
BOOLEAN_T get_database_is_parsed(
  DATABASE_T* database ///< the query database -in 
  )
{
  return database->is_parsed;
}

/**
 *\returns the total number of proteins of the database
 */
int get_database_num_proteins(
  DATABASE_T* database ///< the query database -in 
  )
{
  return database->num_proteins;
}

/**
 *\returns the src FILE* of the database
 */
FILE* get_database_file(
  DATABASE_T* database ///< the query database -in 
  )
{
  return database->file;
}

/**
 * sets the src FILE* of the database
 */
void set_database_file(
  DATABASE_T* database, ///< the database to set it's fields -out
  FILE* file ///< the src file to add -in
  )
{
  database->file = file;
}



/***********************************************
 * Iterators
 ***********************************************/

/**
 * Instantiates a new database_protein_iterator from a database.
 * \returns a DATABASE_PROTEIN_ITERATOR_T object.
 */
DATABASE_PROTEIN_ITERATOR_T* new_database_protein_iterator(
  DATABASE_T* database ///< the database to create a protein iterator -in
  )
{
  //if database is parsed, if not do so..
  if(!database->is_parsed){
    //failed to parse database
    if(!parse_database(database)){
      carp(CARP_FATAL, "failed to parse database, cannot create iterator");
      exit(1);
    }
  }
  
  //create new protein iterator
  DATABASE_PROTEIN_ITERATOR_T* iterator = 
    (DATABASE_PROTEIN_ITERATOR_T*)mycalloc(1, sizeof(DATABASE_PROTEIN_ITERATOR_T));
  iterator->database = database;
  iterator->cur_protein = 0;
  
  return iterator;
}        


/**
 * Frees an allocated database_protein_iterator object.
 */
void free_database_protein_iterator(
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator ///< the iterator to free -in
  )
{
  free(database_protein_iterator);
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional proteins to iterate over, FALSE if not.
 */
BOOLEAN_T database_protein_iterator_has_next(
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator ///< the query iterator -in
  )
{
  if(database_protein_iterator->cur_protein <
     get_database_num_proteins(database_protein_iterator->database)){
    return TRUE;
  }
  
  return FALSE;
}

/**
 * \returns The next protein in the database.
 */
PROTEIN_T* database_protein_iterator_next(
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator  ///< the query iterator -in
  )
{
  ++database_protein_iterator->cur_protein;

  //print number of protein generated to STDERR for every 500 protein reached
  if(database_protein_iterator->cur_protein % 500 == 0){
    carp(CARP_INFO, "reached protein %d out of %d", 
         database_protein_iterator->cur_protein,
         database_protein_iterator->database->num_proteins);
  }

  return database_protein_iterator->database->proteins[database_protein_iterator->cur_protein-1];
}

/**
 * \returns the protein to the corresponding protein_idx in the database.
 */
PROTEIN_T* database_protein_iterator_protein_idx(
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator, ///< the iterator of interest -in
  int protein_idx ///< protein_idx to which protein to return -in
  )
{
  if(!database_protein_iterator_has_next(database_protein_iterator)){
    die("no proteins to return in protein iterator");
  }
  else if(protein_idx < 1 || //might have to be 1 
          protein_idx > database_protein_iterator->database->num_proteins){
    die("protein_idx index out of bounds");
  }
  return database_protein_iterator->database->proteins[protein_idx-1];
}

/***********************************************
 * database_peptide_Iterators
 ***********************************************/

/**
 * Instantiates a new database_peptide_iterator from a database.
 * \returns a DATABASE_PEPTIDE_ITERATOR_T object.
 */
DATABASE_PEPTIDE_ITERATOR_T* new_database_peptide_iterator(
  DATABASE_T* database, ///< the database of interest -in
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide_constraint to filter peptides -in
  )
{
  PROTEIN_T* next_protein;
  
  DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator =
    (DATABASE_PEPTIDE_ITERATOR_T*)mycalloc(1, sizeof(DATABASE_PEPTIDE_ITERATOR_T));

  //set a new protein iterator
  database_peptide_iterator->database_protein_iterator =
    new_database_protein_iterator(database);

  //set peptide constraint
  database_peptide_iterator->peptide_constraint = peptide_constraint;
  
  //check if there's any proteins to create peptides from
  if(database_protein_iterator_has_next(database_peptide_iterator->database_protein_iterator)){
    next_protein =
      database_protein_iterator_next(database_peptide_iterator->database_protein_iterator);

    //set new protein peptide iterator
    database_peptide_iterator->cur_protein_peptide_iterator =
      new_protein_peptide_iterator(next_protein, database_peptide_iterator->peptide_constraint);
 
    //if first protein does not contain a match peptide, reinitailize
    while(!protein_peptide_iterator_has_next(database_peptide_iterator->cur_protein_peptide_iterator)){
      //end of list of peptides for database_peptide_iterator
      if(!database_protein_iterator_has_next(database_peptide_iterator->database_protein_iterator)){
        break;
      }
      else{ //create new protein_peptide_iterator for next protein
        free_protein_peptide_iterator(database_peptide_iterator->cur_protein_peptide_iterator);
        database_peptide_iterator->cur_protein_peptide_iterator =
          new_protein_peptide_iterator(
            database_protein_iterator_next(database_peptide_iterator->database_protein_iterator), 
            database_peptide_iterator->peptide_constraint);
      }
    }
  }
  else{ //no proteins to create peptides from
    carp(CARP_FATAL, "failed to create a database_peptide_iterator, no proteins in database");
    free_database_protein_iterator(database_peptide_iterator->database_protein_iterator);
    free(database_peptide_iterator);
    exit(1);
  }
  return database_peptide_iterator;
}


/**
 * Frees an allocated database_peptide_iterator object.
 */
void free_database_peptide_iterator(
  DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator ///< the iterator to free -in
  )
{
  free_protein_peptide_iterator(database_peptide_iterator->cur_protein_peptide_iterator);
  free_database_protein_iterator(database_peptide_iterator->database_protein_iterator);
  free(database_peptide_iterator);
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T database_peptide_iterator_has_next(
  DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator  ///< the iterator of interest -in
  )
{
  if(protein_peptide_iterator_has_next(database_peptide_iterator->cur_protein_peptide_iterator)){ 
    return TRUE;
  }
  return FALSE;
}

/**
 * \returns The next peptide in the database.
 */
PEPTIDE_T* database_peptide_iterator_next(
  DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator ///< the iterator of interest -in
  )
{
  PEPTIDE_T* next_peptide =
    protein_peptide_iterator_next(database_peptide_iterator->cur_protein_peptide_iterator);
  
  //reset database_peptide_iterator if needed
  while(!protein_peptide_iterator_has_next(database_peptide_iterator->cur_protein_peptide_iterator)){
    //end of list of peptides for database_peptide_iterator
    if(!database_protein_iterator_has_next(database_peptide_iterator->database_protein_iterator)){
      break;
    }
    else{ //create new protein_peptide_iterator for next protein
      free_protein_peptide_iterator(database_peptide_iterator->cur_protein_peptide_iterator);
      database_peptide_iterator->cur_protein_peptide_iterator =
        new_protein_peptide_iterator(
          database_protein_iterator_next(database_peptide_iterator->database_protein_iterator), 
          database_peptide_iterator->peptide_constraint);
    }
  }
  return next_peptide;
}

/***********************************
 * database sorted peptide iterator
 ***********************************/

/**
 *wrap the peptide up in the a peptide_wrapper object
 *\returns a peptide_wrapper object that contains the peptide
 */
PEPTIDE_WRAPPER_T* wrap_peptide(
  PEPTIDE_T* peptide ///< peptide to be wrapped -in
  )
{
  PEPTIDE_WRAPPER_T* new_wrapper = (PEPTIDE_WRAPPER_T*)mycalloc(1, sizeof(PEPTIDE_WRAPPER_T));
  new_wrapper->peptide = peptide;
  return new_wrapper;    
}


/**
 *Initialize the peptide iterator
 *\returns FALSE if there are no peptides to return, TRUE if peptide iterator has at least one peptide
 */
BOOLEAN_T initialize_peptide_iterator(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator, ///< initialize iterator -in
  DATABASE_T* database, ///< the database of interest -in
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide_constraint to filter peptides -in
)
{
  PROTEIN_T* next_protein = NULL;

  //set a new protein iterator
   database_sorted_peptide_iterator->database_protein_iterator = new_database_protein_iterator(database);

  //set peptide constraint
  database_sorted_peptide_iterator->peptide_constraint = peptide_constraint;

  //are there any proteins to create peptides?
  if(database_protein_iterator_has_next(database_sorted_peptide_iterator->database_protein_iterator)){
    next_protein =
      database_protein_iterator_next(database_sorted_peptide_iterator->database_protein_iterator);

    //set new protein peptide iterator
    database_sorted_peptide_iterator->cur_protein_peptide_iterator =
      new_protein_peptide_iterator(next_protein, database_sorted_peptide_iterator->peptide_constraint);
 
    //if first protein does not contain a match peptide, reinitailize
    while(!protein_peptide_iterator_has_next(database_sorted_peptide_iterator->cur_protein_peptide_iterator)){
      //end of list of peptides for database_peptide_iterator
      if(!database_protein_iterator_has_next(database_sorted_peptide_iterator->database_protein_iterator)){
        break;
      }
      else{ //create new protein_peptide_iterator for next protein
        free_protein_peptide_iterator(database_sorted_peptide_iterator->cur_protein_peptide_iterator);
        database_sorted_peptide_iterator->cur_protein_peptide_iterator =
          new_protein_peptide_iterator(
            database_protein_iterator_next(database_sorted_peptide_iterator->database_protein_iterator), 
            peptide_constraint);
      }
    }
    if(protein_peptide_iterator_has_next(database_sorted_peptide_iterator->cur_protein_peptide_iterator)){
      return TRUE;
    }
  }
  
  return FALSE;
}

/**
 * compares two peptides with the give sort type (length, mass, lexical)
 * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
 */
int compareTo(
  PEPTIDE_T* peptide_one, ///< peptide to compare one -in
  PEPTIDE_T* peptide_two, ///< peptide to compare two -in
  SORT_TYPE_T sort_type  ///< sort type(LENGTH, MASS, LEXICAL) -in
  )
{
  //length order
  if(sort_type == LENGTH){
    if(get_peptide_length(peptide_one) >
       get_peptide_length(peptide_two)){
      return 1;
    }
    else if(get_peptide_length(peptide_one) ==
            get_peptide_length(peptide_two)){
      return 0;
    }
    else{
      return -1;
    }
  }
  //mass order
  else if(sort_type == MASS){
    return compare_float(get_peptide_peptide_mass(peptide_one), 
                         get_peptide_peptide_mass(peptide_two));
  }
  //lexicographic order
  else if(sort_type == LEXICAL){
    char* peptide_one_sequence = get_peptide_sequence_pointer(peptide_one);
    char* peptide_two_sequence = get_peptide_sequence_pointer(peptide_two);
    int peptide_one_length = get_peptide_length(peptide_one);
    int peptide_two_length = get_peptide_length(peptide_two);
    int current_idx = 0;
    
    //check if all alphabetically identical
    while(current_idx < peptide_one_length &&
          current_idx < peptide_two_length){
      if(peptide_one_sequence[current_idx] > peptide_two_sequence[current_idx]){
        return 1;
      }
      else if(peptide_one_sequence[current_idx] < peptide_two_sequence[current_idx]){
        return -1;
      }
      ++current_idx;
    }
    
    //alphabetically identical, check if same length
    if(peptide_one_length > peptide_two_length){
      return 1;
    }
    else if(peptide_one_length < peptide_two_length){
      return -1;
    }
    else{
      return 0;
    }
  }
  
  die("ERROR: no matching sort_type");
  //quiet compiler
  return 0;
}

/**
 *check if any duplicate peptides exist to the first peptide on the wrapper
 *then, merge duplicate peptides to the first peptide 
 *assumes that the list is sorted by mass
 *\returns the wrapper list which first peptide is unique in the list
 */
PEPTIDE_WRAPPER_T* merge_duplicates_same_list(
  PEPTIDE_WRAPPER_T* wrapper_list  ///< peptide wrapper list examine -mod
  )
{
  PEPTIDE_WRAPPER_T* current = wrapper_list;
  //for all peptides that have equal mass
  while(current->next_wrapper != NULL &&
        compareTo(wrapper_list->peptide, current->next_wrapper->peptide, MASS)==0){  
    //check identical peptide
    if(compareTo(wrapper_list->peptide, current->next_wrapper->peptide, LEXICAL)==0){
      PEPTIDE_WRAPPER_T* wrapper_to_delete = current->next_wrapper;
      //merge peptides
      merge_peptides(wrapper_list->peptide, wrapper_to_delete->peptide);
      current->next_wrapper = current->next_wrapper->next_wrapper;
      free(wrapper_to_delete);
    }
    else{
      current = current->next_wrapper;
    }
  }
  return wrapper_list;
}

/**
 *check if any duplicate peptides exist in the wrapper_list_second
 *compared to the first peptide of the wrapper_list_first
 *then, remove the peptide from seond list and merge 
 *duplicate peptides to the first peptide in the first list
 *assumes that the list is sorted by mass
 *\returns the modified second wrapper list
 */
PEPTIDE_WRAPPER_T* merge_duplicates_different_list(
  PEPTIDE_WRAPPER_T* wrapper_list_first,  ///< fist peptide wrapper list examine -mod
  PEPTIDE_WRAPPER_T* wrapper_list_second  ///<second peptide wrapper list examine -mod
  )
{
  PEPTIDE_WRAPPER_T* current = wrapper_list_second;
  PEPTIDE_WRAPPER_T* previous = wrapper_list_second;
  
  while(current != NULL &&
        compareTo(wrapper_list_first->peptide, current->peptide, MASS) == 0){
    
    if(compareTo(wrapper_list_first->peptide, current->peptide, LEXICAL) == 0){
      PEPTIDE_WRAPPER_T* wrapper_to_delete = current;
      //merge peptides
      merge_peptides(wrapper_list_first->peptide, wrapper_to_delete->peptide);
      //for first wrapper in the list 
      if(previous == current){
        current = current->next_wrapper;
        previous = current;
        wrapper_list_second = current;
      }
      else{
        current = current->next_wrapper;
        previous->next_wrapper = current;
      }
      free(wrapper_to_delete);
    }
    else{
      if(previous == current){
        current = current->next_wrapper;
      }
      else{
        previous = current;
        current = current->next_wrapper;
      }
    }
  }
  return wrapper_list_second;
}

/**
 * merge wrapper list one and list two by sort type
 * assumes that each list one and two are sorted by the same sort type
 *\returns a sorted wrapper list result of merging list one and two
 */   
PEPTIDE_WRAPPER_T* merge(
  PEPTIDE_WRAPPER_T* wrapper_one, ///< fist peptide wrapper list -in
  PEPTIDE_WRAPPER_T* wrapper_two, ///< second peptide wrapper list -in
  SORT_TYPE_T sort_type, ///< the sort type of the new merged list -in
  BOOLEAN_T unique ///< do you merge two lists into a uniqe list? -in
  )
{
  PEPTIDE_WRAPPER_T* wrapper_final = NULL;
  PEPTIDE_WRAPPER_T* wrapper_current = wrapper_final;
  int compared;

  //loop around until there are no more wrappers to merge
 LOOP:

  if(wrapper_one != NULL || wrapper_two != NULL){
    if(wrapper_one == NULL){
      //there are wrappers on the current list
      if(wrapper_current != NULL){
        wrapper_current->next_wrapper = wrapper_two;
      }
      //there are no wrappers on the current list
      else{
        wrapper_final = wrapper_two;
      }
      wrapper_two = NULL;
    }
    else if(wrapper_two == NULL){
      //there are wrappers on the current list
      if(wrapper_current != NULL){
        wrapper_current->next_wrapper = wrapper_one;        
      }
      //there are no wrappers on the current list
      else{
        wrapper_final = wrapper_one;
      }
      wrapper_one = NULL;
    }
    
    else if((compared = 
             compareTo(wrapper_one->peptide, wrapper_two->peptide, sort_type)) == 1){
      //there are wrappers on the current list
      if(wrapper_current != NULL){
        wrapper_current->next_wrapper = wrapper_two;        
        wrapper_current = wrapper_current->next_wrapper;
      }
      //there are no wrappers on the current list
      else{
        wrapper_current = wrapper_two;
        wrapper_final = wrapper_current;
      }
      wrapper_two = wrapper_two->next_wrapper;
    }
    else{
      int compared_mass = -1;

      //if duplicate peptide, merge into one peptide only if unique is TRUE
      if(compared == 0 && unique){
        //only check if same peptide if mass is identical
        if(sort_type == MASS ||
           (compared_mass = compareTo(wrapper_one->peptide, wrapper_two->peptide, MASS))==0){
             //must be identical peptide, since mass & lexicographically same
          if(sort_type == LEXICAL ||
             compareTo(wrapper_one->peptide, wrapper_two->peptide, LEXICAL)==0){
            
            PEPTIDE_WRAPPER_T* wrapper_to_delete = wrapper_two;
            //merge peptides
            merge_peptides(wrapper_one->peptide, wrapper_two->peptide);
            wrapper_two = wrapper_two->next_wrapper;
            free(wrapper_to_delete);
          }
        }
        //merge all other instances of the same peptide in the list before adding to master list
        if(sort_type == MASS || sort_type == LENGTH){
          wrapper_one = merge_duplicates_same_list(wrapper_one);
          wrapper_two = merge_duplicates_different_list(wrapper_one, wrapper_two);
        }
      }

      //when sorting by length and unique, sort also by mass to speed up the unique check
      //thus, when identical length add the peptide with smaller mass to the list first
      if(unique && sort_type == LENGTH && compared == 0 && compared_mass == 1){
        //there are wrappers on the current list
        if(wrapper_current != NULL){
          wrapper_current->next_wrapper = wrapper_two;        
          wrapper_current = wrapper_current->next_wrapper;
        }
        //there are no wrappers on the current list
        else{
          wrapper_current = wrapper_two;
          wrapper_final = wrapper_current;
        }
        wrapper_two = wrapper_two->next_wrapper;
      }
      else{
        //add wrapper to the merged list
        //there are wrappers on the current list
        if(wrapper_current != NULL){
          wrapper_current->next_wrapper = wrapper_one; 
          wrapper_current = wrapper_current->next_wrapper;
        }
        //there are no wrappers on the current list
        else{
          wrapper_current = wrapper_one;
          wrapper_final = wrapper_current;
        }
        wrapper_one = wrapper_one->next_wrapper;
      }
    }
    goto LOOP;
  }
  return wrapper_final;
}


/**
 * spilt a wrapper list into two equal or almost equal size lists
 * \returns array of two pointers of each split array
 * user must free the array, heap allocated
 */

PEPTIDE_WRAPPER_T** split(
  PEPTIDE_WRAPPER_T* wrapper_src ///< wrapper list to split -in
  )
{
  PEPTIDE_WRAPPER_T** split_wrapper = 
    (PEPTIDE_WRAPPER_T**)mycalloc(2, sizeof(PEPTIDE_WRAPPER_T*));

  PEPTIDE_WRAPPER_T* split1 = NULL;
  PEPTIDE_WRAPPER_T* split1_current = split1;
  PEPTIDE_WRAPPER_T* split2 = NULL;
  PEPTIDE_WRAPPER_T* split2_current = split2;
  BOOLEAN_T one = TRUE;

  //split wrapper until no more wrapper left in wrapper_src
  while(wrapper_src != NULL){
    //do i add it to the first list?
    if(one){
      //first time adding wrapper
      if(split1_current == NULL){
        split1_current = wrapper_src;
        split1 = split1_current;       
      }
      //add one wrapper to the end of the first list 
      else{
        split1_current->next_wrapper = wrapper_src;
        split1_current = split1_current->next_wrapper;
      }
      one = FALSE;
    }
    //add to second list
    else{
      //first time adding wrapper
      if(split2_current == NULL){
        split2_current = wrapper_src;
        split2 = split2_current;
      }
      //add one wrapper to the end of the second list 
      else{
        split2_current->next_wrapper = wrapper_src;
        split2_current = split2_current->next_wrapper;
      }
      one = TRUE;
    }
    wrapper_src = wrapper_src->next_wrapper;
  }
  
  //add list one and two to thesplit_wrapper to return
  if(split1 != NULL){
    split1_current->next_wrapper = NULL;
    split_wrapper[0] = split1;
  }

  if(split2 != NULL){
    split2_current->next_wrapper = NULL;
    split_wrapper[1] = split2;
  }

  return split_wrapper;
}

/**
 * merge sort wrapper list by the sort type
 *\returns a sorted wrapper list
 */   
PEPTIDE_WRAPPER_T* merge_sort(
  PEPTIDE_WRAPPER_T* wrapper_list, ///< the list of wrappers to sort -in
  SORT_TYPE_T sort_type, ///<the sort type (length, mass, lexicographical) -in
  BOOLEAN_T unique ///< return a list of unique peptides? -in
  )
{
  PEPTIDE_WRAPPER_T* final_sorted_list = NULL;
  PEPTIDE_WRAPPER_T* first_half = NULL;
  PEPTIDE_WRAPPER_T* second_half = NULL;
  
  //split wrapper list into two equal lists
  PEPTIDE_WRAPPER_T** wrapper_split = 
    split(wrapper_list);

  first_half = wrapper_split[0];
  second_half = wrapper_split[1];
  free(wrapper_split);


  //recursivelly sort each half
  if(first_half != NULL && second_half != NULL){
    first_half = merge_sort(first_half, sort_type, unique);
    second_half = merge_sort(second_half, sort_type, unique);
    final_sorted_list = merge(first_half, second_half, sort_type, unique);
  }
  else if(first_half != NULL || second_half != NULL){
    final_sorted_list = merge(first_half, second_half, sort_type, unique);
  }
  return final_sorted_list;
}

/**
 * Frees an allocated database_sorted_peptide_iterator object.
 * does not free the peptide it contains
 */
void free_peptide_wrapper(
  PEPTIDE_WRAPPER_T* peptide_wrapper ///< the wrapper to free -in
  )
{
  free(peptide_wrapper);
}

/**
 * Frees an allocated database_sorted_peptide_iterator object.
 * This frees the peptide the wrapper contains as well
 */
void free_peptide_wrapper_all(
  PEPTIDE_WRAPPER_T* peptide_wrapper ///< the wrapper to free -in
  )
{
  free_peptide(peptide_wrapper->peptide);
  free(peptide_wrapper);
}

unsigned long total_number_peptide = 0;
/**
 * Instantiates a new database_sorted_peptide_iterator from a database.
 * \returns a DATABASE_SORTED_PEPTIDE_ITERATOR_T object.
 */
DATABASE_SORTED_PEPTIDE_ITERATOR_T* new_database_sorted_peptide_iterator(
  DATABASE_T* database, ///< the database of interest -in
  PEPTIDE_CONSTRAINT_T* peptide_constraint, ///< the peptide_constraint to filter peptides -in
  SORT_TYPE_T sort_type, ///< the sort type for this iterator -in
  BOOLEAN_T unique ///< only return unique peptides? -in
  )
{
  PEPTIDE_WRAPPER_T* master_list_wrapper = NULL;
  PEPTIDE_WRAPPER_T* list_wrapper = NULL;
  PEPTIDE_WRAPPER_T* current_wrapper = NULL;
  BOOLEAN_T start = TRUE;

  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator =
    (DATABASE_SORTED_PEPTIDE_ITERATOR_T*)mycalloc(1, sizeof(DATABASE_SORTED_PEPTIDE_ITERATOR_T));
  
  //The protein iterator 
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator = NULL;
  
 /// The peptide iterator for the current protein
  PROTEIN_PEPTIDE_ITERATOR_T* cur_protein_peptide_iterator = NULL;

  //initialize, return true if there's at least one peptide to return
  if(initialize_peptide_iterator(
    database_sorted_peptide_iterator, 
    database,
    peptide_constraint)){
    
    //strip all sub iterators from database, just to avoid mutiple references
    database_protein_iterator = database_sorted_peptide_iterator->database_protein_iterator;
    cur_protein_peptide_iterator = database_sorted_peptide_iterator->cur_protein_peptide_iterator;
    database_sorted_peptide_iterator->database_protein_iterator = NULL;
    database_sorted_peptide_iterator->cur_protein_peptide_iterator = NULL;

    do{
      //move to next protein if no peptides from current protein
      if(!protein_peptide_iterator_has_next(cur_protein_peptide_iterator)){
        //if there's more proteins in database...
        if(database_protein_iterator_has_next(database_protein_iterator)){
          //create new protein_peptide_iterator for next protein
          free_protein_peptide_iterator(cur_protein_peptide_iterator);
          cur_protein_peptide_iterator =
            new_protein_peptide_iterator(
               database_protein_iterator_next(database_protein_iterator), peptide_constraint);
        }
      }
      //iterate over all peptides in a protein
      while(protein_peptide_iterator_has_next(cur_protein_peptide_iterator)){
         //debug purpuse
        ++total_number_peptide;
        if(total_number_peptide % 1000000 == 0){
           carp(DEBUG, "number of peptides(not unique): %u", total_number_peptide); 
        }

        if(start){
          start = FALSE;
          current_wrapper =
            wrap_peptide(protein_peptide_iterator_next(cur_protein_peptide_iterator));
          list_wrapper = current_wrapper;
        }
        else{
          //wrap the next protein
          current_wrapper->next_wrapper =
            wrap_peptide(protein_peptide_iterator_next(cur_protein_peptide_iterator));
          current_wrapper = current_wrapper->next_wrapper;
        }
        
        //end of peptides from one parent protein merge list into master list
        if(!protein_peptide_iterator_has_next(cur_protein_peptide_iterator)){
          //add all peptides to the complied master list
          if(master_list_wrapper == NULL){
            master_list_wrapper = list_wrapper;
          }
          else{
            current_wrapper->next_wrapper = master_list_wrapper;
            master_list_wrapper =  list_wrapper;
            current_wrapper = NULL;
          }
          list_wrapper = NULL;      
          start = TRUE;
        }
      }
    } //iterate over all proteins in database
    while(database_protein_iterator_has_next(database_protein_iterator));
    
    carp(DEBUG, "total number of peptides(not unique): %u", total_number_peptide); 

    //sort the master list using merge sort
    master_list_wrapper = merge_sort(master_list_wrapper, sort_type, unique);

    database_sorted_peptide_iterator->peptide_wrapper = master_list_wrapper;
  }
  else{ //no peptides to be created
    free_database_protein_iterator(database_sorted_peptide_iterator->database_protein_iterator);
    free_protein_peptide_iterator(database_sorted_peptide_iterator->cur_protein_peptide_iterator);
    database_sorted_peptide_iterator->peptide_wrapper = NULL;
    return database_sorted_peptide_iterator;
  }

  free_database_protein_iterator(database_protein_iterator);
  free_protein_peptide_iterator(cur_protein_peptide_iterator);
  return database_sorted_peptide_iterator;
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T database_sorted_peptide_iterator_has_next(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator ///< the iterator of interest -in
  )
{
  if(database_sorted_peptide_iterator->peptide_wrapper != NULL){
    return TRUE;
  }
  return FALSE;
}

/**
 * returns each peptide in sorted order
 * \returns The next peptide in the database.
 */
PEPTIDE_T* database_sorted_peptide_iterator_next(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator ///< the iterator of interest -in
  )
{
  //strip the peptide out of the peptide wrapper
  PEPTIDE_T* next_peptide = database_sorted_peptide_iterator->peptide_wrapper->peptide;
  //free the empty peptide wrapper
  PEPTIDE_WRAPPER_T* old_wrapper =  database_sorted_peptide_iterator->peptide_wrapper;
  database_sorted_peptide_iterator->peptide_wrapper = 
    database_sorted_peptide_iterator->peptide_wrapper->next_wrapper;
  free_peptide_wrapper(old_wrapper);
  return next_peptide;
}

/**
 * Frees an allocated database_sorted_peptide_iterator object.
 */
void free_database_sorted_peptide_iterator(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator ///< the iterator to free -in
  )
{
  PEPTIDE_WRAPPER_T* old_wrapper = NULL;
  //free all peptide wrappers the iterator contains
  while(database_sorted_peptide_iterator->peptide_wrapper != NULL){
    old_wrapper = database_sorted_peptide_iterator->peptide_wrapper;
    database_sorted_peptide_iterator->peptide_wrapper = 
      database_sorted_peptide_iterator->peptide_wrapper->next_wrapper;
    free_peptide_wrapper_all(old_wrapper);
  }
  free(database_sorted_peptide_iterator);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

