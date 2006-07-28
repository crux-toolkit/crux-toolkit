/*****************************************************************************
 * \file database.c
 * $Revision: 1.12 $
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

#define MAX_PROTEINS 30000 ///< The maximum number of proteins in a database.

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
  PROTEIN_T* proteins[MAX_PROTEINS];   ///< Proteins in this database (not yet needed.)
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
 * specified sorted order.
 */
struct database_sorted_peptide_iterator {
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator; ///<The protein iterator. 
  PROTEIN_PEPTIDE_ITERATOR_T* 
    cur_protein_peptide_iterator; ///< The peptide iterator for the current protein.
  PEPTIDE_CONSTRAINT_T* peptide_constraint; ///< The constraints for the kind of peptide to iterate over.
  SORT_TYPE_T sort_type; ///< The sort type for this iterator (MASS, LENGTH);
  PEPTIDE_WRAPPER_T* peptide_wrapper; ///< a linklist of peptide wrappers
};

/**
 * \struct peptide_wrapper
 * \brief A database of protein sequences
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
  char*         filename ///< The file from which to parse the database.
  )
{
  DATABASE_T* database = allocate_database();
  set_database_filename(database, filename);
  return database;
}  

//FIXME think about what you're going to do for FILE* file
/**
 * Frees an allocated protein object.
 */
void free_database(
  DATABASE_T* database ///< An allocated database -in
  )
{
  int protein_idx = 0;

  free(database->filename);
  fclose(database->file);
  
  for(; protein_idx < database->num_proteins; ++protein_idx){
    free_protein(database->proteins[protein_idx]);
  }
  free(database);
}

/**
 * Prints a database object to file.
 */
void print_database(
  DATABASE_T* database, 
  FILE* file
  )
{
 
  fprintf(file, "filename:%s\n", database->filename);
  fprintf(file, "is_parsed:");
  

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

// START HERE

/*
 * Scans the database for start positions of protein sequences (using the
 * '>' character) and stores the locations of the starts in the database 
 * member variable starts. Set the is_parsed member variable to true.
 * IF false to parse protein, then frees all existing proteins, closes FILE* and resets num_protein to 0;
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
 * Additional get and set methods
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

/*
 * Scans the database for start positions of protein sequences (using the
 * '>' character) and stores the locations of the starts in the database 
 * member variable starts. Set the is_parsed member variable to true.
 */
BOOLEAN_T parse_database(
  DATABASE_T* database ///< An allocated database
  );

/***********************************
 * database sorted peptide iterator
 ***********************************/

/**
 * \struct database_sorted_peptide_iterator
 * \brief Object to iterate over the peptides within a database, in an
 * specified sorted order.
 */
struct database_sorted_peptide_iterator {
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator; ///<The protein iterator. 
  PROTEIN_PEPTIDE_ITERATOR_T* 
    cur_protein_peptide_iterator; ///< The peptide iterator for the current protein.
  PEPTIDE_CONSTRAINT_T* peptide_constraint; ///< The constraints for the kind of peptide to iterate over.
  SORT_TYPE_T sort_type; ///< The sort type for this iterator (MASS, LENGTH);
  PEPTIDE_WRAPPER_T* peptide_wrapper; ///< a linklist of peptide wrappers
};

/**
 * \struct peptide_wrapper
 * \brief A database of protein sequences
 */
struct peptide_wrapper{
  PEPTIDE_WRAPPER_T* next_wrapper; ///< the next peptide wrapper
  PEPTIDE_T* peptide;   ///< the core, the peptide
};    

/**
 * Instantiates a new database_sorted_peptide_iterator from a database.
 * \returns a DATABASE_SORTED_PEPTIDE_ITERATOR_T object.
 */
DATABASE_SORTED_PEPTIDE_ITERATOR_T* new_database_sorted_peptide_iterator(
  DATABASE_T* database, ///< the database of interest -in
  PEPTIDE_CONSTRAINT_T* peptide_constraint, ///< the peptide_constraint to filter peptides -in
  SORT_TYPE_T sort_type ///< the sort type for this iterator
  )
{
  PROTEIN_T* next_protein;
  PEPTIDE_T* next_peptide;
  PEPTIDE_WRAPPER_T* master_list_wrapper = NULL;
  PEPTIDE_WRAPPER_T* list_wrapper = NULL;
  PEPTIDE_WRAPPER_T* current_wrapper = NULL;
  PEPTIDE_WRAPPER_T* new_wrapper = NULL;
  BOOLEAN_T start = TRUE;

  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator =
    (DATABASE_SORTED_PEPTIDE_ITERATOR_T*)mycalloc(1, sizeof(DATABASE_SORTED_PEPTIDE_ITERATOR_T));

  //initialize
  if(initialize_peptide_iterator(database_sorted_peptide_iterator, database)){
    do{
      //skip first time
      if(!protein_peptide_iterator_has_next(database_peptide_iterator->cur_protein_peptide_iterator)){
        //if there's more proteins in database...
        if(database_protein_iterator_has_next(database_sorted_peptide_iterator->database_protein_iterator)){
          //create new protein_peptide_iterator for next protein
          free_protein_peptide_iterator(database_sorted_peptide_iterator->cur_protein_peptide_iterator);
          database_sorted_peptide_iterator->cur_protein_peptide_iterator =
            new_protein_peptide_iterator(
               database_protein_iterator_next(database_sorted_peptide_iterator->database_protein_iterator), 
               database_sorted_peptide_iterator->peptide_constraint);
        }
      }

      
      while(protein_peptide_iterator_has_next(database_peptide_iterator->cur_protein_peptide_iterator)){
        if(start){
          start = FALSE;
          current_wrapper =
            wrap_peptide(protein_peptide_iterator_next(database_peptide_iterator->cur_protein_peptide_iterator));
          list_wrapper = current_wrapper;
        }
        else{
          //wrap the next protein
          new_wrapper =
            wrap_peptide(protein_peptide_iterator_next(database_peptide_iterator->cur_protein_peptide_iterator));
          current_wrapper->next_wrapper = new_wrapper;
          current_wrapper = new_wrapper;
        }
        
        //end of peptides from one protein
        if(!protein_peptide_iterator_has_next(database_sorted_peptide_iterator->cur_protein_peptide_iterator)){
          if(master_list_wrapper == NULL){
            master_list_wrapper = list_wrapper;
          }
          else{
            master_list_wrapper = merge_sort(master_list_wrapper, list_wrapper, sort_type);
          }
          list_wrapper = NULL;      
          
        }
      }
    }
    while(database_protein_iterator_has_next(database_sorted_peptide_iterator->database_protein_iterator));

    database_sorted_peptide_iterator->peptide_wrapper = master_list_wrapper;
  }
  else{ //no proteins to create peptides from or peptides.
    carp(CARP_FATAL, "failed to create a database_peptide_iterator, no proteins in database");
    free_database_protein_iterator(database_sorted_peptide_iterator->database_protein_iterator);
    free(database_sorted_peptide_iterator);
    exit(1);
  }
  return database_sorted_peptide_iterator;
}

//merge sort
PEPTIDE_WRAPPER_T* merge_sort(
  PEPTIDE_WRAPPER_T* wrapper_one, ///<
  PEPTIDE_WRAPPER_T* wrapper_two, ///<
  SORT_TYPE_T* sort_type ///<
  )
{

  return wrapper_one;
}

//wrap the peptide up
PEPTIDE_WRAPPER_T* wrap_peptide(
  PEPTIDE_T* peptide
  )
{
  PEPTIDE_WRAPPER_T* new_wrapper = (PEPTIDE_WRAPPER_T*)mycalloc(1, sizeof(PEPTIDE_WRAPPER_T));
  new_wrapper->peptide = peptide;
  return new_wrapper;    
}



//check if there are any peptides at all to return????? i think so...
BOOLEAN_T initialize_peptide_iterator(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator,
  DATABASE_T* database, ///< the database of interest -in
)
{
  //set a new protein iterator
  database_sorted_peptide_iterator->database_protein_iterator =
    new_database_protein_iterator(database);

  //set peptide constraint
  database_sorted_peptide_iterator->peptide_constraint = peptide_constraint;

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
            database_sorted_peptide_iterator->peptide_constraint);
      }
    }
    if(protein_peptide_iterator_has_next(database_sorted_peptide_iterator->cur_protein_peptide_iterator)){
      return TRUE;
    }
  }
  
  return FALSE;
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T database_sorted_peptide_iterator_has_next(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_peptide_iterator ///< the iterator of interest -in
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
PEPTIDE_T* database_sorted_peptide_iterator_next(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_peptide_iterator ///< the iterator of interest -in
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

/**
 * Frees an allocated database_sorted_peptide_iterator object.
 */
void free_database_sorted_peptide_iterator(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_peptide_iterator ///< the iterator to free -in
  );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

