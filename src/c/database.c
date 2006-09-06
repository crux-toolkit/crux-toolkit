/*****************************************************************************
 * \file database.c
 * $Revision: 1.26 $
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
#include "sorter.h"

#define MAX_PROTEINS 3300000 ///< The maximum number of proteins in a database.

/**
 * \struct database
 * \brief A database of protein sequences
 */
struct database{
  char*        filename;      ///< Original database filename.
  FILE*        file;          ///< Open filehandle for this database.
                                ///  A database has only one
                                ///  associated file.
  unsigned int num_proteins;             ///< Number of proteins in this database.
  BOOLEAN_T is_parsed;          ///< Has this database been parsed yet.
  PROTEIN_T* proteins[MAX_PROTEINS];   ///< Proteins in this database 
  unsigned long int size; ///< The size of the database in bytes (convenience)
  BOOLEAN_T use_light_protein; //should I use the light/heavy protein option
};    

/**
 * \struct database_protein_iterator
 * \brief Object to iterate over the proteins within a database
 */
struct database_protein_iterator {
  DATABASE_T* database;  ///< The database whose proteins to iterate over 
  unsigned int cur_protein;      ///< The index of the current protein
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
  PROTEIN_T* prior_protein; ///< the protein that was used before the current working protein
  BOOLEAN_T first_passed; ///< is it ok to convert prior_protein to light?
};

/**
 * \struct database_sorted_peptide_iterator
 * \brief Object to iterate over the peptides within a database, in an
 * specified sorted order.(mass, length, lexical)
 */
struct database_sorted_peptide_iterator {
  SORTED_PEPTIDE_ITERATOR_T* sorted_peptide_iterator;
  //DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator; ///<The protein iterator. 
  //PROTEIN_PEPTIDE_ITERATOR_T* 
  //  cur_protein_peptide_iterator; ///< The peptide iterator for the current protein.
  //PEPTIDE_CONSTRAINT_T* peptide_constraint; ///< The constraints for the kind of peptide to iterate over.
  //SORT_TYPE_T sort_type; ///< The sort type for this iterator (MASS, LENGTH, LEXICAL);
  //PEPTIDE_WRAPPER_T* peptide_wrapper; ///< a linklist of peptide wrappers
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
  char*         filename, ///< The file from which to parse the database. -in
  BOOLEAN_T use_light_protein //should I use the light/heavy protein option
  )
{
  DATABASE_T* database = allocate_database();
  set_database_filename(database, filename);
  database->use_light_protein = use_light_protein;
  return database;
}  

/**
 * Frees an allocated protein object.
 */
void free_database(
  DATABASE_T* database ///< An allocated database -in
  )
{
  unsigned int protein_idx = 0;
  
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
  PROTEIN_T* protein = NULL;

  fprintf(file, "filename:%s\n", database->filename);
  fprintf(file, "is_parsed:");
  
  //has the database been parsed?
  if(database->is_parsed){
    fprintf(file, "TRUE\n");
    DATABASE_PROTEIN_ITERATOR_T* iterator = new_database_protein_iterator(database);
 
    while(database_protein_iterator_has_next(iterator)){
      protein = database_protein_iterator_next(iterator);
      //if the database uses light/heavy functionality
      if(database->use_light_protein){
        protein_to_heavy(protein);
      }
      print_protein(protein, stdout);
      //if the database uses light/heavy functionality
      if(database->use_light_protein){
        protein_to_light(protein);
      }
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
 * IF using light_protein functionality will not read in the sequence or id.
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
  unsigned int protein_idx = 0;

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
      if(database->num_proteins == MAX_PROTEINS){
        fclose(file);
        free(new_line);
        carp(CARP_ERROR, "exceeds protein index array size");
        return FALSE;
      }
      //the new protein to be added
      new_protein = allocate_protein();

      //do not parse the protein sequence if using light/heavy functionality
      if(database->use_light_protein){
        //set light and offset
        set_protein_offset(new_protein, working_index);
        set_protein_is_light(new_protein, TRUE);
      }
      else{
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
        set_protein_is_light(new_protein, FALSE);
      }
      
      //add protein to database
      database->proteins[database->num_proteins] = new_protein;
      ++database->num_proteins;
      //set protein index, database
      set_protein_protein_idx(new_protein, database->num_proteins);
      set_protein_database(new_protein, database);
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
 * sets the use_light_protein of the database
 */
void set_database_use_light_protein(
  DATABASE_T* database, ///< the database to set it's fields -out
  BOOLEAN_T use ///< should I use the light/heavy functionality?
  )
{
  database->use_light_protein = use;
}

/**
 *\returns TRUE|FALSE whether the database uses light/heavy
 */
BOOLEAN_T get_database_use_light_protein(
  DATABASE_T* database ///< the query database -in 
  )
{
  return database->use_light_protein;
}

/**
 *\returns the total number of proteins of the database
 */
unsigned int get_database_num_proteins(
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

/**
 *\returns the nth protein of the database
 */
PROTEIN_T* get_database_protein_at_idx(
  DATABASE_T* database, ///< the query database -in 
  unsigned int protein_idx ///< The index of the protein to retrieve -in
  )
{
  return database->proteins[protein_idx-1];
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
/*
PROTEIN_T* database_protein_iterator_protein_idx(
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator, ///< the iterator of interest -in
  unsigned int protein_idx ///< protein_idx to which protein to return -in
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
*/

/***********************************************
 * database_peptide_Iterators - can use the light protein functionality to save space
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
  PROTEIN_T* next_protein = NULL;
  
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

    //if using light/heavy functionality parse the light protein
    if(database->use_light_protein && get_protein_is_light(next_protein)){
      if(!protein_to_heavy(next_protein)){
        carp(CARP_FATAL, "failed to create a database_peptide_iterator, no proteins in database");
        free_database_protein_iterator(database_peptide_iterator->database_protein_iterator);
        free(database_peptide_iterator);
        exit(1);
      }
    }

    //set new protein peptide iterator
    database_peptide_iterator->cur_protein_peptide_iterator =
      new_protein_peptide_iterator(next_protein, database_peptide_iterator->peptide_constraint);
 
    //if first protein does not contain a match peptide, reinitailize
    while(!protein_peptide_iterator_has_next(database_peptide_iterator->cur_protein_peptide_iterator)){
      // covert the heavy back to light
      if(database->use_light_protein && !get_protein_is_light(next_protein)){
        protein_to_light(next_protein);
      }

      //end of list of peptides for database_peptide_iterator
      if(!database_protein_iterator_has_next(database_peptide_iterator->database_protein_iterator)){
        break;
      }
      else{ //create new protein_peptide_iterator for next protein
        //free old iterator
        free_protein_peptide_iterator(database_peptide_iterator->cur_protein_peptide_iterator);
        
        //get next protein
        next_protein = 
          database_protein_iterator_next(database_peptide_iterator->database_protein_iterator);

         //if using light/heavy functionality parse the light protein
        if(database->use_light_protein && get_protein_is_light(next_protein)){
          if(!protein_to_heavy(next_protein)){
            carp(CARP_FATAL, "failed to create a database_peptide_iterator, no proteins in database");
            free_database_protein_iterator(database_peptide_iterator->database_protein_iterator);
            free(database_peptide_iterator);
            exit(1);
          }
        }
        //creat new protein_peptide_iterator
        database_peptide_iterator->cur_protein_peptide_iterator =
          new_protein_peptide_iterator(next_protein, 
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
  //set the current working protein
  database_peptide_iterator->prior_protein = next_protein;
  
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
   //did you reset working protein?
  BOOLEAN_T reset = FALSE;
  
  //the ppeptide to return
  PEPTIDE_T* next_peptide =
    protein_peptide_iterator_next(database_peptide_iterator->cur_protein_peptide_iterator);
  
  DATABASE_T* database = database_peptide_iterator->database_protein_iterator->database;
  
  //reset database_peptide_iterator if needed
  while(!protein_peptide_iterator_has_next(database_peptide_iterator->cur_protein_peptide_iterator)){
    reset = TRUE;
    PROTEIN_T* next_protein = NULL; 
    
    // covert the heavy back to light
    if(database->use_light_protein && next_protein != NULL && !get_protein_is_light(next_protein)){
      protein_to_light(next_protein);
    }
    
    //end of list of peptides for database_peptide_iterator
    if(!database_protein_iterator_has_next(database_peptide_iterator->database_protein_iterator)){
      break;
    }
    else{ //create new protein_peptide_iterator for next protein
      //free old iterator
      free_protein_peptide_iterator(database_peptide_iterator->cur_protein_peptide_iterator);
      
      //get next protein
      next_protein = 
        database_protein_iterator_next(database_peptide_iterator->database_protein_iterator);
      
      //if using light/heavy functionality parse the light protein
      if(database->use_light_protein && get_protein_is_light(next_protein)){
        if(!protein_to_heavy(next_protein)){
          carp(CARP_FATAL, "failed to create a database_peptide_iterator, no proteins in database");
          free_database_protein_iterator(database_peptide_iterator->database_protein_iterator);
          free(database_peptide_iterator);
          exit(1);
        }
      }
      //creat new protein_peptide_iterator
      database_peptide_iterator->cur_protein_peptide_iterator =
        new_protein_peptide_iterator(next_protein, 
                                     database_peptide_iterator->peptide_constraint);
    }        
  }

  //are we using the light functionality?
  if(database->use_light_protein){
    //get the current working protein
    PROTEIN_T* protein_bye = get_protein_peptide_iterator_portein(database_peptide_iterator->cur_protein_peptide_iterator);
    //set first passed, shows that we extraced at least one protein since we moved on to the next protein
    if(!reset && !database_peptide_iterator->first_passed){
      database_peptide_iterator->first_passed = TRUE;
    }
    //convert prior_protein to light
    else if(!reset && database_peptide_iterator->first_passed && (protein_bye != database_peptide_iterator->prior_protein)){
      protein_to_light(database_peptide_iterator->prior_protein);
      database_peptide_iterator->prior_protein = protein_bye;
    }
  }
  
  return next_peptide;
}

/***********************************
 * database sorted peptide iterator
 ***********************************/

/**
 * Instantiates a new database_sorted_peptide_iterator from a database.
 * uses a sorted_peptide_iterator as it's engine
 * \returns a DATABASE_SORTED_PEPTIDE_ITERATOR_T object.
 */
DATABASE_SORTED_PEPTIDE_ITERATOR_T* new_database_sorted_peptide_iterator(
  DATABASE_T* database, ///< the database of interest -in
  PEPTIDE_CONSTRAINT_T* peptide_constraint, ///< the peptide_constraint to filter peptides -in
  SORT_TYPE_T sort_type, ///< the sort type for this iterator -in
  BOOLEAN_T unique ///< only return unique peptides? -in
  )
{
  //create database sorted peptide iterator
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator =
    (DATABASE_SORTED_PEPTIDE_ITERATOR_T*)mycalloc(1, sizeof(DATABASE_SORTED_PEPTIDE_ITERATOR_T));

  //create the database peptide iterator
  DATABASE_PEPTIDE_ITERATOR_T* db_peptide_iterator =
    new_database_peptide_iterator(database, peptide_constraint);

  //create a sorted peptide iterator that will sort all the peptides from db peptide iterator
  SORTED_PEPTIDE_ITERATOR_T* sorted_peptide_iterator =
    new_sorted_peptide_iterator_database(db_peptide_iterator, sort_type, unique);

  //set sorted_peptide_iterator
  database_sorted_peptide_iterator->sorted_peptide_iterator = sorted_peptide_iterator;
  
  free_database_peptide_iterator(db_peptide_iterator); //CHECK ME might wanna check this...

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
  return sorted_peptide_iterator_has_next(database_sorted_peptide_iterator->sorted_peptide_iterator);
}

/**
 * returns each peptide in sorted order
 * \returns The next peptide in the database.
 */
PEPTIDE_T* database_sorted_peptide_iterator_next(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator ///< the iterator of interest -in
  )
{
  return sorted_peptide_iterator_next(database_sorted_peptide_iterator->sorted_peptide_iterator);
}

/**
 * Frees an allocated database_sorted_peptide_iterator object.
 */
void free_database_sorted_peptide_iterator(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator ///< the iterator to free -in
  )
{
  free_sorted_peptide_iterator(database_sorted_peptide_iterator->sorted_peptide_iterator);
  free(database_sorted_peptide_iterator);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

