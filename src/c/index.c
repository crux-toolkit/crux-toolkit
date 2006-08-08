/*****************************************************************************
 * \file index.c
 * $Revision: 1.5 $
 * \brief: Object for representing an index of a database
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include "utils.h"
#include "crux-utils.h"
#include "peptide.h"
#include "protein.h"
#include "index.h"
#include "carp.h"
#include "objects.h"
#include "peptide_constraint.h"
#include "database.h"


#define INDEX_MAX_FILES 10000


/* 
 * How does the create_index (the routine) work?
 *
 * - Creates a database object from a the fasta file member variable.
 *
 * - Depending on the size of the database, determines how many passes it
 *   will need to make over the database peptides (first implementation
 *   will only make one pass, later implementation can make multiple passes
 *   with multiple iterators)
 *
 * - From the directory and mass resolution member variables, 
 *   creates the list of filenames necessary for storing the peptides 
 *   
 * - Creates a database_peptide_iterator object from the database 
 *   and peptide constraint objects
 *
 * - Then starts iterating through peptides 
 *
 *    - From the peptide mass determines which file to place it in, and
 *      then calls the serialize_peptide method
 *
 *      - The serialize peptide method is an object method specifically
 *        designed for this purpose, and writes a peptide in a format from
 *        which it can be reconstructed
 *
 *        - Serialize peptide needs to serialize the peptide_src objects. To
 *        do that, it needs the idx of the peptide_src's protein objects in
 *        the database. This is retrieved from the protein idx member
 *        variable (i.e. this protein is the 0th protein, 1st protein etc.)
 *        which is set at database creation. Note, that this won't have to 
 *        change when we move to light proteins.
 *
 * LATER At some point we will index each of the peptide files too (in a TOC), 
 * so that we can retrieve them rapidly later. 
 *
 * LATER We implement light proteins and possibly an 
 * LATER create_index for protein locations in the database allowing rapid
 * creation of light protein objects.
 */

/* 
 * How does generate_peptides (the executable) with --from-index work?
 *
 * - Given a fasta file, looks for the corresponding index on disk. Fail if it
 *   can't find the corresponding index files.
 *
 * - Instantiates an index object from the fasta file with new_index.
 *
 * - Attempts to parse_index the index.
 *
 * - Instantiates an index_peptide_iterator from the index, according to
 *   constraints passed on the command line to generate_peptides.
 *
 * - Then starts iterating through peptides, which are being loaded from
 *   disk, and outputs them as before
 * 
 * LATER We implement light proteins.
 * LATER use an index for protein offsets in the database allowing rapid
 * creation of light protein objects.
 *
 * LATER Develop a conversion from light to heavy and heavy to light protein
 * objects to avoid excessive memory use.
 */

/**
 * \struct index
 * \brief A index of a database
 */
struct index{
  DATABASE_T* database; ///< The database that has been indexed.
  char* directory; ///< The directory containing the indexed files
  //char* filenames[INDEX_MAX_FILES]; ///< The files that contain the peptides
  PEPTIDE_CONSTRAINT_T* constraint; ///< Constraint which these peptides satisfy
  BOOLEAN_T on_disk; ///< Does this index exist on disk yet?
  float mass_range;  ///< the range of mass that each index file should be partitioned into
  unsigned int max_size;  ///< maximum limit of each index file
};    

/**
 * \struct index files
 * \brief A struct that contains the information of each file
 */
struct index_file{
  char* filenames;  ///< The file name that contain the peptides
  float start_mass; ///< the start mass limit in this file
  float interval;   ///< the interval of the peptides in this file
}INDEX_FILE_T;


/**
 * \returns An (empty) index object.
 */
INDEX_T* allocate_index(void){
  INDEX_T* index = (INDEX_T*)mycalloc(1, sizeof(INDEX_T));
  return index;
}

/**
 * given a fasta_file name it returns the index file directory name
 * format: myfasta_crux_index
 * \returns A heap allocated index file directory name of the given fasta file
 */
char* generate_directory_name(
  char* fasta_filename
  )
{
  int len = strlen(fasta_filename);
  int end_idx = len;
  int end_path = len;  //index of where the "." is located in the file
  char* dir_name = NULL;
  char* dir_name_tag =  "_crux_index";
  
  //cut off the ".fasta" if needed
  for(; end_idx > 0; --end_idx){
    if(strcmp(&fasta_filename[end_idx - 1], ".") == 0){
      end_path = end_idx - 1;
      break;
    }
  }
  
  dir_name = (char*)mycalloc(end_path + strlen(dir_name_tag) + 1, sizeof(char));
  strncpy(dir_name, fasta_filename, end_path);
  strcat(dir_name, dir_name_tag);
  return dir_name;
}
  


/**
 * Assumes that the fasta file is always in the ??? directory
 * \returns A new index object.
 */
INDEX_T* new_index(
  char* fasta_filename,  ///< The fasta file
  PEPTIDE_CONSTRAINT_T* constraint,  ///< Constraint which these peptides satisfy
  float mass_range,  ///< the range of mass that each index file should be partitioned into
  unsigned int max_size  ///< maximum limit of each index file
  )
{
  char** filename_and_path = NULL;
  char* working_dir = NULL;
  INDEX_T* index = allocate_index();
  DATABASE_T* database = new_database(fasta_filename);
    
  filename_and_path = parse_filename_path(fasta_filename);
  working_dir = generate_directory_name(filename_and_path[0]);
  
  //check if the index files are on disk
  //are we currently in the crux dircetory
  if(filename_and_path[1] == NULL){
    if(opendir(working_dir) != NULL){ //maybe memleak
      set_index_on_disk(index, TRUE);
    }
    else{
      set_index_on_disk(index, FALSE);
    }
  }
  else{//we are not in crux directory
    char* full_path = cat_string(filename_and_path[1],  working_dir);
    if(opendir(full_path) != NULL){ //string cat..might not work
      set_index_on_disk(index, TRUE);
    }
    else{
      set_index_on_disk(index, FALSE);
    }
    free(full_path);
  }
  
  //set each field
  set_index_directory(index, working_dir);
  set_index_constraint(index, constraint);
  set_index_database(index, database);
  set_index_mass_range(index, mass_range);
  set_index_max_size(index, max_size);

  //free filename and path string array
  free(filename_and_path[0]);
  free(filename_and_path[1]);
  free(filename_and_path);
  
  return index;
}         

/**
 * Frees an allocated index object.
 */
void free_index(
  INDEX_T* index
  )
{
  free_database(index->database);
  free(index->directory);
  free_peptide_constraint(index->constraint);
  free(index);
}


/**
 * The main index method. Does all the heavy lifting, creating files
 * serializing peptides, etc. The index directory itself should have 
 * a standard suffix (e.g. cruxidx), so that a given fasta file will have
 * an obvious index location.
 *
 * assumes that the current directory is the crux directory where the fasta file is located
 *
 * Note: create an index .info file as to map masses to files and to store 
 * all information that was used to create this index.
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T create_index(
  INDEX_T* index ///< An allocated index
  )
{
  FILE* output = NULL;
  //FILE* info_out = NULL;
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* sorted_iterator = NULL;
  PEPTIDE_T* peptide = NULL;
  int num_peptides = 0; //current number of peptides index file 
  int num_file = 1; //the ith number of index file created
  float current_mass_limit = index->mass_range;
  char* filename_tag = "crux_index_";
  char* file_num = NULL;

  //check if already created index
  if(index->on_disk){
    return TRUE;
  }
  //create temporary directory
  if(mkdir("crux_temp", S_IRWXG) !=0 &&  chdir("crux_temp") != 0){
    carp(CARP_WARNING, "cannot create temporary directory");
    return FALSE;
  }
  
  //create peptide iterator
  sorted_iterator = 
    new_database_sorted_peptide_iterator(index->database, index->constraint, MASS, TRUE);
     
  //check if any peptides are found
  if(!database_sorted_peptide_iterator_has_next(sorted_iterator)){
    carp(CARP_WARNING, "no matches found");
    return FALSE;
  }
 
  do{ 
    char* filename;

    //open the next index file
    if(num_peptides == 0){
      file_num = int_to_char(num_file);
      filename = cat_string(filename_tag, file_num);
      output = fopen(filename, "w" );
      free(file_num);
      free(filename);
    }
    
    peptide = database_sorted_peptide_iterator_next(sorted_iterator);
    
    //set the index file to the correct interval
    while(get_peptide_peptide_mass(peptide) > current_mass_limit ||
       num_peptides > index->max_size){
      fclose(output);
      ++num_file;
      num_peptides = 0;
      file_num = int_to_char(num_file);
      filename = cat_string(filename_tag, file_num);
      output = fopen(filename, "w");
      free(file_num);
      free(filename);
    }
    
    serialize_peptide(peptide, output);
    free_peptide(peptide);
    
    ++num_peptides;
      
  } //serialize the peptides into index files
  while(database_sorted_peptide_iterator_has_next(sorted_iterator));

  fclose(output);
  //free iterator
  free_database_sorted_peptide_iterator(sorted_iterator);


  chdir("..");
  //rename crux_temp to final directory name
  if(rename("crux_temp", index->directory) != 0){
    carp(CARP_WARNING, "cannot rename directory");
    return FALSE;
  }

  index->on_disk = TRUE;
  return TRUE;
}

/**
 * Does this index exist on disk?
 *
 * \returns TRUE if it does. FALSE if it does not.
 */
BOOLEAN_T index_exists(
  INDEX_T* index ///< An allocated index
  )
{
  return index->on_disk;
}

/*
 * Private methods
 */

/*
 * Returns the index filename appropriate for this peptide
 */
char* get_peptide_file_name(
    INDEX_T* index,
    PEPTIDE_T* peptide
    );



/*********************************************
 * set and get methods for the object fields
 *********************************************/

/**
 *\returns the directory of the index
 * returns a heap allocated new copy of the directory
 * user must free the return directory name
 */
char* get_index_directory(
  INDEX_T* index ///< The index -in
  )
{
  return my_copy_string(index->directory);
}

/**
 * sets the directory of the index
 * index->directory must been initiailized
 */
void set_index_directory(
  INDEX_T* index, ///< The index -in
  char* directory ///< the directory to add -in
  )
{
  free(index->directory);
  index->directory = my_copy_string(directory);
}

/**
 *\returns a pointer to the database
 */
DATABASE_T* get_index_database(
  INDEX_T* index ///< The index -in
  )
{
  return index->database;
}

/**
 * sets the database of the index
 */
void set_index_database(
  INDEX_T* index, ///< The index -in
  DATABASE_T* database ///< The database that has been indexed. -in
  )
{
  index->database = database;
}

/**
 *\returns a pointer to the peptides constraint
 */
PEPTIDE_CONSTRAINT_T* get_index_constraint(
  INDEX_T* index ///< The index -in
  )
{
  return index->constraint;
}

/**
 * sets the peptides constraint
 */
void set_index_constraint(
  INDEX_T* index, ///< The index -in
  PEPTIDE_CONSTRAINT_T* constraint ///< Constraint which these peptides satisfy -in
  )
{
  index->constraint = constraint;
}

/**
 *\returns TRUE if index files are on disk else FALSE
 */
BOOLEAN_T get_index_on_disk(
  INDEX_T* index ///< The index -in
  )
{
  return index->on_disk;
}

/**
 * sets the on disk field of index
 */
void set_index_on_disk(
  INDEX_T* index, ///< The index -in
  BOOLEAN_T on_disk ///< Does this index exist on disk yet? -in
  )
{
  index->on_disk = on_disk;
}

/**
 *\returns the range of mass that each index file should be partitioned into
 */
float get_index_mass_range(
  INDEX_T* index ///< The index -in
  )
{
  return index->mass_range;
}

/**
 * sets the mass_range field of index
 */
void set_index_mass_range(
  INDEX_T* index, ///< The index -in
  float mass_range  ///< the range of mass that each index file should be partitioned into -in
  )
{
  index->mass_range = mass_range;
}

/**
 *\returns maximum limit of each index file
 */
unsigned int get_index_max_size(
  INDEX_T* index ///< The index -in
  )
{
  return index->max_size;
}

/**
 * sets the maximum limit of each index file for the index object
 */
void set_index_max_size(
  INDEX_T* index, ///< The index -in
  unsigned int max_size  ///< maximum limit of each index file -in
  )
{
  index->max_size = max_size;
}


/**************************
 * Iterators
 **************************/

/**
 * \struct index_peptide_iterator
 * \brief An iterator to iterate over the peptides in a database
 */
struct index_peptide_iterator{
  INDEX_T* index; ///< The index object which we are iterating over
  PEPTIDE_CONSTRAINT_T* constraint; ///< The constraints which peptides should satisfy
  unsigned int peptide_idx; ///< The (non-object) index of the current peptide.
  char* current_index_filename; ///< The current file that we are reading from
  BOOLEAN_T has_next; ///< Is there another peptide?
};    

/// TODO new_index_peptide_iterator should to see if the index exists on disk

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
