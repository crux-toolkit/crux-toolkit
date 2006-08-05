/*****************************************************************************
 * \file index.c
 * $Revision: 1.3 $
 * \brief: Object for representing an index of a database
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"
#include "crux-utils.h"
#include "peptide.h"
#include "protein.h"
#include "index.h"
#include "carp.h"
#include "objects.h"
#include "peptide_constraint.h"

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
 * LATER At some point we will index each of the peptide files too (in a TOC), 
 * so that we can retrieve them rapidly later. 
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
    if(strcmp(fasta_filename[end_idx - 1], ".") == 0){
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
 * \returns A new index object.
 */
INDEX_T* new_index(
  char* fasta_filename,  ///< The fasta file
  PEPTIDE_CONSTRAINT_T* constraint,  ///< Constraint which these peptides satisfy
  float mass_range,  ///< the range of mass that each index file should be partitioned into
  unsigned int max_size  ///< maximum limit of each index file
  )
{
  char* working_dir = NULL;
  INDEX_T* index = allocate_index();
  DATABASE_T* database = new_database(fasta_filename);
  

  set_index_peptide_constraint(index, constraint);
  
  working_dir = parse_file_path(fasta_filename);

  //are we currently in the crux dircetory
  if(working_dir == NULL){
    opendir("fasta_filename" )
    //scandir(".", &namelist, 0, alphasort);
  }
  set_index_on_disk(index, FALSE);
  set_index_database(index, database);
  
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
 * Note: create an index .info file as to map masses to files and to store 
 * all information that was used to create this index.
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T create_index(
  INDEX_T* index ///< An allocated index
  )
{
  FILE* output = NULL;

  //check if already created index
  if(index->on_disk){
    return TRUE;
  }
  
  
  
  
  


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
/*
 * Iterators
 */

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
