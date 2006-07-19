/*****************************************************************************
 * \file database.c
 * $Revision: 1.6 $
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

#define MAX_PROTEINS 10000 ///< The maximum number of proteins in a database.

/**
 * \struct database
 * \brief A database of protein sequences
 */
struct database{
  char*          filename;      ///< Original database filename.
  FILE_T*        file;          ///< Open filehandle for this database.
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
  DATABASE_PROTEIN_ITERATOR_T* 
  database_protein_iterator; ///<The protein iterator. 
  PROTEIN_PEPTIDE_ITERATOR_T* 
  cur_protein_peptide_iterator; ///< The peptide iterator for the current protein.
  PEPTIDE_CONSTRAINT_T* peptide_constraint; ///< The constraints for the kind of peptide to iterate over.
  };

// FIXME I think all of these fields are necessary? But maybe some can go
// in the PROTEIN_PEPTIDE_ITERATOR?



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
  fclose(databse->file);
  
  for(; protein_idx < database->num_proteins; ++protein_idx){
    free_protein(database->proteins[protein_idx]);
  }
  free(database);
}

/**
 * Prints a database object to file.
 */
void print_database(DATABASE_T* database, FILE* file);

// START HERE

/**
 * Copies database object src to dest.
 */
void copy_database(
  DATABASE_T* src,
  DATABASE_T* dest);


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
  
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  PROTEIN_T* new_protein;

  //check if already parsed
  if(database->is_parsed){
    return TRUE;
  }

  //set memory for all proteins
  database->proteins = (PROTEIN_T*)mycalloc(MAX_PROTEINS, sizeof(PROTEIN_T));
  

  //open file and 
  FILE* file = fopen(database->filename, "r");

  working_index = ftell(file);
  //check each line until reach 'S' line
  while((line_length =  getline(&new_line, &buf_length, file)) != -1){
    if(new_line[0] == '>'){
      if(database->num_proteins == MAX_PROTEINS+1){
        fclose(file);
        free(new_line);
        fprintf(stderr, "ERROR: exceeds protein index array size\n");
        return FALSE;
      }
      
      new_protein = allocate_protein();

      //rewind to the begining of the protein to include ">" line
      fseek(file, working_idex, SEEK_SET);
      
      //failed to parse the protein from fasta file
      if(!parse_protein_fasta_file(new_protein ,file)){
        fclose(file);
        free_protein(new_protein);
        for(; protein_idx < database->num_proteins; ++protein_idx){
          free_protein(database->proteins[protein_idx]);
        }
        database->num_proteins = 0;
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
  free(new_line);
  database->file = file;
  return TRUE;
}


/**
 * \returns FALSE if database has not yet been parsed or if the nth protein
 * cannot be parsed.
 */
BOOLEAN_T get_database_protein_at_idx(
    DATABASE_T* database, ///< A parsed database object -in
    int protein_idx,      ///< The index of the protein to retrieve -in
    PROTEIN_T* protein   ///< A pointer to a pointer to a PROTEIN object -out
    );




/**
 *\returns the filename of the database
 * returns a heap allocated new copy of the filename
 * user must free the return filename
 */
char* get_database_filename(
  DATABASE_T* database ///< the query database -in 
  )
{
  return copy_string(database->filename);
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
  database->filename = copy_string(filename);
}


/**
 * \struct database
 * \brief A database of protein sequences
 */
struct database{
  char*          filename;      ///< Original database filename.
  FILE_T*        file;          ///< Open filehandle for this database.
                                ///  A database has only one
                                ///  associated file.
  int num_proteins;             ///< Number of proteins in this database.
  BOOLEAN_T is_parsed;          ///< Has this database been parsed yet.
  PROTEIN_T* proteins[MAX_PROTEINS];   ///< Proteins in this database (not yet needed.)
};    


/**
 *\returns the FILE* of the database
 * returns a heap allocated new copy of the filename
 * user must free the return filename
 */
char* get_database_filename(
  DATABASE_T* database ///< the query database -in 
  )
{
  return copy_string(database->filename);
}

/**
 * sets the filename of the database
 * protein->sequence must been initiailized
 */
void set_database_filename(
  DATABASE_T* database, ///< the database to set it's fields -out
  char* filename ///< the filename to add -in
  )


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

