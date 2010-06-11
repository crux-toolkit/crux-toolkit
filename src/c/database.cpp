/*************************************************************************//**
 * \file database.cpp
 * \brief Object for representing a database of protein sequences.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include "utils.h"
#include "crux-utils.h"
#include "peptide.h"
#include "protein.h"
#include "database.h"
#include "hash.h"
#include "carp.h"
#include "objects.h"
#include "peptide_constraint.h"
#include "sorter.h"
#include "protein_index.h"

#include <map>

using namespace std;

#define MAX_PROTEINS 3300000 ///< The maximum number of proteins in a database.

//Comparator function for c type strings.
struct cmp_str {

  bool operator()(char const *a, char const *b) {
    return strcmp(a, b) < 0;
  }
};

/**
 * \struct database
 * \brief A database of protein sequences
 */
struct database{
  char*        filename; ///< Original database filename.
  FILE*        file;     ///< Open filehandle for this database.
                         ///  A database has only one associated file.
  unsigned int num_proteins; ///< Number of proteins in this database.
  BOOLEAN_T is_parsed;  ///< Has this database been parsed yet.
  PROTEIN_T* proteins[MAX_PROTEINS]; ///< Proteins in this database. 
  map<char*, PROTEIN_T*, cmp_str> protein_map; //map for proteins 
  BOOLEAN_T is_hashed; //Indicator of whether the database has been hashed/mapped.
  unsigned long int size; ///< The size of the database in bytes (convenience)
  BOOLEAN_T use_light_protein; ///< should I use the light/heavy protein option
  BOOLEAN_T is_memmap; ///< Are we using a memory mapped fasta file? 
  void* data_address; ///< pointer to the beginning of the memory mapped data, 
  unsigned int pointer_count; ///< number of pointers referencing this database. 
  long file_size; ///< the size of the binary fasta file, when memory mapping
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
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator; 
    ///< The protein iterator. 
  PROTEIN_PEPTIDE_ITERATOR_T* cur_protein_peptide_iterator; 
    ///< The peptide iterator for the current protein.
  PEPTIDE_CONSTRAINT_T* peptide_constraint; 
    ///< The constraint for the kind of peptide to iterate over.
  PROTEIN_T* prior_protein; 
    ///< the protein that was used before the current working protein
  BOOLEAN_T first_passed; 
    ///< is it ok to convert prior_protein to light?
  PEPTIDE_T* cur_peptide; ///< the peptide to return by next()
};

/**
 * \struct database_sorted_peptide_iterator
 * \brief Object to iterate over the peptides within a database, in an
 * specified sorted order.(mass, length, lexical)
 */
struct database_sorted_peptide_iterator {
  SORTED_PEPTIDE_ITERATOR_T* sorted_peptide_iterator; 
    ///< the peptide iterator that sorts the peptides
};

/* Private Functions */
PEPTIDE_T* database_peptide_iterator_next_from_file(
  DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator
  );
BOOLEAN_T database_peptide_iterator_has_next_from_file(
  DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator
  );
void generate_all_peptides(DATABASE_PEPTIDE_ITERATOR_T* iter);


/**
 * \returns An (empty) database object.
 */
DATABASE_T* allocate_database(void){
  DATABASE_T* database = (DATABASE_T*)mycalloc(1,sizeof(DATABASE_T));
  database->filename = NULL;
  database->file = NULL;
  database->num_proteins = 0; 
  database->is_parsed = FALSE;
  database->size = 0; 
  database->use_light_protein = FALSE; 
  database->is_memmap = FALSE;
  database->data_address = NULL;
  database->pointer_count = 1;
  database->file_size = 0;
  database->is_hashed = FALSE;
  database->protein_map =  map<char*, PROTEIN_T*, cmp_str>();
  // fprintf(stderr, "Free: Allocation: %i\n", database->pointer_count);
  return database;
}

/**
 * \returns A new database object.
 */
DATABASE_T* new_database(
  const char*         filename, ///< The file from which to parse the database. 
  ///< either text fasta file or binary fasta file -in
  BOOLEAN_T is_memmap ///< are we using a memory mapped binary fasta file? 
  ///< If so, all proteins are memory mapped -in
  )         
{
  carp(CARP_DEBUG, "Creating new database from '%s'", filename);
  DATABASE_T* database = allocate_database();
  set_database_filename(database, filename);
  database->is_memmap = is_memmap;
  database->num_proteins = 0;

  return database;
}  

/**
 * Frees an allocated protein object.
 */
void free_database(
  DATABASE_T* database ///< An allocated database -in
  )
{
  
  // decrement database pointer counter
  --database->pointer_count;
  carp(CARP_DETAILED_DEBUG, "Database pointer count %i",  
      database->pointer_count);

  // DEBUG show the databse pointer count
  // printf("Free: After free: %s: %d\n", database->pointer_count);

  // only free up memory when pointer count is zero
  if(database->pointer_count > 0){ 
    // TODO maybe change this to number of proteins? 
    // since I think we now have circular references
    return;
  }
  
  free(database->filename);
  
  // only free proteins if been parsed and file has been opened
  if(database->is_parsed){
    carp(CARP_INFO, "Freeing database.");
    
    // free each protein in the array
    unsigned int protein_idx;
    for(protein_idx=0; protein_idx < database->num_proteins; ++protein_idx){
      free_protein(database->proteins[protein_idx]);
    }
    // free(database->proteins);
    
    // free memory mapped binary file from memory
    if(database->is_memmap){
      // un map the memory!!
      if(munmap(database->data_address, database->file_size) != 0){
        carp(CARP_ERROR, "failed to unmap the memory of binary fasta file");
      }
    }
    // not memory mapped
    else{
      // close file handle
      carp(CARP_INFO, "Closing database filehandle");
      fclose(database->file);
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
  
  // has the database been parsed?
  if(database->is_parsed){
    fprintf(file, "TRUE\n");
    DATABASE_PROTEIN_ITERATOR_T* iterator 
      = new_database_protein_iterator(database);
 
    while(database_protein_iterator_has_next(iterator)){
      protein = database_protein_iterator_next(iterator);
      // if the database uses light/heavy functionality
      if(database->use_light_protein){
        protein_to_heavy(protein);
      }
      print_protein(protein, stdout);
  
      // if the database uses light/heavy functionality
      /** 
       * uncomment this code if you want to restore a protein to 
       * light after converted to heavy
      if(database->use_light_protein){
        protein_to_light(protein);
      }
      */
    }
    free_database_protein_iterator(iterator);
  }
  else{
    fprintf(file, "FALSE\n");
  }
}

/**
 * Parses a database from the text based fasta file in the filename
 * member variable
 * reads in all proteins in the fasta file and creates a protein object
 * and adds them to the database protein array
 * total proteins in fasta file must not exceed MAX_PROTEIN constant
 * IF using light_protein functionality will not read in the sequence or id.
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_database_text_fasta(
  DATABASE_T* database ///< An allocated database -in
  )
{
  unsigned long working_index;
  FILE* file = NULL;
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  PROTEIN_T* new_protein;
  unsigned int protein_idx;

  carp(CARP_DEBUG, "Parsing text fasta file '%s'", database->filename);
  // check if already parsed
  if(database->is_parsed){
    return TRUE;
  }
  
  // open file and 
  file = fopen(database->filename, "r");
  
  // check if succesfully opened file
  if(file == NULL){
    carp(CARP_ERROR, "Failed to open fasta file %s", database->filename);
    return FALSE;
  }
  
  // check if use light protein and parse thos light proteins fomr protein index
  if(database->use_light_protein && protein_index_on_disk(database->filename, FALSE)){
    // let the user know that protein index file is being used
    carp(CARP_INFO, "using protein index file");

    // create a protein index iterator
    PROTEIN_INDEX_ITERATOR_T* protein_index_iterator =
      new_protein_index_iterator(database->filename);

    // iterate over all proteins in protein index
    while(protein_index_iterator_has_next(protein_index_iterator)){
      // check if there's space for more proteins
      if(database->num_proteins == MAX_PROTEINS){
        free_protein_index_iterator(protein_index_iterator);
        carp(CARP_ERROR, "exceeds protein index array size");
        return FALSE;
      }
      
      new_protein = protein_index_iterator_next(protein_index_iterator);
      set_protein_database(new_protein, database);
      
      // add protein to database
      database->proteins[database->num_proteins] = new_protein;
      ++database->num_proteins;
    }
    // job well done..free iterator
    free_protein_index_iterator(protein_index_iterator);
  }
  else{  
    working_index = ftell(file);
    // check each line until reach '>' line
    while((line_length =  getline(&new_line, &buf_length, file)) != -1){
      if(new_line[0] == '>'){
        if(database->num_proteins == MAX_PROTEINS){
          fclose(file);
          free(new_line);
          carp(CARP_ERROR, "exceeds protein index array size");
          return FALSE;
        }
        // the new protein to be added
        new_protein = allocate_protein();
        
        // do not parse the protein sequence if using light/heavy functionality
        if(database->use_light_protein){
          // set light and offset
          set_protein_offset(new_protein, working_index);
          set_protein_is_light(new_protein, TRUE);
        }
        else{
          // rewind to the beginning of the protein to include ">" line
          fseek(file, working_index, SEEK_SET);
          
          // failed to parse the protein from fasta file
          // protein offset is set in the parse_protein_fasta_file method
          if(!parse_protein_fasta_file(new_protein ,file)){
            fclose(file);
            free_protein(new_protein);
            for(protein_idx=0;protein_idx<database->num_proteins;protein_idx++){
              free_protein(database->proteins[protein_idx]);
            }
            database->num_proteins = 0;
            carp(CARP_ERROR, "failed to parse fasta file");
            return FALSE;
          }
          set_protein_is_light(new_protein, FALSE);
        }
        
        // add protein to database
        database->proteins[database->num_proteins] = new_protein;
        // set protein index, database
        set_protein_protein_idx(new_protein, database->num_proteins);
        set_protein_database(new_protein, database);
        ++database->num_proteins;
      }
      working_index = ftell(file);
    }
    free(new_line);
  }
  
  // yes the database is paresed now..!!
  database->is_parsed = TRUE;
  database->file = file;
  return TRUE;
}

/**
 * memory maps the binary fasta file for the database
 *\return TRUE if successfully memory map binary fasta file, else FALSE
 */
BOOLEAN_T memory_map_database(
  DATABASE_T* database, ///< An allocated database -in/out
  int file_d  ///<  file descriptor -in
  )
{
  struct stat file_info;
  
  // get information of the binary fasta file
  if (stat(database->filename, &file_info) == -1) {
    carp(CARP_ERROR,
         "Failed to retrieve information of binary fasta file: %s",
         database->filename);
    return FALSE;
  }
  
  // set size of the binary fasta file in database
  // this is used later to know how much to unmap
  database->file_size = file_info.st_size;
  
  // memory map the entire binary fasta file!
  database->data_address = mmap((caddr_t)0, file_info.st_size, PROT_READ, MAP_PRIVATE /*MAP_SHARED*/, file_d, 0);

  // check if memory mapping has succeeded
  if ((caddr_t)(database->data_address) == (caddr_t)(-1)){
    carp(CARP_ERROR, "failed to use mmap function for binary fasta file: %s", database->filename);
    return FALSE;
  }
  
  return TRUE;
}

/**
 * Assumes that there is a 1 at the very end after all the proteins in binary file
 *\return TRUE successfully populates the proteins from memory mapped binary fasta file, else FALSE
 */
BOOLEAN_T populate_proteins_from_memmap(
  DATABASE_T* database ///< An allocated database -in/out
  )
{
  PROTEIN_T* new_protein;
  unsigned int protein_idx = 0;
  char* data = (char*)database->data_address;
  
  // parse proteins until the end of list
  while((int)data[0] != 1){
    // check if anymore space for protein
    if(database->num_proteins == MAX_PROTEINS){
      carp(CARP_ERROR, "exceeds protein index array size");
      // free all proteins before return
      for(protein_idx = 0; protein_idx < database->num_proteins; ++protein_idx){
        free_protein(database->proteins[protein_idx]);
      }
      database->num_proteins = 0;
      return FALSE;
    }
    
    // the new protein to be added
    new_protein = allocate_protein();
    
    // parse protein from memory map
    if(!parse_protein_binary_memmap(new_protein, &data)){
      // failed to parse the protein from memmap
      // free all proteins, and return FALSE
      free_protein(new_protein);
      for(; protein_idx < database->num_proteins; ++protein_idx){
        free_protein(database->proteins[protein_idx]);
      }
      database->num_proteins = 0;
      carp(CARP_ERROR, "failed to parse fasta file");
      return FALSE;
    }
    set_protein_is_light(new_protein, FALSE);
    
    // add protein to database
    database->proteins[database->num_proteins] = new_protein;
    // set protein index, database
    set_protein_protein_idx(new_protein, database->num_proteins);
    set_protein_database(new_protein, database);
    ++database->num_proteins;
  }

  return TRUE;
}

/**
 * \brief Parses a database from the binary fasta file in the filename
 * member variable.
 *
 * Memory maps the binary fasta file into memory. The protein
 * sequences are not copied, but just pointed to the memory mapped
 * location. 
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_database_memmap_binary(
  DATABASE_T* database ///< An allocated database -in
  )
{
  int file_d = -1;
  carp(CARP_DEBUG, "Parsing binary fasta file '%s'", database->filename);
 
  // check if already parsed
  if(database->is_parsed){
    return TRUE;
  }
  
  // open file and 
  file_d = open(database->filename, O_RDONLY);
  
  // check if succesfully opened file
  if(file_d == -1){
    carp(CARP_FATAL, "Failed to open file to parse database");
    return FALSE;
  }

  // FIXME, if what to use some light protein for binary file change here...
  // check if user request light protein
  // When using a binary file in memory map, cannot use light protein
  // change to FALSE on light protein useage
  if(database->use_light_protein){
    carp(CARP_WARNING, 
         "memory mapping does not support light protein,changing settings to use heavy protein");
    database->use_light_protein = FALSE;;
  }

  // memory map the binary fasta file into memory
  if(!memory_map_database(database, file_d)){
    carp(CARP_ERROR, "Failed to memory map binary fasta file into memory");
    return FALSE;
  }

  // populate the proteins from the memory mapped fasta file
  if(!populate_proteins_from_memmap(database)){
    carp(CARP_ERROR, "Failed to populate the proteins from memory mapped fasta file");
    return FALSE;
  }
   
  // yes the database is paresed now..!!
  database->is_parsed = TRUE;
  return TRUE;
}


/**
 * Parses a database from the file in the filename member variable
 * The is_memmap field in the database struct determines whether the
 * input file is a binary fasta file or normal text fasta file.
 *
 * IF is_memmap is true, memory maps the entire binary fasta file into memory
 * and then creates protein objects that point to the memory mapped binary file
 *
 * IF is_memmap is false, uses the traditional text fasta file which
 * it parses out the various peptides for each protein. Only when
 * using text fasta file can you use light/heavy protein, in which 
 * if using light_protein functionality will not read in the sequence
 * or id. Will parse sequence if protein  
 * is needed, lazy parsing.
 *
 * For Both cases, reads in all proteins in file and creates a protein object
 * and adds them to the database protein array
 * total proteins in fasta file must not exceed MAX_PROTEIN constant
 *
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_database(
  DATABASE_T* database ///< An allocated database -in
  )
{
  // should we parse the database using memory mapped binary fasta file?
  if(database->is_memmap){
    return parse_database_memmap_binary(database);   
  }
  else{ // parse database from normal text fasta file, no memory mapping!
    return parse_database_text_fasta(database);
  }
  
  // succeeded in parsing database
  return TRUE;
}


/**
 * \brief Changes a database from one that reads from a fasta file to
 * one that reads from a binary/memmory mapped protein file.
 *
 * If database already has binary source (i.e. is_memmap == TRUE), 
 * returns TRUE.  
 * Opens the fasta file pointed to by filename for reading.  Creates an
 * output file with the name given.  Reads in each protein from the
 * text file and serializes it to the output file.  Closes both files.
 * Changes filename to point to new output file and sets is_memmap to
 * true. Parses the database.
 * \returns TRUE if all processes succeed, else FALSE.
 */
BOOLEAN_T transform_database_text_to_memmap(
  DATABASE_T* database,
  char* output_dir
  ){

  BOOLEAN_T success = FALSE;

  // from output_dir name and database filename, get binary_filename
  char* binary_filename = generate_name_path( database->filename, ".fasta",
                                           "-binary-fasta", output_dir);


  carp(CARP_DEBUG, "Transforming text file '%s' to binary file '%s'",
       database->filename, binary_filename);

  // create binary fasta
  success = create_binary_fasta_here(database->filename, 
                                     binary_filename);

  if(! success ){
    carp(CARP_ERROR, 
         "Could not create binary fasta file '%s' from text fasta file '%s'", 
         binary_filename, database->filename);
    return FALSE;
  }
  // change database filename to new binary fasta
  char* binary_filename_no_path = parse_filename(binary_filename);
  set_database_filename(database, binary_filename);

  // set is_memmap to true
  set_database_memmap(database, TRUE);

  // parse the binary fasta
  success = parse_database(database);

  free(binary_filename);
  free(binary_filename_no_path);
  return success;
}


/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 */

/**
 *sets TRUE,FALSE whether the database uses memory mapped
 */
void set_database_memmap(
  DATABASE_T* database, ///< the query database -in 
  BOOLEAN_T is_memmap  ///< is the database memory mapped?
  )
{
  database->is_memmap = is_memmap;
}

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
 *\returns the pointer to the filename of the database
 * user must not free or change the filename
 */
char* get_database_filename_pointer(
  DATABASE_T* database ///< the query database -in 
  )
{
  return database->filename;
}

/**
 * sets the filename of the database
 * protein->sequence must been initiailized
 */
void set_database_filename(
  DATABASE_T* database, ///< the database to set it's fields -out
  const char* filename ///< the filename to add -in
  )
{
  free(database->filename);
  // MEMLEAK below
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
 * \returns the nth protein of the database
 * 
 */
PROTEIN_T* get_database_protein_at_idx(
  DATABASE_T* database,    ///< the query database -in 
  unsigned int protein_idx ///< The index of the protein to retrieve -in
  )
{
  //carp(CARP_DETAILED_DEBUG, "Getting db protein idx = %i, num proteins %i", 
  //     protein_idx, database->num_proteins);
  if( protein_idx >= database->num_proteins ){
    carp(CARP_FATAL, 
         "Protein index %i out of bounds.  %i proteins in the database",
         protein_idx, database->num_proteins);
  }

  return database->proteins[protein_idx];
}

/**
 *\returns the protein designated by protein id of the database
 */
PROTEIN_T* get_database_protein_by_id_string(
  DATABASE_T* database, ///< the query database -in
  const char* protein_id ///< The id string for this protein -in
  ) {

  //TODO - Implement as a hashtable rather than a map to make 
  //this even faster if needed.
  PROTEIN_T* protein = NULL;
  if (database->is_hashed) {
    map<char*, PROTEIN_T*>::iterator find_iter;
    find_iter = database->protein_map.find((char*)protein_id);

    if (find_iter != database->protein_map.end()) {
      protein = find_iter->second;
    }
  } else {
    //create the hashtable of protein ids
    for (unsigned int protein_idx = 0;
      protein_idx < database->num_proteins;
      protein_idx++) {

      PROTEIN_T* current_protein = database->proteins[protein_idx];
      char* current_id = get_protein_id_pointer(current_protein);
      database->protein_map[current_id] = current_protein;

      if (strcmp(current_id, protein_id)==0) {
        protein = current_protein;
      }
        
    }
    database->is_hashed = TRUE;
  }
  return protein;
}

/**
 * increase the pointer_count produced by this database.
 * \returns database pointer
 */
DATABASE_T* copy_database_ptr(
  DATABASE_T* database ///< the query database -in/out
  )
{
  if( database == NULL ){
    return NULL;
  }
  ++database->pointer_count;
  return database;
}




/***********************************************
 * Iterators
 ***********************************************/

/**
 * Instantiates a new database_protein_iterator from a database.
 * 
 * \returns a DATABASE_PROTEIN_ITERATOR_T object.
 */
DATABASE_PROTEIN_ITERATOR_T* new_database_protein_iterator(
  DATABASE_T* database ///< the database to create a protein iterator -in
  )
{
  // if database is parsed, if not do so..
  if(!database->is_parsed){
    // failed to parse database
    if(!parse_database(database)){
      carp(CARP_FATAL, "Failed to parse database, cannot create iterator");
    }
  }
  
  // create new protein iterator
  DATABASE_PROTEIN_ITERATOR_T* iterator = (DATABASE_PROTEIN_ITERATOR_T*)
    mycalloc(1, sizeof(DATABASE_PROTEIN_ITERATOR_T));
  iterator->database = copy_database_ptr(database);
  iterator->cur_protein = 0;

  return iterator;
}        


/**
 * Frees an allocated database_protein_iterator object.
 */
void free_database_protein_iterator(
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator 
    ///< the iterator to free -in
  )
{
  // subtract pointer count
  free_database(database_protein_iterator->database);
  free(database_protein_iterator);
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional proteins to iterate over, FALSE if not.
 */
BOOLEAN_T database_protein_iterator_has_next(
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator 
    ///< the query iterator -in
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
  DATABASE_PROTEIN_ITERATOR_T* database_protein_iterator  
    ///< the query iterator -in
  )
{
  ++database_protein_iterator->cur_protein;

  // print number of protein generated to STDERR for every 500 protein reached
  if(database_protein_iterator->cur_protein % 500 == 0){
    carp(CARP_DETAILED_DEBUG, "Reached protein %d out of %d", 
         database_protein_iterator->cur_protein,
         database_protein_iterator->database->num_proteins);
  }
  
  return database_protein_iterator->database->proteins[database_protein_iterator->cur_protein-1];
}

/***********************************************
 * database_peptide_Iterators - can use the light protein functionality 
 * to save space
 ***********************************************/

/**
 * Instantiates a new database_peptide_iterator from a database.
 * \returns a DATABASE_PEPTIDE_ITERATOR_T object.
 */
DATABASE_PEPTIDE_ITERATOR_T* new_database_peptide_iterator(
  DATABASE_T* database, 
    ///< the database of interest -in
  PEPTIDE_CONSTRAINT_T* peptide_constraint 
    ///< the peptide_constraint with which to filter peptides -in
  )
{
  // set peptide implementation to linklist peptide_src
  // this determines which peptide free method to use
  set_peptide_src_implementation(TRUE);

  PROTEIN_T* next_protein = NULL;
  
  DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator =
    (DATABASE_PEPTIDE_ITERATOR_T*)
    mycalloc(1, sizeof(DATABASE_PEPTIDE_ITERATOR_T));
  
  // set a new protein iterator
  database_peptide_iterator->database_protein_iterator =
    new_database_protein_iterator(database);
  if( database_peptide_iterator->database_protein_iterator == NULL){
    carp(CARP_ERROR, 
         "Could not create protein iterator for database peptide iterator.");
    return NULL;
  }

  // set peptide constraint
  database_peptide_iterator->peptide_constraint 
    = copy_peptide_constraint_ptr(peptide_constraint);
  
  // check if there are any proteins to create peptides from
  if(database_protein_iterator_has_next(
        database_peptide_iterator->database_protein_iterator)){

    next_protein = database_protein_iterator_next(
        database_peptide_iterator->database_protein_iterator);

    // if using light/heavy functionality parse the light protein
    if(database->use_light_protein && get_protein_is_light(next_protein)){
      if(!protein_to_heavy(next_protein)){
        carp(CARP_FATAL, "failed to create a database_peptide_iterator,"
                         "no proteins in database");
      }
    }

    // set new protein peptide iterator
    database_peptide_iterator->cur_protein_peptide_iterator =
      new_protein_peptide_iterator(next_protein, peptide_constraint);
 
    // if first protein does not contain a match peptide, reinitailize
    while(!protein_peptide_iterator_has_next(
          database_peptide_iterator->cur_protein_peptide_iterator)){
      // covert the heavy back to light
      /** 
       * uncomment this code if you want to restore a protein to 
       * light after converted to heavy
      if(database->use_light_protein && !get_protein_is_light(next_protein)){
        protein_to_light(next_protein);
      }
      */

      // end of list of peptides for database_peptide_iterator
      if(!database_protein_iterator_has_next(
            database_peptide_iterator->database_protein_iterator)){
        break;
      }
      else{ // create new protein_peptide_iterator for next protein
        // free old iterator
        free_protein_peptide_iterator(
            database_peptide_iterator->cur_protein_peptide_iterator);
        
        // get next protein
        next_protein = database_protein_iterator_next(
            database_peptide_iterator->database_protein_iterator);

         // if using light/heavy functionality parse the light protein
        if(database->use_light_protein && get_protein_is_light(next_protein)){
          if(!protein_to_heavy(next_protein)){
            carp(CARP_FATAL, "failed to create a database_peptide_iterator"
                              " no proteins in database");
          }
        }
        // create new protein_peptide_iterator
        database_peptide_iterator->cur_protein_peptide_iterator =
          new_protein_peptide_iterator(next_protein, peptide_constraint);
      }
    }
  }
  else{ // no proteins to create peptides from
    carp(CARP_FATAL, "failed to create a database_peptide_iterator,"
                     "no proteins in database");
  }
  // set the current working protein
  database_peptide_iterator->prior_protein = next_protein;
  
  database_peptide_iterator->cur_peptide = 
      database_peptide_iterator_next_from_file(database_peptide_iterator);

  return database_peptide_iterator;
}

/**
 * Frees an allocated database_peptide_iterator object.
 */
void free_database_peptide_iterator(
  DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator 
  ///< the iterator to free -in
  )
{
  free_protein_peptide_iterator(
                     database_peptide_iterator->cur_protein_peptide_iterator);
  free_database_protein_iterator(
                     database_peptide_iterator->database_protein_iterator);
  free_peptide_constraint(database_peptide_iterator->peptide_constraint);
  free(database_peptide_iterator);

}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over,
 * FALSE if not. 
 */
BOOLEAN_T database_peptide_iterator_has_next_from_file(
  DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator  
  ///< the iterator of interest -in
  )
{
  if(protein_peptide_iterator_has_next(database_peptide_iterator->cur_protein_peptide_iterator)){ 
    return TRUE;
  }
  return FALSE;
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over,
 * FALSE if not. 
 */
BOOLEAN_T database_peptide_iterator_has_next(
  DATABASE_PEPTIDE_ITERATOR_T* iter ///< the iterator of interest -in
  )
{
  if( iter == NULL ){
    return FALSE;
  }

  if( iter->cur_peptide == NULL ){
    return FALSE;
  } else {
    return TRUE;
  }
}

PEPTIDE_T* database_peptide_iterator_next(
  DATABASE_PEPTIDE_ITERATOR_T* iter){

  if( iter == NULL ){
    return NULL;
  }
  
  PEPTIDE_T* return_peptide = iter->cur_peptide;
  if( database_peptide_iterator_has_next_from_file(iter) ){
    iter->cur_peptide = database_peptide_iterator_next_from_file(iter);
  } else {
    iter->cur_peptide = NULL;
  }
  return return_peptide;
}

/**
 * \returns The next peptide in the database.
 */
PEPTIDE_T* database_peptide_iterator_next_from_file(
//PEPTIDE_T* database_peptide_iterator_next(
  DATABASE_PEPTIDE_ITERATOR_T* database_peptide_iterator
  ///< the iterator of interest -in
  )
{
  /*BF: Could this be simplified?  if next peptide, return it
    if not, look for next protein, if not return NULL
  */

   // did you reset working protein?
  BOOLEAN_T reset = FALSE;
  
  // the peptide to return
  PEPTIDE_T* next_peptide =
    protein_peptide_iterator_next(
                   database_peptide_iterator->cur_protein_peptide_iterator);
  
  DATABASE_T* database = 
    database_peptide_iterator->database_protein_iterator->database;
  
  // reset database_peptide_iterator if needed
  while(!protein_peptide_iterator_has_next(
               database_peptide_iterator->cur_protein_peptide_iterator)){
    reset = TRUE;
    PROTEIN_T* next_protein = NULL; 
    
    /** 
     * uncomment this code if you want to restore a protein to 
     * light after converted to heavy
    // covert the heavy back to light
    if(database->use_light_protein && next_protein != NULL && !get_protein_is_light(next_protein)){
      protein_to_light(next_protein);
    }
    */

    // end of list of peptides for database_peptide_iterator
    if(!database_protein_iterator_has_next(database_peptide_iterator->database_protein_iterator)){
      break;
    }
    else{ // create new protein_peptide_iterator for next protein
      // free old iterator
      free_protein_peptide_iterator(database_peptide_iterator->cur_protein_peptide_iterator);
      
      // get next protein
      next_protein = 
        database_protein_iterator_next(database_peptide_iterator->database_protein_iterator);
      
      // if using light/heavy functionality parse the light protein
      if(database->use_light_protein && get_protein_is_light(next_protein)){
        if(!protein_to_heavy(next_protein)){
          carp(CARP_FATAL, "Failed to create a database_peptide_iterator, " 
                            "no proteins in database");
        }
      }
      // create new protein_peptide_iterator
      database_peptide_iterator->cur_protein_peptide_iterator =
        new_protein_peptide_iterator(next_protein, 
                                     database_peptide_iterator->peptide_constraint);
    }        
  }

  // are we using the light functionality?
  if(database->use_light_protein){
    // get the current working protein
    PROTEIN_T* protein_bye = get_protein_peptide_iterator_portein(database_peptide_iterator->cur_protein_peptide_iterator);
    // set first passed, shows that we extraced at least one protein since we moved on to the next protein
    if(!reset && !database_peptide_iterator->first_passed){
      database_peptide_iterator->first_passed = TRUE;
    }
    // convert prior_protein to light
    else if(!reset && database_peptide_iterator->first_passed && (protein_bye != database_peptide_iterator->prior_protein)){
      /** 
       * uncomment this code if you want to restore a protein to 
       * light after converted to heavy
      protein_to_light(database_peptide_iterator->prior_protein);
      */
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
  PEPTIDE_CONSTRAINT_T* peptide_constraint, 
    ///< the peptide_constraint to filter peptides -in
  SORT_TYPE_T sort_type, ///< the sort type for this iterator -in
  BOOLEAN_T unique ///< only return unique peptides? -in
  )
{
  // create database sorted peptide iterator
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator =
    (DATABASE_SORTED_PEPTIDE_ITERATOR_T*)
    mycalloc(1, sizeof(DATABASE_SORTED_PEPTIDE_ITERATOR_T));

  // create the database peptide iterator
  DATABASE_PEPTIDE_ITERATOR_T* db_peptide_iterator =
    new_database_peptide_iterator(database, peptide_constraint);

  // create a sorted peptide iterator from db peptide iterator
  SORTED_PEPTIDE_ITERATOR_T* sorted_peptide_iterator =
    new_sorted_peptide_iterator_database(
        db_peptide_iterator, sort_type, unique);

  // set sorted_peptide_iterator
  database_sorted_peptide_iterator->sorted_peptide_iterator 
    = sorted_peptide_iterator;
  
  free_database_peptide_iterator(db_peptide_iterator); 

  return database_sorted_peptide_iterator;
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides, FALSE if not.
 */
BOOLEAN_T database_sorted_peptide_iterator_has_next(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator 
    ///< the iterator of interest -in
  )
{
  return sorted_peptide_iterator_has_next(
      database_sorted_peptide_iterator->sorted_peptide_iterator);
}

/**
 * returns each peptide in sorted order
 * \returns The next peptide in the database.
 */
PEPTIDE_T* database_sorted_peptide_iterator_next(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator 
    ///< the iterator of interest -in
  )
{
  return sorted_peptide_iterator_next(
      database_sorted_peptide_iterator->sorted_peptide_iterator);
}

/**
 * Frees an allocated database_sorted_peptide_iterator object.
 */
void free_database_sorted_peptide_iterator(
  DATABASE_SORTED_PEPTIDE_ITERATOR_T* database_sorted_peptide_iterator 
    ///< the iterator to free -in
  )
{
  free_sorted_peptide_iterator(
      database_sorted_peptide_iterator->sorted_peptide_iterator);
  free(database_sorted_peptide_iterator);
}


/**********************************************************************
 * wrapper, for generate_peptides_iterator, cast back to original type
 ***********************************************************************/

/**
 * Frees an allocated database_peptide_iterator object.
 */
void void_free_database_peptide_iterator(
  void* database_peptide_iterator ///< the iterator to free -in
  )
{
  free_database_peptide_iterator(
      (DATABASE_PEPTIDE_ITERATOR_T*)database_peptide_iterator);
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides, FALSE if not.
 */
BOOLEAN_T void_database_peptide_iterator_has_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  )
{
  return database_peptide_iterator_has_next(
      (DATABASE_PEPTIDE_ITERATOR_T*)database_peptide_iterator);
}

/**
 * \returns The next peptide in the database.
 */
PEPTIDE_T* void_database_peptide_iterator_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  )
{
  return database_peptide_iterator_next(
      (DATABASE_PEPTIDE_ITERATOR_T*)database_peptide_iterator);
}

/**********************************************************************
 * wrapper, for generate_peptides_iterator, cast back to original type
 ***********************************************************************/

/**
 * Frees an allocated database_sorted_peptide_iterator object.
 */
void void_free_database_sorted_peptide_iterator(
  void* database_peptide_iterator ///< the iterator to free -in
  )
{
  free_database_sorted_peptide_iterator(
      (DATABASE_SORTED_PEPTIDE_ITERATOR_T*)database_peptide_iterator);
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T void_database_sorted_peptide_iterator_has_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  )
{
  return database_sorted_peptide_iterator_has_next(
      (DATABASE_SORTED_PEPTIDE_ITERATOR_T*)database_peptide_iterator);
}

/**
 * returns each peptide in sorted order
 * \returns The next peptide in the database.
 */
PEPTIDE_T* void_database_sorted_peptide_iterator_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  )
{
  return database_sorted_peptide_iterator_next(
      (DATABASE_SORTED_PEPTIDE_ITERATOR_T*)database_peptide_iterator);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

