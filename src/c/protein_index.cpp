/*************************************************************************//**
 * \file protein_index.cpp
 * \brief Object for creating a protein index or binary fasta file
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
#include "protein_index.h"


/*** type def. ***/

/**
 * \struct protein_index
 * \brief Object to store the protein relation to the fasta file
 */
struct protein_index{
  unsigned long int offset;  ///< The offset of the protein in the fasta file
  unsigned int protein_idx;   ///< The protein idx of the protein in the fasta file
};

/**
 * \struct protein_index_iterator
 * \brief Object to iterate over the proteins in the protein index file
 */
struct protein_index_iterator{
  FILE* file;  ///< The file handler of the fasta file
  PROTEIN_T* next_protein; ///< the next protein index to return
  BOOLEAN_T has_next; ///< is there a new protein to return?
};

/***************/

/**
 *
 *\returns the file handler for the output file
 */
FILE* get_output_file(
  char* fasta_file,  ///< input fasta file -in
  BOOLEAN_T is_binary_file  ///< Are we creating a binary fasta file?
  )
{
  char* name = NULL;
  
  // create a binary fasta file?
  if(is_binary_file){
    name = generate_name(fasta_file, "_binary_fasta", ".fasta", NULL);
  }
  else{// create a normal protein index file
    name = generate_name(fasta_file, "_protein_index", ".fasta", NULL);
  }
  
  FILE* file = fopen(name, "w");
  free(name);
  return file;
}


/**
 * creates a protein index on to the output_file
 * \returns TRUE if successfully creates a protein index, else false
 */
BOOLEAN_T create_protein_index(
  char* fasta_file  ///< input fasta file -in
  )
{
  unsigned long working_index;
  FILE* file = NULL;
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  unsigned int protein_idx = 0;
  FILE* output_file = NULL;

  // open file and 
  file = fopen(fasta_file, "r");

  // check if succesfully opened file
  if(file == NULL){
    carp(CARP_FATAL, "Failed to open fasta file '%s'", fasta_file);
    return FALSE;
  }

  // get output file
  output_file = get_output_file(fasta_file, FALSE);

  // check if succesfully created file
  if(output_file == NULL){
    carp(CARP_ERROR, "Failed to create protein index file");
    fclose(file);
    return FALSE;
  }

  working_index = ftell(file);
  // check each line until reach '>' line
  while((line_length =  getline(&new_line, &buf_length, file)) != -1){
    if(new_line[0] == '>'){ 
      ++protein_idx;
      fprintf(output_file, "* %d ",protein_idx);
      fprintf(output_file, "%ld\n", working_index);
    }    
    working_index = ftell(file);
  }
  
  free(new_line);
  fclose(file);
  fclose(output_file);

  return TRUE;
}

/***************************************************************
 * FIXME, implement this if needed
 *

**
 * \returns An (empty) protein_index object.
 *
PROTEIN_INDEX_T* allocate_protein_index(void);

**
 * creates a protein_index that contains the offset and protein index of the protein
 * in the fasta file.
 *\returns a new protein_index object
 *
PROTEIN_INDEX_T* new_protein_index(
  unsigned long int offset, ///< The file location in the source file in the database
  unsigned int protein_idx ///< The index of the protein in it's database.
  );

**
 *
 * free a protein index object
 *
void free_protein_index(
  PROTEIN_INDEX_T* protein_index  ///< the protein index to free
  );                        

*****************************************************************/

/**
 * input is the fasta file name which the protein index
 * should have been created.
 * or if creating binary fasta file, is that already on disk?
 *
 *\returns TRUE if protein index or binary fasta file is on disk, else FALSE
 */
BOOLEAN_T protein_index_on_disk(
  char* fasta_file, ///< input fasta file -in
  BOOLEAN_T is_binary ///< are we looking for the binary fasta file? or preotein index
  )
{
  char* name = NULL;
  
  // create a binary fasa file?
  if(is_binary){
    name = generate_name(fasta_file, "_binary_fasta", ".fasta", NULL);
  }
  else{// create a normal protein index file
    name = generate_name(fasta_file, "_protein_index", ".fasta", NULL);
  }
  
  // check if can open file
  if(access(name, F_OK)){
    free(name);
    return FALSE;
  }
  
  free(name);
  return TRUE;
}

/**
 * protein index iterator
 * the protein index iterator parses the protein index file
 * for each protein index one by one and returns a new light protein.
 */

/**
 *
 *\returns TRUE if successfully sets the protein_index_iterator, else FALSE
 */
BOOLEAN_T setup_protein_index_iterator(
  PROTEIN_INDEX_ITERATOR_T* protein_index_iterator ///< the iterator to setup -in
  )
{
  // used to parse each line from file
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  FILE* file = protein_index_iterator->file;

  // protein fields
  char star[2] = "";
  unsigned long int offset;
  unsigned int protein_idx;
  
  PROTEIN_T* protein = NULL;
  BOOLEAN_T found = FALSE;

  while((line_length =  getline(&new_line, &buf_length, file)) != -1){
    // begining of the protein feilds
    if(new_line[0] == '*'){
      // read the crux_index_file information
      if(sscanf(new_line,"%s %d %ld", 
                star, &protein_idx, &offset) < 2){
        free(new_line);
        carp(CARP_WARNING, "incorrect file format");
        fclose(file);
        return FALSE;
      }
      found = TRUE;
      break;
    }
    // skip header lines
    else if(new_line[0] == '#'){
      continue;
    }
  }
  
  // there is a next protein to return
  if(found){
    protein = new_light_protein(offset, protein_idx);
    protein_index_iterator->next_protein = protein;
    protein_index_iterator->has_next = TRUE;
  }
  // no more proteins..
  else{
    protein_index_iterator->has_next = FALSE;
  }
  
  free(new_line);
  return TRUE;
}

/**
 * There must be the correct protein index file present
 *\returns a new heap allocated protein index iterator
 */
PROTEIN_INDEX_ITERATOR_T* new_protein_index_iterator(
  char* fasta_file ///< input fasta file -in
  )
{
  char* name = generate_name(fasta_file, "_protein_index", ".fasta", NULL);
  FILE* file = fopen(name, "r");

  if(file == NULL){
    carp(CARP_FATAL, "failed to open protein index file: %s", name);
  }

  free(name);

  PROTEIN_INDEX_ITERATOR_T* iterator = 
    (PROTEIN_INDEX_ITERATOR_T*)mycalloc(1, sizeof(PROTEIN_INDEX_ITERATOR_T));
  
  iterator->file = file;
  
  // set up the protein_index_iterator
  if(!setup_protein_index_iterator(iterator)){
    carp(CARP_FATAL, "failed to setup protein_index_iterator");
  }
  
  return iterator;
}

/**
 * Frees the allocated protein index iterator
 */
void free_protein_index_iterator(
  PROTEIN_INDEX_ITERATOR_T* protein_index_iterator ///< the iterator to free -in
  )
{
  // free the file handler
  fclose(protein_index_iterator->file);

  // free unused protein
  if(protein_index_iterator->next_protein != NULL){
    free_protein(protein_index_iterator->next_protein);
  }
  // free iterator
  free(protein_index_iterator);
}

/**
 *
 *\returns TRUE if there is another protein to return, else FALSE
 */
BOOLEAN_T protein_index_iterator_has_next(
  PROTEIN_INDEX_ITERATOR_T* protein_index_iterator ///< the iterator of interest -in
  )
{
  return protein_index_iterator->has_next;
}

/**
 *
 *\return the next protein in the protein index file
 */
PROTEIN_T* protein_index_iterator_next(
  PROTEIN_INDEX_ITERATOR_T* protein_index_iterator ///< the iterator of interest -in
  )
{
  PROTEIN_T* protein = protein_index_iterator->next_protein;
  protein_index_iterator->next_protein = NULL;

  // set up the protein_index_iterator
  if(!setup_protein_index_iterator(protein_index_iterator)){
    carp(CARP_FATAL, "failed to setup protein_index_iterator");
  }
  
  return protein;
}

/***********************************************
 *
 * Create binary fasta file
 *
 ***********************************************/

/**
 * \brief Writes serialized proteins from text fasta file to the
 * output_file stream.
 *
 * Expects open filestream and closes the file before returning.
 * \returns TRUE if successfully creates a binary fasta file, else false
 */
/*
 * TODO (BF 26-Feb-2008): I would prefer to see this either take two
 * file streams (input and output) and close neither or take two file
 * names (char*), open both files and close both files.  Is there a
 * reason it is being done this way?
 */
BOOLEAN_T create_binary_fasta_file(
  char* fasta_file,  ///< input fasta file -in
  FILE* output_file  ///< the output filestream -out
  )
{
  unsigned long working_index;
  FILE* file = NULL;
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  unsigned int protein_idx = 0;
  PROTEIN_T* new_protein = NULL;

  carp(CARP_DEBUG, "Creating binary fasta");
  // open file and 
  file = fopen(fasta_file, "r");

  // check if succesfully opened file
  if(file == NULL){
    carp(CARP_ERROR, "Failed to open fasta file '%s'", fasta_file);
    return FALSE;
  }

  // check if succesfully created file
  if(output_file == NULL){
    carp(CARP_FATAL, "Failed to create protein index file");
  }
  
  working_index = ftell(file);
  // check each line until reach '>' line
  while((line_length =  getline(&new_line, &buf_length, file)) != -1){
    if(new_line[0] == '>'){
      // the new protein to be serialize
      new_protein = allocate_protein();
      
      // rewind to the begining of the protein to include ">" line
      fseek(file, working_index, SEEK_SET);
          
      // failed to parse the protein from fasta file
      // protein offset is set in the parse_protein_fasta_file method
      if(!parse_protein_fasta_file(new_protein ,file)){
        fclose(file);
        free_protein(new_protein);
        carp(CARP_ERROR, "Failed to parse fasta file");
        return FALSE;
      }
      set_protein_is_light(new_protein, FALSE);
      
      // serialize protein as binary to output file
      serialize_protein(new_protein, output_file);

      // update protein count
      ++protein_idx;

      // free this protein
      free_protein(new_protein);
    }
    
    // print status
    if(protein_idx % 1000 == 0){
      carp(CARP_INFO, "Reached protein %d", protein_idx);
    }

    working_index = ftell(file);
  }

  // write the end character to binary fasta file
  char term_char = 1;  // use 1 and not '*' as the terminal
                       // character for the file b/c id length is 
                       // stored in same field and id len == 42
                       // is the smae as '*'
  fwrite(&term_char, sizeof(char), 1, output_file);

  // print final status
  //  carp(CARP_INFO, "Serialized total protein: %d", protein_idx);
  carp(CARP_INFO, "Total proteins found: %d", protein_idx);
  
    
  free(new_line);
  fclose(file);
  fclose(output_file);

  return TRUE;
}

/**
 * creates a binary fasta file in the output_file in the same directory as 
 * the fasta file
 * \returns TRUE if successfully creates a binary fasta file, else false
 */
BOOLEAN_T create_binary_fasta(
  char* fasta_file  ///< input fasta file -in
  )
{
  // get output file
  FILE* output_file = get_output_file(fasta_file, TRUE);
  
  return create_binary_fasta_file(fasta_file, output_file);
}

/**
 * wrapper for create_binary_fasta_file so that two filenames are
 * passed instead of a filename and a filestream.  Eventually should
 * merge to one method.  This should be putting binary fasta's in temp
 * dirs, so overwriting not an issue.  In the future, fix it (write a
 * create_file_in_path method with read/write info)
 */
BOOLEAN_T create_binary_fasta_here(
  char* fasta_filename,
  char* binary_filename
){

  FILE* output_file = fopen(binary_filename, "w");
  return create_binary_fasta_file(fasta_filename, output_file);
}

/**
 * creates a binary fasta file on to the output_file in currenty directory
 * sets the output file name to the pointer passed in as argument
 * \returns TRUE if successfully creates a binary fasta file, else false
 */

BOOLEAN_T create_binary_fasta_in_cur(
  char* fasta_file_w_path, ///< input fasta file with full path -in
  char* fasta_filename, ///< input fasta a file, only filename -in
  char** output_file_name ///< get output filename -out
  )
{
  // get output filename
  *output_file_name = generate_name(fasta_filename, "_binary_fasta", ".fasta", NULL);
  
  // open output file
  FILE* file = fopen(*output_file_name, "w");
  
  return create_binary_fasta_file(fasta_file_w_path, file);
}


/**
 * Heap allocated char*, user must free
 *\returns the binary fasta name which was created from the given fasta file
 */
char* get_binary_fasta_name(
  char* fasta_file  ///< input fasta file -in                            
  )
{
  // separate path from file name
  char** path_filename = parse_filename_path(fasta_file);

  // get binary fasta name
  char* binary_fasta_name = generate_name(path_filename[0], "_binary_fasta", ".fasta", NULL);

  // free path and filename
  free(path_filename[0]);
  if(path_filename[1] != NULL){
    free(path_filename[1]);
  }
  free(path_filename);
  
  carp(CARP_DETAILED_DEBUG, "binary fasta name: %s", binary_fasta_name);

  return binary_fasta_name;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

