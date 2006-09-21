/*****************************************************************************
 * \file protein_index.c
 * $Revision: 1.1 $
 * \brief: Object for creating a protein index
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
  char* fasta_file  ///< input fasta file -in
  )
{
  char* name = generate_name(fasta_file, "_protein_index");
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

  //open file and 
  file = fopen(fasta_file, "r");

  //check if succesfully opened file
  if(file == NULL){
    carp(CARP_FATAL, "failed to open fasta file");
    return FALSE;
  }

  //get output file
  output_file = get_output_file(fasta_file);

  //check if succesfully created file
  if(output_file == NULL){
    carp(CARP_FATAL, "failed to create protein index file");
    fclose(file);
    return FALSE;
  }

  working_index = ftell(file);
  //check each line until reach '>' line
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
 *\returns TRUE if protein index is on disk, else FALSE
 */
BOOLEAN_T protein_index_on_disk(
  char* fasta_file
  )
{
  char* name = generate_name(fasta_file, "_protein_index");

  //check if can open file
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
  //used to parse each line from file
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  FILE* file = protein_index_iterator->file;

  //protein fields
  char star[2] = "";
  unsigned long int offset;
  unsigned int protein_idx;
  
  PROTEIN_T* protein = NULL;
  BOOLEAN_T found = FALSE;

  while((line_length =  getline(&new_line, &buf_length, file)) != -1){
    //begining of the protein feilds
    if(new_line[0] == '*'){
      //read the crux_index_file information
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
    //skip header lines
    else if(new_line[0] == '#'){
      continue;
    }
  }
  
  //there is a next protein to return
  if(found){
    protein = new_light_protein(offset, protein_idx);
    protein_index_iterator->next_protein = protein;
    protein_index_iterator->has_next = TRUE;
  }
  //no more proteins..
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
  char* name = generate_name(fasta_file, "_protein_index");
  FILE* file = fopen(name, "r");

  if(file == NULL){
    carp(CARP_FATAL, "failed to open protein index file: %s", name);
    exit(1);
  }

  free(name);

  PROTEIN_INDEX_ITERATOR_T* iterator = 
    (PROTEIN_INDEX_ITERATOR_T*)mycalloc(1, sizeof(PROTEIN_INDEX_ITERATOR_T));
  
  iterator->file = file;
  
  //set up the protein_index_iterator
  if(!setup_protein_index_iterator(iterator)){
    carp(CARP_FATAL, "failed to setup protein_index_iterator");
    exit(1);
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
  //free the file handler
  fclose(protein_index_iterator->file);

  //free unused protein
  if(protein_index_iterator->next_protein != NULL){
    free_protein(protein_index_iterator->next_protein);
  }
  //free iterator
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

  //set up the protein_index_iterator
  if(!setup_protein_index_iterator(protein_index_iterator)){
    free_protein_index_iterator(protein_index_iterator);
    carp(CARP_FATAL, "failed to setup protein_index_iterator");
    exit(1);
  }
  
  return protein;
}
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

