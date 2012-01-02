/*************************************************************************//**
 * \file ProteinIndex.cpp
 * \brief Object for creating a protein index or binary fasta file
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"
#include "crux-utils.h"
#include "Peptide.h"
#include "Protein.h"
#include "Database.h"
#include "carp.h"
#include "objects.h"
#include "PeptideConstraint.h"
#include "ProteinIndex.h"
#include "WinCrux.h"


/***************/

/**
 *
 *\returns the file handler for the output file
 */
FILE* get_output_file(
  char* fasta_file,  ///< input fasta file -in
  bool is_binary_file  ///< Are we creating a binary fasta file?
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
  
  FILE* file = fopen(name, "wb");
  free(name);
  return file;
}


/**
 * creates a protein index on to the output_file
 * \returns true if successfully creates a protein index, else false
 */
bool ProteinIndex::create(
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
  file = fopen(fasta_file, "rb");

  // check if succesfully opened file
  if(file == NULL){
    carp(CARP_FATAL, "Failed to open fasta file '%s'", fasta_file);
    return false;
  }

  // get output file
  output_file = get_output_file(fasta_file, false);

  // check if succesfully created file
  if(output_file == NULL){
    carp(CARP_ERROR, "Failed to create protein index file");
    fclose(file);
    return false;
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

  return true;
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
 *\returns TRUE if protein index or binary fasta file is on disk, else false
 */
bool ProteinIndex::onDisk(
  const char* fasta_file, ///< input fasta file -in
  bool is_binary ///< are we looking for the binary fasta file? or preotein index
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
    return false;
  }
  
  free(name);
  return true;
}

/**
 * protein index iterator
 * the protein index iterator parses the protein index file
 * for each protein index one by one and returns a new light protein.
 */


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

