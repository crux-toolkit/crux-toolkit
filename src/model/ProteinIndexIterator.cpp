/*************************************************************************//**
 * \file ProteinIndexIterator.cpp
 * \brief Object for iterating a protein index file
 ****************************************************************************/
#include "ProteinIndexIterator.h"

using namespace Crux;

/**
 *
 *\returns TRUE if successfully sets the protein_index_iterator, else FALSE
 */
bool ProteinIndexIterator::setup() {

  // used to parse each line from file
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  FILE* file = file_;

  // protein fields
  char star[2] = "";
  unsigned long int offset;
  unsigned int protein_idx;
  
  Protein* protein = NULL;
  bool found = false;

  while((line_length =  getline(&new_line, &buf_length, file)) != -1){
    // begining of the protein feilds
    if(new_line[0] == '*'){
      // read the crux_index_file information
      if(sscanf(new_line,"%s %d %ld", 
                star, &protein_idx, &offset) < 2){
        free(new_line);
        carp(CARP_WARNING, "incorrect file format");
        fclose(file);
        return false;
      }
      found = true;
      break;
    }
    // skip header lines
    else if(new_line[0] == '#'){
      continue;
    }
  }
  
  // there is a next protein to return
  if(found){
    protein = Protein::newLightProtein(offset, protein_idx);
    next_protein_ = protein;
    has_next_ = true;
  }
  // no more proteins..
  else{
    has_next_ = false;
  }
  
  free(new_line);
  return true;
}

/**
 * There must be the correct protein index file present
 *\returns a new heap allocated protein index iterator
 */
ProteinIndexIterator::ProteinIndexIterator(
  const char* fasta_file ///< input fasta file -in
  )
{
  char* name = generate_name(fasta_file, "_protein_index", ".fasta", NULL);
  FILE* file = fopen(name, "rb");

  if(file == NULL){
    carp(CARP_FATAL, "failed to open protein index file: %s", name);
  }

  free(name);
  
  next_protein_ = NULL;
  has_next_ = false;
  file_ = file;
  
  // set up the protein_index_iterator
  if(!setup()){
    carp(CARP_FATAL, "failed to setup ProteinIndexIterator");
  }
}

/**
 * Frees the allocated protein index iterator
 */
ProteinIndexIterator::~ProteinIndexIterator() {

  // free the file handler
  fclose(file_);

  // free unused protein
  if(next_protein_ != NULL){
    delete (next_protein_);
  }
}

/**
 *
 *\returns TRUE if there is another protein to return, else FALSE
 */
bool ProteinIndexIterator::hasNext() {

  return has_next_;
}

/**
 *\return the next protein in the protein index file
 */
Protein* ProteinIndexIterator::next() {

  Protein* protein = next_protein_;
  next_protein_ = NULL;

  // set up the protein_index_iterator
  if(!setup()){
    carp(CARP_FATAL, "failed to setup protein_index_iterator");
  }
  
  return protein;
}
