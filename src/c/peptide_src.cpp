/*************************************************************************//**
 * \file peptide_src.cpp
 * \brief Object for mapping a peptide to its parent protein.
 ****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "carp.h"
#include "utils.h"
#include "mass.h"
#include "objects.h"
#include "peptide.h"
#include "protein.h"
#include "peptide_src.h"
#include "peptide_constraint.h"

#include <vector>
#include <string>

#include "DelimitedFile.h"
#include "MatchFileReader.h"

using namespace std;

/**
 * \struct peptide_src
 * \brief object for mapping a peptide to it's parent protein.
 */
struct peptide_src{
  //PEPTIDE_TYPE_T peptide_type;///< the peptide type for the corresponding protein
  DIGEST_T digestion; ///< how specific the ends are relative to the enzyme
  PROTEIN_T* parent_protein; ///< the parent of this preptide
  int start_idx; ///< start index of the peptide in the protein sequence, first residue is 1 
  PEPTIDE_SRC_T* next_association; ///< a linklist of peptide_src     
};

/**
 * \returns An (empty) peptide_src object.
 */
PEPTIDE_SRC_T* allocate_peptide_src(void){
  PEPTIDE_SRC_T* peptide_src =
    (PEPTIDE_SRC_T*)mymalloc(sizeof(PEPTIDE_SRC_T));

  peptide_src->next_association = NULL;
  return peptide_src;
}

/**
 *\returns a PROTEIN_PEPTIDE_ASSOCIATION object, populated with user specified parameters
 */
PEPTIDE_SRC_T* new_peptide_src(
//PEPTIDE_TYPE_T peptide_type, ///< the peptide type for the corresponding protein -in
  DIGEST_T digest,
  PROTEIN_T* parent_protein, ///< the parent of this preptide -in
  int start_idx ///< start index of the peptide in the protein sequence -in
  )
{
  PEPTIDE_SRC_T* new_association = allocate_peptide_src();
  //set_peptide_src_peptide_type(new_association, peptide_type);
  set_peptide_src_digest(new_association, digest);
  set_peptide_src_parent_protein(new_association, parent_protein);
  set_peptide_src_start_idx(new_association, start_idx);
  return new_association;
}

/**
 *\returns an array of PROTEIN_PEPTIDE_SRC object
 * only used in index.c, when the peptide src count for  peptide is known
 */
PEPTIDE_SRC_T* new_peptide_src_array(
  int size ///< the size of the peptide_src array -in
  )
{
  int array_idx = 0;
  PEPTIDE_SRC_T* src_array = (PEPTIDE_SRC_T*)mycalloc(size,
                                                      sizeof(PEPTIDE_SRC_T));
  
  // set all next peptide src pointers
  for(array_idx = 0; array_idx < size - 1; ++array_idx){
    ((PEPTIDE_SRC_T*)(&(src_array[array_idx])))->next_association =
      &src_array[array_idx + 1];
  }
  ((PEPTIDE_SRC_T*)(&(src_array[array_idx])))->next_association = NULL;
  return src_array;
}

/**
 * \brief Fill in the values from the original array into the new array.
 * Assumes that the new array has been allocated by new_peptide_src_array().
 */
void copy_peptide_src_array(PEPTIDE_SRC_T* original_array, 
                            PEPTIDE_SRC_T* new_array, 
                            int array_size){
  int src_idx = 0;
  for(src_idx =0; src_idx < array_size; src_idx++){
    new_array[src_idx].digestion = original_array[src_idx].digestion;
    new_array[src_idx].parent_protein = original_array[src_idx].parent_protein;
    new_array[src_idx].start_idx = original_array[src_idx].start_idx;
  }

}

/**
 *\returns a linklist of PROTEIN_PEPTIDE_SRC object
 * only used in index.c, when the peptide src count for peptide is known
 */
PEPTIDE_SRC_T* new_peptide_src_linklist(
  int size ///< the size of the peptide_src array -in
  )
{
  int src_idx = 1;
  // create one peptide src
  PEPTIDE_SRC_T* src_list = (PEPTIDE_SRC_T*)mycalloc(1, sizeof(PEPTIDE_SRC_T));
  PEPTIDE_SRC_T* curr_src = src_list;

  // set all next peptide src pointers, if size is greater than 1
  for(src_idx = 0; src_idx < size - 1; ++src_idx){
    curr_src->next_association = (PEPTIDE_SRC_T*)mycalloc(1, sizeof(PEPTIDE_SRC_T));
    curr_src = curr_src->next_association;
  }

  // set last association to null
  curr_src->next_association = NULL;
  return src_list;
}

/**
 *\returns the PROTEIN_PEPTIDE_SRC object in the array with the index
 * index starts at 0.
 * only used in index.c, when the peptide src count for  peptide is known
 */
void set_peptide_src_array(
  PEPTIDE_SRC_T* src_array , ///< the working peptide src_arry -out
  int array_idx, ///< array index of the peptide_src to set
  //PEPTIDE_TYPE_T peptide_type, ///< the peptide type for the corresponding protein -in
  DIGEST_T digest,
  PROTEIN_T* parent_protein, ///< the parent of this preptide -in
  int start_idx ///< start index of the peptide in the protein sequence -in
  )
{
  // set all valuse
  PEPTIDE_SRC_T* peptide_src = &src_array[array_idx];
  //  peptide_src->peptide_type = peptide_type;
  peptide_src->digestion = digest;
  peptide_src->parent_protein = parent_protein;
  peptide_src->start_idx = start_idx;
}


/**
 * Frees the entire allocated peptide_src linklist object
 * Assumes that peptide src is Link list implementation
 */
void free_peptide_src(
  PEPTIDE_SRC_T* peptide_src  ///< object to free -in 
  )
{
  PEPTIDE_SRC_T* to_free = peptide_src;
  PEPTIDE_SRC_T* next = peptide_src->next_association;

  // free first peptide_src
  free(to_free);
  
  // iterate over all peptide_srcs
  while(next != NULL){
    to_free = next;
    next = next->next_association;
    free(to_free);
  }

  /*
  if(peptide_src->next_association != NULL){
    free_peptide_src(peptide_src->next_association);
  }  
  free(peptide_src);  
  */
}

/**
 * Frees the an individual allocated peptide_src object
 * assumes that new_association pointer is NULL or some other pointer exist for the rest of the linklist 
 */
void free_one_peptide_src(
  PEPTIDE_SRC_T* peptide_src  ///< object to free -in 
  )
{
  free(peptide_src);
}

// FIXME might need to change how this is printed
/**
 * Prints a peptide object to file.
 */
/*
void print_peptide_src(
  PEPTIDE_SRC_T* peptide_src, ///< object to print -in 
  FILE* file  ///< the out put stream -out
  )
{
  PEPTIDE_TYPE_T peptide_type = get_peptide_src_peptide_type(peptide_src);
  char* sequence = get_protein_sequence(peptide_src->parent_protein);
  fprintf(file, "parent protein:%s\n", sequence);

  fprintf(file, "peptide start: %d\n", peptide_src->start_idx);

  if(peptide_type == TRYPTIC){
    fprintf(file, "peptide type:%s\n", "TRYPTIC");
  }
  else if(peptide_type == PARTIALLY_TRYPTIC){
    fprintf(file, "peptide type:%s\n", "PARTIALLY_TRYPTIC");
  }
  else if(peptide_type == N_TRYPTIC){
    fprintf(file, "%s", "N_TRYPTIC");
  }
  else if(peptide_type == C_TRYPTIC){
    fprintf(file, "%s", "C_TRYPTIC");
  }
  else if(peptide_type == NOT_TRYPTIC){
    fprintf(file, "peptide type:%s\n", "NOT_TRYPTIC");
  }
  else if(peptide_type == ANY_TRYPTIC){
    fprintf(file, "peptide type:%s\n", "ANY_TRYPTIC");
  }
  free(sequence);
}
*/

/**
 * Copies the entire linklist of peptide_src object src to dest.
 * dest must be a heap allocated peptide_src
 */
void copy_peptide_src(
  PEPTIDE_SRC_T* src, ///< source peptide_src -in
  PEPTIDE_SRC_T* dest ///< destination peptide_src -out
  )
{
  PEPTIDE_SRC_T* next_association;
  //set_peptide_src_peptide_type(dest, src->peptide_type);
  set_peptide_src_digest(dest, src->digestion);
  set_peptide_src_parent_protein(dest, src->parent_protein);
  set_peptide_src_start_idx(dest, src->start_idx);
  // check if end of the linklist
  if(get_peptide_src_next_association(src) != NULL){
    next_association = allocate_peptide_src();
    dest->next_association = next_association;
    copy_peptide_src(src->next_association, next_association);
  }
}

/**
 * sets the peptide type
 * peptide type: TRYPTIC, PARTIALLY_TRYPTIC, N_TRYPTIC, C_TRYPTIC, NON_TRYPTIC
 */
/*
void set_peptide_src_peptide_type( 
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  PEPTIDE_TYPE_T peptide_type ///< the type of the peptide -in
  )
{
  new_association->peptide_type = peptide_type;
}
*/
/**
 * \returns the peptide type with association to the parent protein
 * peptide type: TRYPTIC, PARTIALLY_TRYPTIC, N_TRYPTIC, C_TRYPTIC, NON_TRYPTIC
 */
 /*
PEPTIDE_TYPE_T get_peptide_src_peptide_type( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  )
{
  return peptide_src->peptide_type;
}
 */

/**
 * sets the level of digestion
 */
void set_peptide_src_digest( 
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  DIGEST_T digest ///< the type of the peptide -in
  ){
  new_association->digestion = digest;
}

/**
 * \returns the level of digestion
 */
DIGEST_T get_peptide_src_digest( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  ){
  return peptide_src->digestion;
}

/**
 * sets the parent protein
 */
void set_peptide_src_parent_protein(
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  PROTEIN_T* parent_protein ///< the parent of this preptide -in  
  )
{
  new_association->parent_protein = parent_protein;

}

/**
 * \returns a pointer to the parent protein
 */
PROTEIN_T* get_peptide_src_parent_protein( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  )
{
  return peptide_src->parent_protein;
}

/**
 * sets the start index of the peptide in the protein sequence
 */
void set_peptide_src_start_idx(
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  int start_idx ///< start index of the peptide in the protein sequence -in
  )
{
  new_association->start_idx = start_idx;
}

/**
 * \returns the start index of the peptide in the protein sequence
 */
int get_peptide_src_start_idx( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  )
{
  return peptide_src->start_idx;
}

/**
 * sets the next peptide_src on the link list
 * assumes that the src_association's next_association field is NULL
 */
void set_peptide_src_next_association(
  PEPTIDE_SRC_T* src_association, ///< the peptide_src to set -out   
  PEPTIDE_SRC_T* new_association ///< the new peptide_src to add -in   
  )
{
  src_association->next_association = new_association;
}

/**
 * \returns the next peptide_src on the link list
 */
PEPTIDE_SRC_T* get_peptide_src_next_association( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  )
{
  return peptide_src->next_association;
}

/**
 * \returns a pointer to the start of the peptide with in it's parent protein sequence
 */
char* get_peptide_src_sequence_pointer(
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  )
{
  char* start_pointer = get_protein_sequence_pointer(peptide_src->parent_protein);
  return &(start_pointer[peptide_src->start_idx - 1]);
}

/**
 *\returns the peptide_src strct size, value of sizeof function
 */
int get_peptide_src_sizeof(){
  return sizeof(PEPTIDE_SRC_T);
}

/**
 * serialize peptide src in binary
 * The peptide serialization format looks like this:
 *
 *<int: protein index><PEPTIDE_TYPE_T: peptide_type><int: peptide start index>
 * the protein index is the index of the parent protein in the
 *database DATABASE_T 
 *
 */
void serialize_peptide_src(
  PEPTIDE_SRC_T* peptide_src, ///< peptide_src to serialize -in   
  FILE* file  ///< output file -in   
  )
{
  // write protein index in database
  unsigned int protein_idx = 
    get_protein_protein_idx(peptide_src->parent_protein);
  // carp(CARP_DETAILED_DEBUG, "protein idx to write is %i", protein_idx);

  fwrite(&protein_idx, sizeof(int), 1, file);
  carp(CARP_DETAILED_DEBUG, "Serializing protein src of index %d", 
       protein_idx); 
   
  // write peptide src type(tryptic, all, ...)
  //  fwrite(&(peptide_src->peptide_type), sizeof(PEPTIDE_TYPE_T), 1, file);
  fwrite(&(peptide_src->digestion), sizeof(DIGEST_T), 1, file);
  // write start index in protein of peptide in this peptide src
  fwrite(&(peptide_src->start_idx), sizeof(int), 1, file);
  
}

/**
 * Return the number of bytes taken up by one peptide_src when
 * serialized to file.  Used for skipping past peptide_src in an index
 * file. 
 */
int size_of_serialized_peptide_src(){
  //return (sizeof(int)*2 + sizeof(PEPTIDE_TYPE_T));
  return (sizeof(int)*2 + sizeof(DIGEST_T));
}

/**
 * \brief Read in the peptide_src objects from the given file and
 * assosiated them with the given peptide.  
 * Proteins for the pepitde_src are found in the given database.  If
 * database is NULL, does not set proteins.  (This option is used for
 * sorting index files while creating index.)  Either array or 
 * linked list implementation of multiple peptide_src is used based on
 * the value of use_array.
 *
 * \returns TRUE if peptide_src's were successfully parsed, else
 * returns FALSE.
 */
BOOLEAN_T parse_peptide_src_tab_delimited(
  PEPTIDE_T* peptide,   ///< assign peptide_src(s) to this peptide
  MatchFileReader& file,           ///< file to read from
  DATABASE_T* database, ///< database containing proteins
  BOOLEAN_T use_array) ///< use array implementation vs. linked list
{
  //TODO - Implement

  if( peptide == NULL ){
    carp(CARP_ERROR, "Cannot parse peptide src with NULL peptide.");
    return FALSE;
  }

  carp(CARP_DETAILED_DEBUG,"Parsing id line:%s", file.getString(PROTEIN_ID_COL).c_str());

  vector<string> protein_ids;
  file.getStringVectorFromCell(PROTEIN_ID_COL, protein_ids);

  int num_peptide_src = protein_ids.size();
  
  PEPTIDE_SRC_T* peptide_src = NULL;

  // allocate new src based on requested type
  if(use_array){
    peptide_src = new_peptide_src_array(num_peptide_src);
  } else {
    peptide_src = new_peptide_src_linklist(num_peptide_src);
  }

  // give it to the peptide
  add_peptide_peptide_src_array(peptide, peptide_src);

  DIGEST_T digestion = 
    string_to_digest_type((char*)file.getString(CLEAVAGE_TYPE_COL).c_str()); 
  
  PROTEIN_T* parent_protein = NULL;
  int start_index = 1;
  for (vector<string>::iterator iter = protein_ids.begin();
    iter != protein_ids.end();
    ++iter) {

    carp(CARP_DETAILED_DEBUG,"Parsing %s",iter -> c_str());
    // get the protein and peptide index e.g. X(10)
    size_t left_paren_index = iter -> find('(');

    if (left_paren_index == string::npos) {
      //protein id is the string.
      string protein_id_string = *iter;
      
      parent_protein =
        get_database_protein_by_id_string(database, protein_id_string.c_str());
      
      //find the start index by searching the protein sequence.
      string protein_sequence(get_protein_sequence_pointer(parent_protein));

      //if sequence is decoy sequence, recover the position from
      //the unshuffled sequence.
      string sequence;
      if (file.getString(UNSHUFFLED_SEQUENCE_COL).empty()) {
        sequence = file.getString(SEQUENCE_COL);
      } else {
        sequence = file.getString(UNSHUFFLED_SEQUENCE_COL);
      }
      size_t pos = protein_sequence.find(sequence);

      if (pos == string::npos) {
        carp(CARP_FATAL, "Can't find sequence %s in %s",sequence.c_str(), protein_sequence.c_str());
      }
      start_index = (int)pos + 1;

    } else {
      string protein_id_string = iter -> substr(0, left_paren_index);
      string peptide_start_index_string = iter -> substr(left_paren_index+1, 
      iter -> length() - 1);

      //  set fields in new peptide src
      parent_protein =
        get_database_protein_by_id_string(database, protein_id_string.c_str());
     
      if (parent_protein == NULL) {
        carp(CARP_FATAL, "Can't find protein %s", iter -> c_str());
        continue;
      }

      DelimitedFile::from_string<int>(start_index, peptide_start_index_string); 
    }
    // set parent protein of the peptide src
    set_peptide_src_parent_protein(peptide_src, parent_protein);

    // set digest type of peptide src
    set_peptide_src_digest(peptide_src, digestion);

    // set start index of peptide src
    set_peptide_src_start_idx(peptide_src, start_index);

    // set current peptide_src to the next empty peptide src
    peptide_src = get_peptide_src_next_association(peptide_src);
  } // next peptide_src in file

  return TRUE;
}


/**
 * \brief Read in the peptide_src objects from the given file and
 * assosiated them with the given peptide.  
 * Proteins for the pepitde_src are found in the given database.  If
 * database is NULL, does not set proteins.  (This option is used for
 * sorting index files while creating index.)  Either array or 
 * linked list implementation of multiple peptide_src is used based on
 * the value of use_array.
 *
 * \returns TRUE if peptide_src's were successfully parsed, else
 * returns FALSE.
 */
BOOLEAN_T parse_peptide_src(
  PEPTIDE_T* peptide,   ///< assign peptide_src(s) to this peptide
  FILE* file,           ///< file to read from
  DATABASE_T* database, ///< database containing proteins
  BOOLEAN_T use_array) ///< use array implementation vs. linked list
{
  if( peptide == NULL || file == NULL ){
    carp(CARP_ERROR, "Cannot parse peptide src with NULL peptide or file.");
    return FALSE;
  }

  // first field in file should be number of src's
  int num_peptide_src = -1;
  size_t num_read = fread(&num_peptide_src, sizeof(int), 1, file);
  if( num_peptide_src < 1 || num_read != 1 ){
    carp(CARP_ERROR, 
         "Index file corrupted, peptide must have at least one peptide src");
    return FALSE;
  }

  PEPTIDE_SRC_T* peptide_src = NULL;
  // allocate new src based on requested type
  if(use_array){
    peptide_src = new_peptide_src_array(num_peptide_src);
  }else{
    peptide_src = new_peptide_src_linklist(num_peptide_src);
  }

  // give it to the peptide
  add_peptide_peptide_src_array(peptide, peptide_src);

  // read in each peptide_src (prot index, peptide type, start index)
  int src_idx = 0;
  int protein_index = -1;
  PROTEIN_T* parent_protein = NULL;
  //PEPTIDE_TYPE_T peptide_type;
  DIGEST_T digestion;
  int start_index = -1;
  for(src_idx = 0; src_idx < num_peptide_src; src_idx++){
    // get protein index
    size_t num_read = fread(&protein_index, (sizeof(int)), 1, file);
    if( protein_index < 0 || num_read != 1){
      carp(CARP_ERROR, "Index file corrupted could not read protein index");
      free(peptide_src);
      return FALSE;
    }
    carp(CARP_DETAILED_DEBUG, "Peptide src %d has protein idx %i", 
         src_idx, protein_index);

    // read peptide type of peptide src
    /*
    if(fread(&peptide_type, sizeof(PEPTIDE_TYPE_T), 1, file) != 1){
      carp(CARP_ERROR, "Index file corrupted, failed to read peptide type.");
      free(peptide_src);
      return FALSE;
    }
    */

    // read digestion level of peptide src
    if(fread(&digestion, sizeof(DIGEST_T), 1, file) != 1){
      carp(CARP_ERROR, "Index file corrupted, failed to read digestion type.");
      free(peptide_src);
      return FALSE;
    }

    // read start index of peptide in parent protein of thsi peptide src
    if(fread(&start_index, sizeof(int), 1, file) != 1){
      carp(CARP_ERROR, "Index file corrupted, failed to read start index.");
      free(peptide_src);
      return FALSE;
    }

    // set fields in new peptide src
    parent_protein = 
      get_database_protein_at_idx(database, protein_index);
    
    // set parent protein of the peptide src
    set_peptide_src_parent_protein(peptide_src, parent_protein);

    // set peptide type of peptide src
    //    set_peptide_src_peptide_type(peptide_src, peptide_type);
    set_peptide_src_digest(peptide_src, digestion);

    // set start index of peptide src
    set_peptide_src_start_idx(peptide_src, start_index);
    
    // set current_peptide_src to the next empty peptide src
    peptide_src = get_peptide_src_next_association(peptide_src);
  }// next peptide_src in file

  carp(CARP_DETAILED_DEBUG, "Finished parsing peptide src.");
  return TRUE;
}



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
