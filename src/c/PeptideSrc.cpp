/*************************************************************************//**
 * \file PeptideSrc.cpp
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
#include "DatabasePeptideIterator.h"
#include "Peptide.h"
#include "Protein.h"
#include "PeptideSrc.h"
#include "PeptideConstraint.h"
#include "PeptideIterator.h"

#include <vector>
#include <string>

#include "DelimitedFile.h"
#include "MatchFileReader.h"

using namespace std;


/**
 * Static variable definitions
 */
map<string, Peptide* > PeptideSrc::sequence_to_peptide_; ///< Maps a sequence to a peptide object
map<string, Peptide* > PeptideSrc::decoy_sequence_to_peptide_; ///< Maps a decoy sequence to a peptide object


/**
 * Initializes an (empty) PeptideSrc object
 */
void PeptideSrc::init() {

  digestion_ = (DIGEST_T)0;
  parent_protein_ = NULL;
  start_idx_ = 0;
  next_association_ = NULL;
}

/**
 * \returns An (empty) peptide_src object.
 */
PeptideSrc::PeptideSrc() {

  init();

}

/**
 *\returns a PROTEIN_PEPTIDE_ASSOCIATION object, populated with user specified parameters
 */
PeptideSrc::PeptideSrc(
  DIGEST_T digest,
  Protein* parent_protein, ///< the parent of this preptide -in
  int start_idx ///< start index of the peptide in the protein sequence -in
  ) {

  init();
  setDigest(digest);
  setParentProtein(parent_protein);
  setStartIdx(start_idx);
}

/**
 *\returns an array of PROTEIN_PEPTIDE_SRC object
 * only used in index.c, when the peptide src count for  peptide is known
 */
PeptideSrc* PeptideSrc::newArray(
  int size ///< the size of the peptide_src array -in
  ) {

  PeptideSrc* src_array = new PeptideSrc[size];

  src_array[0].init();
  // set all next peptide src pointers
  for(int array_idx = 1; array_idx < size ; array_idx++){
    src_array[array_idx].init();
    src_array[array_idx-1].setNextAssociation(&(src_array[array_idx]));
  }
  return src_array;
}

/**
 * \brief Fill in the values from the original array into the new array.
 * Assumes that the new array has been allocated by new_peptide_src_array().
 */
void PeptideSrc::copyArray(
  PeptideSrc* original_array, 
  PeptideSrc* new_array, 
  int array_size){

  for(int src_idx = 0; src_idx < array_size; src_idx++){
    new_array[src_idx].digestion_ = original_array[src_idx].digestion_;
    new_array[src_idx].parent_protein_ = original_array[src_idx].parent_protein_;
    new_array[src_idx].start_idx_ = original_array[src_idx].start_idx_;
  }

}

/**
 *\returns a linklist of PROTEIN_PEPTIDE_SRC object
 * only used in index.c, when the peptide src count for peptide is known
 */
PeptideSrc* PeptideSrc::newLinklist(
  int size ///< the size of the peptide_src array -in
  ) {

  // create one peptide src
  PeptideSrc* src_list = new PeptideSrc();
  PeptideSrc* curr_src = src_list;

  // set all next peptide src pointers, if size is greater than 1
  for(int src_idx = 0; src_idx < size - 1; ++src_idx){
    curr_src->next_association_ = new PeptideSrc();
    curr_src = curr_src->next_association_;
  }

  // set last association to null
  curr_src->next_association_ = NULL;
  return src_list;
}

/**
 *\returns the PROTEIN_PEPTIDE_SRC object in the array with the index
 * index starts at 0.
 * only used in index.c, when the peptide src count for  peptide is known
 */
void PeptideSrc::setArray(
  PeptideSrc* src_array , ///< the working peptide src_arry -out
  int array_idx, ///< array index of the peptide_src to set
  DIGEST_T digest,
  Protein* parent_protein, ///< the parent of this preptide -in
  int start_idx ///< start index of the peptide in the protein sequence -in
  ) {

  // set all valuse
  PeptideSrc* peptide_src = &src_array[array_idx];
  peptide_src->digestion_ = digest;
  peptide_src->parent_protein_ = parent_protein;
  peptide_src->start_idx_ = start_idx;
}

/**
 * Frees the entire allocated peptide_src linklist object
 * Assumes that peptide src is Link list implementation
 */
void PeptideSrc::free(PeptideSrc* peptide_src) {

  PeptideSrc* to_free = peptide_src;
  PeptideSrc* next = peptide_src->next_association_;

  // free first peptide_src
  delete to_free;
  
  // iterate over all peptide_srcs
  while(next != NULL){
    to_free = next;
    next = next->next_association_;
    delete to_free;
  }
}

/**
 * Frees the an individual allocated peptide_src object
 * assumes that new_association pointer is NULL or some other pointer exist for the rest of the linklist 
 */
PeptideSrc::~PeptideSrc() {
  ;
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
void PeptideSrc::copy(
  PeptideSrc* src, ///< source peptide_src -in
  PeptideSrc* dest ///< destination peptide_src -out
  )
{
  PeptideSrc* next_association;
  //set_peptide_src_peptide_type(dest, src->peptide_type);
  dest->setDigest(src->digestion_);
  dest->setParentProtein(src->parent_protein_);
  dest->setStartIdx(src->start_idx_);
  // check if end of the linklist
  if(src->getNextAssociation() != NULL){
    next_association = new PeptideSrc();
    dest->next_association_ = next_association;
    copy(src->next_association_, next_association);
  }
}

/**
 * sets the level of digestion
 */
void PeptideSrc::setDigest(
  DIGEST_T digest ///< the type of the peptide -in
  ){

  digestion_ = digest;
}

/**
 * \returns the level of digestion
 */
DIGEST_T PeptideSrc::getDigest() {

  return digestion_;
}

/**
 * sets the parent protein
 */
void PeptideSrc::setParentProtein(
  Protein* parent_protein ///< the parent of this preptide -in  
  ) {

  parent_protein_ = parent_protein;
}

/**
 * \returns a pointer to the parent protein
 */
Protein* PeptideSrc::getParentProtein() {

  return parent_protein_;
}

/**
 * sets the start index of the peptide in the protein sequence
 */
void PeptideSrc::setStartIdx(
  int start_idx ///< start index of the peptide in the protein sequence -in
  ) {

  start_idx_ = start_idx;
}

/**
 * \returns the start index of the peptide in the protein sequence
 */
int PeptideSrc::getStartIdx() {

  return start_idx_;
}

/**
 * sets the next peptide_src on the link list
 * assumes that the src_association's next_association field is NULL
 */
void PeptideSrc::setNextAssociation(
  PeptideSrc* new_association ///< the new peptide_src to add -in   
  ) {

  next_association_ = new_association;
}

/**
 * \returns the next peptide_src on the link list
 */
PeptideSrc* PeptideSrc::getNextAssociation() {

  return next_association_;
}

/**
 * \returns a pointer to the start of the peptide with in it's parent protein sequence
 */
char* PeptideSrc::getSequencePointer() {

  char* start_pointer = parent_protein_->getSequencePointer();
  return &(start_pointer[start_idx_ - 1]);
}

/**
 *\returns the peptide_src strct size, value of sizeof function
 */
int PeptideSrc::getSizeOf(){
  return sizeof(PeptideSrc);
}

/**
 * serialize peptide src in binary
 * The peptide serialization format looks like this:
 *
 *<int: protein index><PEPTIDE_TYPE_T: peptide_type><int: peptide start index>
 * the protein index is the index of the parent protein in the
 *database Database 
 *
 */
void PeptideSrc::serialize(
  FILE* file  ///< output file -in   
  ) {

  // write protein index in database
  unsigned int protein_idx = parent_protein_->getProteinIdx();
  // carp(CARP_DETAILED_DEBUG, "protein idx to write is %i", protein_idx);

  fwrite(&protein_idx, sizeof(int), 1, file);
  carp(CARP_DETAILED_DEBUG, "Serializing protein src of index %d", 
       protein_idx); 
   
  // write peptide src type(tryptic, all, ...)
  fwrite(&(digestion_), sizeof(DIGEST_T), 1, file);
  // write start index in protein of peptide in this peptide src
  fwrite(&(start_idx_), sizeof(int), 1, file);
  
}

/**
 * Return the number of bytes taken up by one peptide_src when
 * serialized to file.  Used for skipping past peptide_src in an index
 * file. 
 */
int PeptideSrc::sizeOfSerialized(){

  return (sizeof(int)*2 + sizeof(DIGEST_T));
}

/**
 * fills the sequence_to_peptide_ member variable for use in parseTabDelimited
 * used when the tab delimited file doesn't provide a protein id, but we have
 * sequences and access to the database.
 */
void PeptideSrc::fillPeptides(
  Database* database, ///< the protein database 
  Database* decoy_database ///< the decoy database
  ) {
  
  PeptideConstraint* constraint = PeptideConstraint::newFromParameters();
  PeptideIterator* iterator = new DatabasePeptideIterator(database, constraint, true, false);

  while (iterator->hasNext()) {

    Peptide* peptide = iterator->next();
    char* sequence = peptide->getSequence();
    string sequence_string = string(sequence);
    std::free(sequence);

    if (sequence_to_peptide_.find(sequence_string) == sequence_to_peptide_.end()) {
      sequence_to_peptide_[sequence_string] = peptide;
    }
  }
  delete iterator;
  delete constraint;


  //TODO - map the decoy sequences to Peptides if we need them in the future.
  if (decoy_database == NULL) {
    carp(CARP_INFO, "decoy database is null");
  }

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
 * \returns true if peptide_src's were successfully parsed, else
 * returns false.
 */
bool PeptideSrc::parseTabDelimited(
  Peptide* peptide,   ///< assign peptide_src(s) to this peptide
  MatchFileReader& file,           ///< file to read from
  Database* database, ///< database containing proteins
  bool use_array, ///< use array implementation vs. linked list
  Database* decoy_database ///< database containing decoy proteins
) {

  if( peptide == NULL ){
    carp(CARP_ERROR, "Cannot parse peptide src with NULL peptide.");
    return false;
  }

  carp(CARP_DETAILED_DEBUG,"Parsing id line:%s", 
       file.getString(PROTEIN_ID_COL).c_str());

  //if the protein id field is empty, then we have to search the database...
  if (file.empty(PROTEIN_ID_COL)) {
    carp_once(CARP_WARNING, "empty protein id string in tab delimited file. "
                            "searching database to find proteins to match peptide "
                            "sequence");

    //if we haven't done this already, build a map of sequence strings to peptide
    //objects.
    if (sequence_to_peptide_.size() == 0) {
      fillPeptides(database, decoy_database);
    }

    char* seq = peptide->getUnshuffledModifiedSequence();
    string seq_string(seq);
    std::free(seq);
    
    if (sequence_to_peptide_.find(seq_string) == sequence_to_peptide_.end()) {
      carp(CARP_WARNING, "Cannot find peptide in database!");
      return false;
    }
    Peptide* src_peptide = sequence_to_peptide_[seq_string];

    PeptideSrc* dest = NULL;

    carp_once(CARP_INFO, 
      "adding %d proteins to peptide %s",
      src_peptide->getNumPeptideSrc(), seq_string.c_str());
    if (use_array) {
      dest = newArray(src_peptide->getNumPeptideSrc());
      PeptideSrc::copy(src_peptide->getPeptideSrc(), dest);
    }  else {
      dest = newLinklist(src_peptide->getNumPeptideSrc());
      PeptideSrc::copy(src_peptide->getPeptideSrc(), dest);
      
    }
    peptide->addPeptideSrcArray(dest);
    peptide->getParentProtein()->getProteinIdx();
    return true;


  } else {
    vector<string> protein_ids;
    file.getStringVectorFromCell(PROTEIN_ID_COL, protein_ids);
  
    if (protein_ids.size() == 0) {
      carp(CARP_ERROR, "No protein ids found!");
      return false;
    }

  int num_peptide_src = protein_ids.size();
  
  PeptideSrc* peptide_src = NULL;

  // allocate new src based on requested type
  if(use_array){
    peptide_src = newArray(num_peptide_src);
  } else {
    peptide_src = newLinklist(num_peptide_src);
  }

  // give it to the peptide
  peptide->addPeptideSrcArray(peptide_src);

  DIGEST_T digestion = 
    string_to_digest_type((char*)file.getString(CLEAVAGE_TYPE_COL).c_str()); 
  
  Protein* parent_protein = NULL;
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
        database->getProteinByIdString(protein_id_string.c_str());
      
      if (parent_protein == NULL) {
       if( decoy_database != NULL ){
          parent_protein =
            decoy_database->getProteinByIdString(protein_id_string.c_str());
        }
        if (parent_protein == NULL) {
          carp(CARP_WARNING, "Can't find protein %s",protein_id_string.c_str());
          continue;
        }
      }

      //find the start index by searching the protein sequence.
      string protein_sequence(parent_protein->getSequencePointer());

      //if sequence is decoy sequence, recover the position from
      //the unshuffled sequence.
      string sequence;

      if (file.empty(UNSHUFFLED_SEQUENCE_COL)) {
        sequence = file.getString(SEQUENCE_COL);
      } else {
        sequence = file.getString(UNSHUFFLED_SEQUENCE_COL);
      }
      size_t pos = protein_sequence.find(sequence);

      if (pos == string::npos) {
        carp(CARP_WARNING, "Can't find sequence %s in %s:%s",
          sequence.c_str(),
          protein_id_string.c_str(),
          protein_sequence.c_str());
        pos = 0;
      }
      start_index = (int)pos + 1;

    } else {
      string protein_id_string = iter -> substr(0, left_paren_index);
      string peptide_start_index_string = iter -> substr(left_paren_index+1, 
      iter -> length() - 1);

      //  set fields in new peptide src
      parent_protein =
        database->getProteinByIdString(protein_id_string.c_str());
     
      if (parent_protein == NULL) {
       if( decoy_database != NULL ){
          parent_protein =
            decoy_database->getProteinByIdString(protein_id_string.c_str());
        }
        if (parent_protein == NULL) {
          carp(CARP_WARNING, "Can't find protein %s", iter -> c_str());
          continue;
        }
      }

      from_string<int>(start_index, peptide_start_index_string); 
    }
    // set parent protein of the peptide src
    peptide_src->setParentProtein(parent_protein);

    // set digest type of peptide src
    peptide_src->setDigest(digestion);

    // set start index of peptide src
    peptide_src->setStartIdx(start_index);

    // set current peptide_src to the next empty peptide src
    peptide_src = peptide_src->getNextAssociation();
  } // next peptide_src in file

  carp(CARP_DETAILED_DEBUG, "Done parsing id line:%s", file.getString(PROTEIN_ID_COL).c_str());

  return true;
  } 
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
 * \returns true if peptide_src's were successfully parsed, else
 * returns false.
 */
bool PeptideSrc::parse(
  Peptide* peptide,   ///< assign peptide_src(s) to this peptide
  FILE* file,           ///< file to read from
  Database* database, ///< database containing proteins
  bool use_array) ///< use array implementation vs. linked list
{
  if( peptide == NULL || file == NULL ){
    carp(CARP_ERROR, "Cannot parse peptide src with NULL peptide or file.");
    return false;
  }

  // first field in file should be number of src's
  int num_peptide_src = -1;
  size_t num_read = fread(&num_peptide_src, sizeof(int), 1, file);
  if( num_peptide_src < 1 || num_read != 1 ){
    carp(CARP_ERROR, 
         "Index file corrupted, peptide must have at least one peptide src");
    return false;
  }

  PeptideSrc* peptide_src = NULL;
  // allocate new src based on requested type
  if(use_array){
    peptide_src = newArray(num_peptide_src);
  }else{
    peptide_src = newLinklist(num_peptide_src);
  }

  // give it to the peptide
  peptide->addPeptideSrcArray(peptide_src);

  // read in each peptide_src (prot index, peptide type, start index)
  int src_idx = 0;
  int protein_index = -1;
  Protein* parent_protein = NULL;
  DIGEST_T digestion;
  int start_index = -1;
  for(src_idx = 0; src_idx < num_peptide_src; src_idx++){
    // get protein index
    size_t num_read = fread(&protein_index, (sizeof(int)), 1, file);
    if( protein_index < 0 || num_read != 1){
      carp(CARP_ERROR, "Index file corrupted could not read protein index");
      delete peptide_src;
      return false;
    }
    carp(CARP_DETAILED_DEBUG, "Peptide src %d has protein idx %i", 
         src_idx, protein_index);

    // read peptide type of peptide src
    /*
    if(fread(&peptide_type, sizeof(PEPTIDE_TYPE_T), 1, file) != 1){
      carp(CARP_ERROR, "Index file corrupted, failed to read peptide type.");
      free(peptide_src);
      return false;
    }
    */

    // read digestion level of peptide src
    if(fread(&digestion, sizeof(DIGEST_T), 1, file) != 1){
      carp(CARP_ERROR, "Index file corrupted, failed to read digestion type.");
      delete peptide_src;
      return false;
    }

    // read start index of peptide in parent protein of thsi peptide src
    if(fread(&start_index, sizeof(int), 1, file) != 1){
      carp(CARP_ERROR, "Index file corrupted, failed to read start index.");
      delete peptide_src;
      return false;
    }

    // set fields in new peptide src
    parent_protein = 
      database->getProteinAtIdx(protein_index);
    
    // set parent protein of the peptide src
    peptide_src->setParentProtein(parent_protein);

    // set peptide type of peptide src
    //    set_peptide_src_peptide_type(peptide_src, peptide_type);
    peptide_src->setDigest(digestion);

    // set start index of peptide src
    peptide_src->setStartIdx(start_index);
    
    // set current_peptide_src to the next empty peptide src
    peptide_src = peptide_src->getNextAssociation();
  }// next peptide_src in file

  carp(CARP_DETAILED_DEBUG, "Finished parsing peptide src.");
  return true;
}



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
