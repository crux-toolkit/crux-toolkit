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
#include "MatchCollectionParser.h"


using namespace std;
using namespace Crux;

/**
 * Static variable definitions
 */
map<string, Peptide* > PeptideSrc::sequence_to_peptide_; ///< Maps a sequence to a peptide object
map<string, Peptide* > PeptideSrc::decoy_sequence_to_peptide_; ///< Maps a decoy sequence to a peptide object


/**
 * \returns An (empty) peptide_src object.
 */
PeptideSrc::PeptideSrc() {
  digestion_ = (DIGEST_T)0;
  parent_protein_ = NULL;
  start_idx_ = 0;
  start_idx_original_ = 0;
}

/**
 *\returns a PeptideSrc object, populated with user specified parameters
 */
PeptideSrc::PeptideSrc(
  DIGEST_T digest,
  Protein* parent_protein, ///< the parent of this peptide -in
  int start_idx ///< start index of the peptide in the protein sequence -in
  ) {

  setDigest(digest);
  setParentProtein(parent_protein);
  setStartIdx(start_idx);
}


/**
 * Frees the entire allocated peptide_srcs object
 */
void PeptideSrc::free(vector<PeptideSrc*>& peptide_srcs) {


  for(vector<PeptideSrc*>::iterator iter = peptide_srcs.begin();
     iter != peptide_srcs.end();
     ++iter) {
  delete *iter; 

  }
  peptide_srcs.clear();
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
  vector<PeptideSrc*>& src, ///< source peptide_src -in
  vector<PeptideSrc*>& dest ///< destination peptide_src -out
  )
{

  for (vector<PeptideSrc*>::iterator iter = src.begin();
       iter != src.end();
       ++iter) {
    PeptideSrc* dest_src = new PeptideSrc(*(*iter));
    dest.push_back(dest_src);

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
 * \sets the original start index of the peptide in the protein sequence
 */
void PeptideSrc::setStartIdxOriginal(
  int start_idx ///< start index of the peptide in the original protein sequence -in
) {
  start_idx_original_ = start_idx;
}

/**
 * \returns the original start index of the peptide in the protein sequence
 */
int PeptideSrc::getStartIdxOriginal() {
  return start_idx_original_;
}

/**
 * \returns a pointer to the start of the peptide with in it's parent protein sequence
 */
char* PeptideSrc::getSequencePointer() {

  return parent_protein_->getSequencePointer(start_idx_ - 1);

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
  carp(CARP_DETAILED_DEBUG, "protein idx to write is %i", protein_idx);

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
 * \brief Read in the peptide_src objects from the given file and
 * assosiated them with the given peptide.  
 * Proteins for the pepitde_src are found in the given database.  If
 * database is NULL, does not set proteins.  (This option is used for
 * sorting index files while creating index.) 
 *
 * \returns true if peptide_src's were successfully parsed, else
 * returns false.
 */
bool PeptideSrc::parseTabDelimited(
  Peptide* peptide,   ///< assign peptide_src(s) to this peptide
  MatchFileReader& file,           ///< file to read from
  Database* database, ///< database containing proteins
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
    vector<PeptideSrc*>& src_peptide_srcs = src_peptide->getPeptideSrcVector();

    for (vector<PeptideSrc*>::iterator iter = src_peptide_srcs.begin();
      iter != src_peptide_srcs.end();
      ++iter) {
      PeptideSrc* new_src = *iter;
      peptide->addPeptideSrc(new_src);
    }

    return true;


  } else {

    vector<string> protein_ids;
    file.getStringVectorFromCell(PROTEIN_ID_COL, protein_ids);
  
    if (protein_ids.size() == 0) {
      carp(CARP_ERROR, "No protein ids found!");
      return false;
    }

    vector<string> flanking_aas;
    file.getStringVectorFromCell(FLANKING_AA_COL, flanking_aas);

    if (protein_ids.size() != flanking_aas.size()) {
      carp_once(CARP_DEBUG, 
                "Flanking amino acid count (%d) did not match protein count (%d) for protein %s.", 
                flanking_aas.size(), 
                protein_ids.size(),
                file.getString(PROTEIN_ID_COL).c_str());
      carp_once(CARP_DEBUG, "Only reporting error once; others may exist")
      while(flanking_aas.size() < protein_ids.size()) {
        flanking_aas.push_back("");
      }
    }

    //For every protein id source, create the object and add it to the list.
    for (size_t idx = 0; idx < protein_ids.size(); idx++) {
    
      PeptideSrc* peptide_src = new PeptideSrc();
      DIGEST_T digestion = 
        string_to_digest_type((char*)file.getString(CLEAVAGE_TYPE_COL).c_str()); 
  

      Protein* parent_protein = NULL;
      int start_index = 1;

      string protein_id = protein_ids.at(idx);
      string flanking_aa = flanking_aas.at(idx);
      string prev_aa = "", next_aa = "";
      if (flanking_aa.length() == 2) {
        prev_aa = flanking_aa[0];
        next_aa = flanking_aa[1];
      }

      carp(CARP_DETAILED_DEBUG,"Parsing %s", protein_id.c_str());
      // get the protein and peptide index e.g. X(10)
      size_t left_paren_index = protein_id.find('(');

      if (left_paren_index == string::npos) {
        //protein id is the string.
        bool is_decoy;

        parent_protein=MatchCollectionParser::getProtein(
          database, decoy_database, protein_id, is_decoy);
        if (parent_protein == NULL) {
          carp(CARP_WARNING, "Can't find protein %s",protein_id.c_str());
          continue;
        }


        //find the start index
        MODIFIED_AA_T* mod_seq;
        int seq_length = convert_to_mod_aa_seq(file.getString(SEQUENCE_COL).c_str(), &mod_seq);
        char* unmodified_sequence = modified_aa_to_unmodified_string(mod_seq, seq_length);
        string sequence = unmodified_sequence;
        std::free(unmodified_sequence);
        std::free(mod_seq);

        start_index = parent_protein->findStart(sequence, prev_aa, next_aa);
        if (start_index == -1) {
          carp(CARP_FATAL, "Can't find sequence %s in %s:%s",
            sequence.c_str(),
            protein_id.c_str());
        }
      } else {
        string protein_id_string = protein_id.substr(0, left_paren_index);
        string peptide_start_index_string = protein_id.substr(left_paren_index+1, 
          protein_id.length() - 1);

        bool is_decoy;    
        //  set fields in new peptide src
        parent_protein = MatchCollectionParser::getProtein(
          database, decoy_database, protein_id_string, is_decoy);

        MODIFIED_AA_T* mod_seq;
        int seq_length = convert_to_mod_aa_seq(file.getString(SEQUENCE_COL).c_str(), &mod_seq);
        char* unmodified_sequence = modified_aa_to_unmodified_string(mod_seq, seq_length);
        string sequence = unmodified_sequence;
        std::free(unmodified_sequence);
        std::free(mod_seq);

        //string sequence = file.getString(SEQUENCE_COL);
        start_index = parent_protein->findStart(sequence, prev_aa, next_aa);
      }

      // set parent protein of the peptide src
      peptide_src->setParentProtein(parent_protein);

      // set digest type of peptide src
      peptide_src->setDigest(digestion);

      // set start index of peptide src
      peptide_src->setStartIdx(start_index);

      // add it to the list of protein sources for this peptide
      peptide->addPeptideSrc(peptide_src);

    }
  } // next peptide_src in file

  carp(CARP_DETAILED_DEBUG, "Done parsing id line:%s", file.getString(PROTEIN_ID_COL).c_str());

  return true;
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
 * sorting index files while creating index.)  
 * \returns true if peptide_src's were successfully parsed, else
 * returns false.
 */
bool PeptideSrc::parse(
  Peptide* peptide,   ///< assign peptide_src(s) to this peptide
  FILE* file,           ///< file to read from
  Database* database) ///< database containing proteins
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

  // read in each peptide_src (prot index, peptide type, start index)
  int src_idx = 0;
  int protein_index = -1;
  Protein* parent_protein = NULL;
  DIGEST_T digestion;
  int start_index = -1;
  for(src_idx = 0; src_idx < num_peptide_src; src_idx++){

    PeptideSrc* peptide_src = new PeptideSrc();

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

    peptide->addPeptideSrc(peptide_src);

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
