/**
 * \file DatabasePeptideIterator.cpp
 * \brief Object to iterator over the peptides within a database
 *****************************************************************************/

#include "DatabasePeptideIterator.h"

using namespace std;
using namespace Crux;

/***********************************************
 * database_peptide_Iterators - can use the light protein functionality 
 * to save space
 ***********************************************/

/**
 * Instantiates a new database_peptide_iterator from a database.
 * \returns a DATABASE_PEPTIDE_ITERATOR_T object.
 */

DatabasePeptideIterator::DatabasePeptideIterator(
  Database* database, 
    ///< the database of interest -in
  PeptideConstraint* peptide_constraint,
    ///< the peptide_constraint with which to filter peptides -in
  bool store_all_peptides, ///< for removing duplicates
  bool is_decoy ///< return decoy instead of target peptides
  )
{
  Protein* next_protein = NULL;

  // Initialize
  database_protein_iterator_ = NULL;
  cur_protein_peptide_iterator_ = NULL;
  peptide_constraint_ = NULL;
  prior_protein_ = NULL;
  first_passed_ = false;
  store_all_peptides_ = false;
  is_decoy_ = is_decoy;
  
  
  // set up peptide storage
  store_all_peptides_ = store_all_peptides;
  peptide_map_ = map<char*, Peptide*, cmp_str>();
  cur_map_position_ = peptide_map_.begin();

  // set a new protein iterator
  database_protein_iterator_ = new DatabaseProteinIterator(database);
  if( database_protein_iterator_ == NULL){
    carp(CARP_ERROR, 
         "Could not create protein iterator for database peptide iterator.");
  }

  // set peptide constraint
  peptide_constraint_ = PeptideConstraint::copyPtr(peptide_constraint);
  
  // check if there are any proteins to create peptides from
  if(database_protein_iterator_->hasNext()){

    next_protein = database_protein_iterator_->next();

    // if using light/heavy functionality parse the light protein
    if(database->getUseLightProtein() && next_protein->getIsLight()){
      if(!next_protein->toHeavy()){
        carp(CARP_FATAL, "failed to create a database_peptide_iterator,"
                         "no proteins in database");
      }
    }

    // set new protein peptide iterator
    cur_protein_peptide_iterator_ =
      new ProteinPeptideIterator(next_protein, peptide_constraint);
 
    // if first protein does not contain a match peptide, reinitailize
    while(!cur_protein_peptide_iterator_->hasNext()) {
      // covert the heavy back to light
      /** 
       * uncomment this code if you want to restore a protein to 
       * light after converted to heavy
      if(database->use_light_protein && !get_protein_is_light(next_protein)){
        protein_to_light(next_protein);
      }
      */

      // end of list of peptides for database_peptide_iterator
      if(!database_protein_iterator_->hasNext()){
        break;
      }
      else{ // create new protein_peptide_iterator for next protein
        // free old iterator
            delete (cur_protein_peptide_iterator_);
        
        // get next protein
        next_protein = database_protein_iterator_->next();

         // if using light/heavy functionality parse the light protein
        if(database->getUseLightProtein() && next_protein->getIsLight()){
          if(!next_protein->toHeavy()){
            carp(CARP_FATAL, "failed to create a database_peptide_iterator"
                              " no proteins in database");
          }
        }
        // create new protein_peptide_iterator
        cur_protein_peptide_iterator_ =
          new ProteinPeptideIterator(next_protein, peptide_constraint);
      }
    }
  }
  else{ // no proteins to create peptides from
    carp(CARP_FATAL, "failed to create a database_peptide_iterator,"
                     "no proteins in database");
  }
  // set the current working protein
  prior_protein_ = next_protein;
  
  // queue up first peptide
  if( store_all_peptides ){
    generateAllPeptides(); // fill the map
    queueFirstPeptideFromMap();

  } else {  
    next_peptide_ = nextFromFile();
  }
  // parent class requires a call to queueNextPeptide via initialize(),
  // but the first peptide is already queued.  This is a lazy way to
  // get around it.
  already_initialized_ = false;
  initialize();
}

/**
 * Private function to generate peptides for this iterator, sort them
 * and remove duplicates before returning them to the caller.
 * Peptides are stored in the iterator's member variable peptide_map.
 */
void DatabasePeptideIterator::generateAllPeptides(){

  Peptide* cur_peptide = NULL;
  map<char*, Peptide*, cmp_str>& peptide_map = peptide_map_;
  map<char*, Peptide*, cmp_str>::iterator peptide_map_ptr;

  // populate map with all peptides, combining when duplicates found
  while(hasNextFromFile()){
    cur_peptide = nextFromFile();
    char* sequence = cur_peptide->getSequence();

    // does it already exist in the map?
    peptide_map_ptr = peptide_map.find(sequence);
    if( peptide_map_ptr == peptide_map.end() ){ // not found, add it
      peptide_map[sequence] = cur_peptide;
    } else {  // already exists, combine new peptide with existing
      Peptide::mergePeptidesCopySrc(peptide_map_ptr->second, cur_peptide);
      delete cur_peptide; 
      free(sequence); 
    }
  } // next peptide

}

/**
 * Once the peptide_map has been populated with peptides, point the
 * cur_peptide to first in map, point map pointer to second in map.
 */
void DatabasePeptideIterator::queueFirstPeptideFromMap(){

  if( peptide_map_.empty() ){  // no peptides to return
    next_peptide_ = NULL;
    cur_map_position_ = peptide_map_.end();
    return;
  }

  // set cur_peptide to first peptide in map
  next_peptide_ = peptide_map_.begin()->second;
  // set map pointer to one past first (i.e. next to return)
  cur_map_position_ = peptide_map_.begin();  
  ++(cur_map_position_);
}

/**
 * Frees an allocated database_peptide_iterator object.
 */
DatabasePeptideIterator::~DatabasePeptideIterator() {

  delete cur_protein_peptide_iterator_;
  delete database_protein_iterator_;
  PeptideConstraint::free(peptide_constraint_);

  // free seqs in map
  map<char*, Peptide*, cmp_str>::iterator peptide_iter = 
    peptide_map_.begin();
  for(; peptide_iter != peptide_map_.end(); 
      ++peptide_iter){

    free (peptide_iter->first);  // free sequence
    // peptide freed elsewhere?  segfault if here
  }
  peptide_map_.clear();
}

/**
 * The basic iterator functions.
 * \returns true if there are additional peptides to iterate over,
 * false if not. 
 */
bool DatabasePeptideIterator::hasNextFromFile() {

  return (cur_protein_peptide_iterator_->hasNext());
}

/**
 * Implementation of PeptideIterator's method to prepare the iterator
 * to return the next peptide.
 * \returns True if there is another peptide to return, else false.
 */
bool DatabasePeptideIterator::queueNextPeptide() {

  // parent class will make the first call without returning any peptides
  // but we already queued up the first one in the constructor
  if( already_initialized_ == false ){
    already_initialized_ = true;
    if(is_decoy_ && next_peptide_){
      next_peptide_->transformToDecoy();
    }
    return (next_peptide_ != NULL); // are we starting out with a peptide
  }
  bool has_next = false;

  // fetch next peptide either from file or from map
  if( store_all_peptides_ ){
    if( cur_map_position_ == peptide_map_.end() ){
      next_peptide_ = NULL;
      has_next = false;
    } else {
      next_peptide_ = cur_map_position_->second;
      ++(cur_map_position_);
      has_next = true;
    }
  } else {
    if( hasNextFromFile() ){
      next_peptide_ = nextFromFile();
      has_next = (next_peptide_ != NULL);
    } else {
      next_peptide_ = NULL;
      has_next = false;
    }
  }

  // take care of decoys
  if(is_decoy_ && next_peptide_){
    next_peptide_->transformToDecoy();
  }

  return has_next;
}

/**
 * \returns The next peptide in the database.
 */
Peptide* DatabasePeptideIterator::nextFromFile()
{
  /*BF: Could this be simplified?  if next peptide, return it
    if not, look for next protein, if not return NULL
  */

   // did you reset working protein?
  bool reset = false;
  
  // the peptide to return
  Peptide* next_peptide =
    cur_protein_peptide_iterator_->next();
  
  Database* database = 
    database_protein_iterator_->getDatabase();
  
  // reset database_peptide_iterator if needed
  while(!cur_protein_peptide_iterator_->hasNext()){
    reset = true;
    Protein* next_protein = NULL; 
    
    /** 
     * uncomment this code if you want to restore a protein to 
     * light after converted to heavy
    // covert the heavy back to light
    if(database->use_light_protein && next_protein != NULL && !get_protein_is_light(next_protein)){
      protein_to_light(next_protein);
    }
    */

    // end of list of peptides for database_peptide_iterator
    if(!database_protein_iterator_->hasNext()){
      break;
    }
    else{ // create new protein_peptide_iterator for next protein
      // free old iterator
      delete cur_protein_peptide_iterator_;
      
      // get next protein
      next_protein = 
        database_protein_iterator_->next();
      
      // if using light/heavy functionality parse the light protein
      if(database->getUseLightProtein() && next_protein->getIsLight()){
        if(!next_protein->toHeavy()){
          carp(CARP_FATAL, "Failed to create a database_peptide_iterator, " 
                            "no proteins in database");
        }
      }
      // create new protein_peptide_iterator
      cur_protein_peptide_iterator_ =
        new ProteinPeptideIterator(next_protein, peptide_constraint_);
    }        
  }

  // are we using the light functionality?
  if(database->getUseLightProtein()){
    // get the current working protein
    Protein* protein_bye = cur_protein_peptide_iterator_->getProtein();
    // set first passed, shows that we extraced at least one protein since we moved on to the next protein
    if(!reset && !first_passed_){
      first_passed_ = true;
    }
    // convert prior_protein to light
    else if(!reset && first_passed_ && (protein_bye != prior_protein_)){
      /** 
       * uncomment this code if you want to restore a protein to 
       * light after converted to heavy
      protein_to_light(database_peptide_iterator->prior_protein);
      */
      prior_protein_ = protein_bye;
    }
  }
  
  return next_peptide;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
