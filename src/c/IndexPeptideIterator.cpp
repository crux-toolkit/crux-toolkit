/************************************************************************//**
 * \file IndexPeptideIterator.cpp
 * \brief Object for representing iterating peptides from an index.
 ****************************************************************************/

#include "IndexPeptideIterator.h"
#include <errno.h>
#include "WinCrux.h"

using namespace Crux;
/***********************************************
 *  The basic index_peptide_iterator functions.
 ***********************************************/

/**
 * Instantiates a new index_peptide_iterator from an index.
 * \returns a new heap allocated index_peptide_iterator object
 */
IndexPeptideIterator::IndexPeptideIterator(
  Index* index, ///< The index object which we are iterating over -in
  PeptideConstraint* constraint, ///< the peptide constraint
  bool is_decoy ///< return target or decoy peptides
  )
  : is_decoy_(is_decoy)
{
  carp(CARP_DETAILED_DEBUG, "Creating new %sindex iterator", 
       ((is_decoy) ? "decoy " : ""));


  // initialize a new index_peptide_iterator object
  index_ = NULL;
  current_index_file_ = 0;
  index_file_ = NULL;
  
  // set index
  index_ = index->copyPtr();
  if (constraint != NULL) {
    constraint_ = PeptideConstraint::copyPtr(constraint);
  }

  //Find the index files that match our constraint
  index_->getIndexMap(is_decoy_)->getIndexFiles(constraint_, index_files_);

  // set remaining iterator fields
  current_index_file_ = 0;
  index_file_ = NULL;

  // sets has_next, peptide, index_file
  carp(CARP_DETAILED_DEBUG, "IndexPeptideIterator queueing first peptide");
  initialize();
}

/**
 * Frees an allocated index_peptide_iterator object.
 */
IndexPeptideIterator::~IndexPeptideIterator()
{
  // if did not iterate over all peptides, free the last peptide not returned
  if(hasNext()){
    delete next_peptide_;
  }
  
  // free the index
  Index::free(index_);

  if (constraint_) {
    PeptideConstraint::free(constraint_);
  }

}

/**
 * \brief Prepare an index peptide iterator to have its index files
 * searched.  Changes the state of the iterator if its index_file
 * is NULL and the current file is not the last in the list.  The
 * routine that searches for peptides in a file closes any files it
 * has read to the end and sets the pointer to NULL.
 * \returns true if there is a file ready to be read or false if no
 * more files remain.
 */
bool IndexPeptideIterator::findNextIndexFile(){

  carp(CARP_DETAILED_DEBUG, "Finding file");
  // file is ready to read
  if( index_file_ != NULL ){
    carp(CARP_DETAILED_DEBUG, "Current file ready to read");
    return true;
  }
  // no more files to open
  if( current_index_file_ == index_files_.size() ){
  carp(CARP_DETAILED_DEBUG, "The last file has been opened. no more");
    return false;
  }
  // else, current is NULL and there are more to open
  char* filename=index_files_[current_index_file_]->getFilename();
  carp(CARP_DETAILED_DEBUG, "Opening new file %s", filename);
  index_file_ = fopen(filename, "r");

  if( index_file_ == NULL){
    carp(CARP_ERROR, "Could not open index file %s", filename);
    return false;
  }

  return true;

}
  


/**
 * \brief Search for a peptide matching the peptide constraint in the
 * current index file of the iterator.  If the iterator has no open
 * index file, returns false.  If the current file reaches the end
 * without finding a suitable peptide, the file is closed, file handle
 * set to NULL, and the number of the current file is increemnted.
 *
 * \returns true if the iterator is ready to return a new peptide,
 * false if no peptide meeting the constraint could abe found in the
 * current file. 
 */
bool IndexPeptideIterator::findPeptideInCurrentIndexFile()
{
  carp(CARP_DETAILED_DEBUG, "Looking for peptide in current file");

  FILE* cur_file = index_file_;
  if( cur_file == NULL ){
    carp(CARP_DETAILED_DEBUG, "current file is null");
    return false;
  }

  // peptide to return, reuse this memory while we look
  Peptide* peptide = new Peptide();
  // constraint to meet
  PeptideConstraint* index_constraint = constraint_;
  if (index_constraint == NULL) {
    carp(CARP_WARNING, "NULL constriant, getting from INDEX");
    index_constraint = index_->getSearchConstraint();
  }
  // loop until we get to a peptide that fits the constraint, 
  // a peptide bigger (mass) than the constraint, or reach eof
  bool peptide_fits = false;
  bool file_finished = false;
  long int src_loc = 0;  // in case we need to parse the peptide src

  while( !peptide_fits && !file_finished ){// until pep_fits or file done
    carp(CARP_DETAILED_DEBUG, "Once around the find peptide loop %i, %i", 
         peptide_fits, file_finished);

    // read in next peptide
    bool found_pep = peptide->parseNoSrc(cur_file, &src_loc);
    // returns false if eof
    if( ! found_pep ){
      carp(CARP_DETAILED_DEBUG, "parse peptide returned false");
      file_finished = true;
      continue;
    }

    // check our peptide to see if it fits the constraint
    FLOAT_T peptide_mass = peptide->getPeptideMass();
    int peptide_length = peptide->getLength();

    // if peptide mass larger than constraint, no more peptides to return
    if(peptide_mass > index_constraint->getMaxMass()){
      carp(CARP_DETAILED_DEBUG, "peptide found is bigger than constraint");
      file_finished = true;
      continue;
    }// use else if?

    // does this peptide fit the constraint?

    double min_mass = index_constraint->getMinMass();
    int min_length = index_constraint->getMinLength();
    int max_length = index_constraint->getMaxLength();


    /*
  carp(CARP_DETAILED_DEBUG, "Index search constraints: %i-%i, %.2f-%.2f", get_peptide_constraint_min_length(iterator->index->search_constraint), get_peptide_constraint_max_length(iterator->index->search_constraint), get_peptide_constraint_min_mass(iterator->index->search_constraint), get_peptide_constraint_max_mass(iterator->index->search_constraint));
  carp(CARP_DETAILED_DEBUG, "Index search constraints var: %i-%i, %.2f-2f", min_length, max_length, min_mass); //, max_mass);
    */


    if( peptide_mass >= min_mass ){
      carp(CARP_DETAILED_DEBUG, "peptide mass bigger than min mass.");
    }else{
      carp(CARP_DETAILED_DEBUG, "peptide mass %f smaller than min mass %f.",
           peptide_mass, min_mass);
    }
    if( peptide_length >= min_length ){
      carp(CARP_DETAILED_DEBUG, "peptide length bigger than min length.");
    }
    if( peptide_length <= max_length ){
      carp(CARP_DETAILED_DEBUG, "peptide length smaller than max length.");
    }

    if( peptide_mass >= min_mass
        && peptide_length >= min_length
        && peptide_length <= max_length){
      carp(CARP_DETAILED_DEBUG, "peptide passes constraint.");
      peptide_fits = true;
    }else{
      carp(CARP_DETAILED_DEBUG, "peptide doesn't pass constraint. " \
           "pep len,mass: %i, %.2f, min: %.2f, len: %i-%i", 
           peptide_length, peptide_mass, min_mass, min_length, max_length);
    }

  }// read next peptide

  // we have a peptide to return, get the peptide_src for it
  long int pep_end = 0;
  if( peptide_fits ){
    carp(CARP_DETAILED_DEBUG, "Found a peptide that fits constraint");
    pep_end = ftell(cur_file);
    fseek(cur_file, src_loc, SEEK_SET);
    Database* database = index_->getDatabase(is_decoy_);
    if( ! PeptideSrc::parse(peptide, cur_file, database) ){
      carp(CARP_ERROR, "Could not parse peptide src");
      file_finished = true; // maybe we could read more, but unlikly
    }
  }

  // we broke out of the loop because we don't need to read the file anymore
  if( file_finished  ){
      carp(CARP_DETAILED_DEBUG, "Done with this index file.");
      fclose(cur_file);
      index_file_ = NULL;
      current_index_file_ += 1;
      next_peptide_ = NULL;

      carp(CARP_DETAILED_DEBUG, "about to free peptide.");
      delete peptide;
      carp(CARP_DETAILED_DEBUG, "Done cleaning up.");
      return false;
  }

  // else, everything worked!  finish setting fields

  //return file pointer to end of peptide
  fseek(cur_file, pep_end, SEEK_SET);
  index_file_ = cur_file;

  // add peptide to iterator
  next_peptide_ = peptide;

  return true;
}

/**
 * \brief Find the next peptide for the iterator to return.
 *
 * Called in the process of initializing a new iterator and by
 * next().  Looks in the current index_file for peptides matching the
 * constraint.  If none found, checks remaining index_files.  If no
 * constraint-satisfying peptides are found in any of the remaining
 * index files, next_peptide is set to NULL and has_next is set to false.
 * \returns true if no errors were encountered while reading files
 * (even if there is no peptide to return).
 */
bool IndexPeptideIterator::queueNextPeptide() {

  bool found = false;
  while(findNextIndexFile()){
    found = findPeptideInCurrentIndexFile();
    if(found == true ){
      carp(CARP_DETAILED_DEBUG, "Found returned true");

      // set peptide, set has_next done in find
      return true;
    }
  }// try the next index file

  // no more index files to try
  next_peptide_ = NULL;
  return false;
}

Index* IndexPeptideIterator::getIndex() {
  
  return index_;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
