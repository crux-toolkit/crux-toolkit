/**
 * \file IndexFilteredPeptideIterator.cpp 
 * $Revision: 1.25 $
 * \brief Object for iterating peptides from an index.  Unlike a
 * standard index iterator, it filters by digestion type
 *****************************************************************************/
#include "IndexFilteredPeptideIterator.h"
#include <vector>

using namespace std;
using namespace Crux;

/**
 * Instantiates a new index_filtered_peptide_iterator from a index.
 * \returns a new heap allocated index_filtered_peptide_iterator object
 */
IndexFilteredPeptideIterator::IndexFilteredPeptideIterator(
    Index* index ///< The index object which we are iterating over -in
 ) : IndexPeptideIterator(index, NULL)
{
  setup();
}

/**
 * Frees an allocated index_filtered_peptide_iterator object.
 */
IndexFilteredPeptideIterator::~IndexFilteredPeptideIterator(){
  // if did not iterate over all peptides, free the last peptide not returned
  if(has_next_){
    delete peptide_;
  }
}

/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
Peptide* IndexFilteredPeptideIterator::next(){
  Peptide* peptide_to_return = peptide_;
  setup();
  return peptide_to_return;
}

/**
 * The basic iterator functions.
 * Does the iterator have more peptides to return.
 *\returns True if there are additional peptides, false if not.
 */
bool IndexFilteredPeptideIterator::hasNext(){
  return has_next_;
}

/**
 * Sets up the iterator.
 * \returns True if successfully sets up the iterator, else false.
 */
bool IndexFilteredPeptideIterator::setup()
{
  Peptide* peptide = NULL;
  DIGEST_T required_digestion = index_->getSearchConstraint()->getDigest();
  bool match = false;
  
  // initialize index_filered
  while(IndexPeptideIterator::hasNext()){
    peptide = IndexPeptideIterator::next();
    vector<PeptideSrc*>& srcs = peptide->getPeptideSrcVector();
    // mass, length has been already checked in index_peptide_iterator
    // check if peptide type matches the constraint
    // find at least one peptide_src for which cleavage is correct
    for (size_t idx = 0; idx < srcs.size(); idx++) {
      if(srcs[idx]->getDigest() >= required_digestion){
        match = true;
        break;
      }
    }
    
    // add more filters to the peptides here, if they don't meet
    // requirements change 'match' to false 
    
    // this peptide meets the peptide_type
    if(match){
      peptide_ = peptide;
      has_next_ = true;
      return true;
    }
    delete peptide;
  }// next peptide
  // no peptides meet the constraint
  has_next_ = false;
  peptide_ = NULL;
  return true;
}

/**********************************************************************
 * wrapper, for generate_peptides_iterator, cast back to original type
 ***********************************************************************/
/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
Peptide* void_index_filtered_peptide_iterator_next(
  void* iterator ///< the index_filtered_peptide_iterator to initialize -in
  )
{
  return ((IndexFilteredPeptideIterator*)iterator)->next();

}

/**
 * The basic iterator functions.
 *\returns True if there are additional peptides to iterate over, false if not.
 */
bool void_index_filtered_peptide_iterator_has_next(
  void* iterator ///< the index_filtered_peptide_iterator to initialize -in
  )
{
  return ((IndexFilteredPeptideIterator*)iterator)->hasNext();
}

/**
 * Frees an allocated index_filtered_peptide_iterator object.
 */
void void_free_index_filtered_peptide_iterator(
  void* iterator ///< the index_filtered_peptide_iterator to initialize -in
  )
{
  delete (IndexFilteredPeptideIterator*)iterator;
}

/**
 * Frees an allocated index_peptide_iterator object.
 */
void void_free_index_peptide_iterator(
    void* index_peptide_iterator  ///< the iterator to free -in
    )
{
  delete (IndexPeptideIterator*)index_peptide_iterator;
}

/**
 * The basic iterator functions.
 * \returns True if there are additional peptides to iterate over, false if not.
 */
bool void_index_peptide_iterator_has_next(
    void* index_peptide_iterator ///< the iterator of interest -in
    )
{

  return ((IndexPeptideIterator*)index_peptide_iterator)->hasNext();
}

/**
 * \returns The next peptide in the index.
 */
Peptide* void_index_peptide_iterator_next(
    void* index_peptide_iterator ///< the iterator of interest -in
    )
{

  return ((IndexPeptideIterator*)index_peptide_iterator)->next();
}



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
