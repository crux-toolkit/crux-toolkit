/**
 * \file MatchIterator.cpp 
 * $Revision: 1.38 $
 * \brief An object that iterates over the match objects in the
 * specified match_collection for the specified score type (SP, XCORR)
 */
#include "MatchIterator.h"


/**
 * match_iterator routines!
 *
 */


/**
 * Initializes a match iterator object.
 */
void MatchIterator::init() {

  match_collection_ = NULL;
  match_mode_ = (SCORER_TYPE_T)0;
  match_idx_ = 0;
  match_total_ = 0;
}

MatchIterator::MatchIterator(
  MatchCollection* match_collection
  ) {
    // TODO (BF 06-Feb-08): Could we pass back an iterator with has_next==False
  if (match_collection == NULL){
    carp(CARP_FATAL, "Null match collection passed to match iterator");
  }
  // is there an existing iterator?
  if(match_collection->iterator_lock_){
    carp(CARP_FATAL, 
         "Can only have one match iterator instantiated at a time");
  }
  
  init();
  // set items
  match_collection_ = match_collection;
  match_idx_ = 0;
  match_total_ = match_collection->match_total_;



  // ok lock up match collection
  match_collection->iterator_lock_ = true;
}


/**
 * create a new memory allocated match iterator, which iterates over
 * match iterator only one iterator is allowed to be instantiated per
 * match collection at a time 
 *\returns a new memory allocated match iterator
 */
MatchIterator::MatchIterator(
  MatchCollection* match_collection,
  ///< the match collection to iterate -out
  SCORER_TYPE_T score_type,
  ///< the score type to iterate (LOGP_EXP_SP, XCORR) -in
  bool sort_match  ///< should I return the match in sorted order?
  )
{
  // TODO (BF 06-Feb-08): Could we pass back an iterator with has_next==False
  if (match_collection == NULL){
    carp(CARP_FATAL, "Null match collection passed to match iterator");
  }
  // is there an existing iterator?
  if(match_collection->iterator_lock_){
    carp(CARP_FATAL, 
         "Can only have one match iterator instantiated at a time");
  }
  
  // has the score type been populated in match collection?
  if(!match_collection->scored_type_[score_type]){
    const char* score_str = scorer_type_to_string(score_type);
    carp(CARP_ERROR, "New match iterator for score type %s.", score_str);
    carp(CARP_FATAL, 
         "The match collection has not been scored for request score type.");
  }
  
  init();
  // set items
  match_collection_ = match_collection;
  match_mode_ = score_type;
  match_idx_ = 0;
  match_total_ = match_collection->match_total_;

  // only sort if requested and match collection is not already sorted
  if(sort_match){
    match_collection->sort(score_type);
  }


  // ok lock up match collection
  match_collection->iterator_lock_ = true;
}

/**
 * \brief Create a match iterator to return matches from a collection
 * grouped by spectrum and sorted by given score type.
 *
 * \returns A heap-allocated match iterator.
 */
/*
MatchIterator::MatchIterator(
  MatchCollection* match_collection,  ///< for iteration -in
  SCORER_TYPE_T scorer ///< the score type to sort by -in
){

  init();

  // set up fields
  match_collection_ = match_collection;
  match_mode_ = scorer;
  match_idx_ = 0;
  match_total_ = match_collection->match_total_;

  match_collection->spectrumSort(scorer);

  match_collection->iterator_lock_ = true;
}
*/
/**
 * free the memory allocated iterator
 */
MatchIterator::~MatchIterator()
{
  if (match_collection_ != NULL){
    match_collection_->iterator_lock_ = false;
  }
}



/**
 * Does the match_iterator have another match struct to return?
 *\returns true, if match iter has a next match, else False
 */
bool MatchIterator::hasNext()
{
  return (match_idx_ < match_total_);
}

/**
 * return the next match struct!
 *\returns the match in decreasing score order for the match_mode(SCORER_TYPE_T)
 */
Crux::Match* MatchIterator::next()
{
  return match_collection_->match_[match_idx_++];
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
