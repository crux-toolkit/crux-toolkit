/**
 * \file MatchIterator.h 
 * $Revision: 1.38 $
 * \brief An object that iterates over the match objects in the
 * specified match_collection for the specified score type (SP, XCORR)
 */
#ifndef MATCHITERATOR_H
#define MATCHITERATOR_H

#include "MatchCollection.h"

class MatchIterator {
 protected:
  MatchCollection* match_collection_; ///< the match collection to iterate -out
  SCORER_TYPE_T match_mode_; ///< the current working score (SP, XCORR)
  int match_idx_;            ///< current match to return
  int match_total_;          ///< total_match_count

  /**
   * Initializes a match iterator object.
   */
  void init();

 public:

  MatchIterator(
    MatchCollection* match_collection ///< the match collection to iterate -in
    );

  /**
   * create a new memory allocated match iterator
   * creates a new the generate_peptides_iterator inside the match_iterator
   *\returns a new memory allocated match iterator
   */
  MatchIterator(
    MatchCollection* match_collection, ///< the match collection to iterate -in
    SCORER_TYPE_T match_mode, ///< the mode to set (MATCH_SP, MATCH_XCORR) -in
    bool sort_match = false  ///< should I return the match in sorted order? (default false)
    );

  /**
   * \brief Create a match iterator to return matches from a collection
   * grouped by spectrum and sorted by given score type.
   *
   * \returns A heap-allocated match iterator.
   */
  /*
  MatchIterator(
    MatchCollection* match_collection, ///< match collection to iterate -in
    SCORER_TYPE_T scorer ///< the score to sort by (MATCH_SP, MATCH_XCORR) -in
  );
  */
  /**
   * free the memory allocated iterator
   */
  ~MatchIterator();

  /**
   * Does the match_iterator have another match struct to return?
   * MUST set the iterator to correct mode before calling this method
   *\returns true, if match iter has a next match, else False
   */
  bool hasNext();

  /**
   * return the next match struct!
   * MUST set the iterator to correct mode before initialially calling this method
   *\returns the match in decreasing score order for the match_mode(SCORER_TYPE_T)
   */
  Crux::Match* next();

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
