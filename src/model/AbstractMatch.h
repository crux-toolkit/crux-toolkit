/**
 * \file AbstractMatch.h
 * $Revision: 1.00 $ 
 * \brief Object for holding match scores
 ****************************************************************************/
#ifndef ABSTRACTMATCH_H_
#define ABSTRACTMATCH_H_

#include "match_objects.h"

class AbstractMatch {

 protected:
  ScoreMap scores_;  ///< map of scores
  RankMap ranks_; ///< map of ranks for scores

 public:

  /**
   * \returns a new memory allocated abstract match
   */
  AbstractMatch();

  /**
   * default destructor
   */
  virtual ~AbstractMatch();

  /**
   * \returns the match score for a particular score type
   */
  virtual FLOAT_T getScore(
    SCORER_TYPE_T type ///<score type desired
    ) const;

  /**
   * \returns the match rank for a particular score type
   */
  virtual int getRank(
    SCORER_TYPE_T type ///<score type desired
    ) const;

  /**
   * sets the match score for particular score type
   */
  virtual void setScore(
    SCORER_TYPE_T type, ///< score to set
    FLOAT_T score ///< score value
  );

  /**
   * sets the match rank for particular score type
   */
  virtual void setRank(
    SCORER_TYPE_T type, ///< rank to set
    int rank ///< rank value
  );
  
  /**
   * \returns whether the match has a particular score type assigned
   */
  virtual bool hasScore(
    SCORER_TYPE_T type ///< score to test
  ) const;

  virtual bool hasRank(
    SCORER_TYPE_T type ///< rank to test
  ) const;

  /**
   * \returns whether the match is a decoy or not (default false).
   */
  virtual bool isDecoy();  

  /**
   * \returns the beginning iterator for all set scores in the match
   */
  ScoreMapIterator scoresBegin();

  /**
   * \returns the end iterator for all set scores in the match
   */
  ScoreMapIterator scoresEnd();

  /**
   * \returns the beginning iterator for all set ranks in the match
   */
  RankMapIterator ranksBegin();

  /**
   * \returns the end iterator for all set ranks in the match
   */
  RankMapIterator ranksEnd();

};

#endif //ABSTRACTMATCH_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
