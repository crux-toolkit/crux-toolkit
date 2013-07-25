/*************************************************************************//**
 * \file AbstractMatch.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 1/26 2013
 * \brief Object for holding match scores.
 ****************************************************************************/
#include "AbstractMatch.h"
#include "carp.h"
  
/**
 * \returns a new memory allocated abstract match
 */
AbstractMatch::AbstractMatch() {
}

/**
 * default destructor
 */
AbstractMatch::~AbstractMatch() {
}


/**
 * \returns the match score for a particular score type
 */
FLOAT_T AbstractMatch::getScore(
  SCORER_TYPE_T type ///< score type desired
  ) const {

  ScoreMap::const_iterator findScore = scores_.find(type);
  if (findScore == scores_.end()) {
    carp(CARP_WARNING, "Score not set!");
    return 0.0;
  }

  return findScore->second;

}

/**
 * \returns the match rank for a particular score type
 */
int AbstractMatch::getRank(
  SCORER_TYPE_T type ///< score type desired
  ) const {

  RankMap::const_iterator findRank = ranks_.find(type);
  if (findRank == ranks_.end()) {
    carp(CARP_WARNING, "Rank not set!");
    return 0;
  }

  return findRank->second;

}

/**
 * sets the match score for particular score type
 */
void AbstractMatch::setScore(
  SCORER_TYPE_T type, ///< score to set
  FLOAT_T score ///< score value
  ) {

  scores_[type] = score;  

}

/**
 * sets the match rank for particular score type
 */
void AbstractMatch::setRank(
  SCORER_TYPE_T type, ///< rank to set
  int rank ///< rank value
  ) {

  ranks_[type] = rank;

}

/**
 * \returns whether the match has a particular score type assigned
 */
bool AbstractMatch::hasScore(
  SCORER_TYPE_T type ///< score to test
  ) const {

  return (scores_.find(type) != scores_.end());
}

/**
 * \returns whether the match has a particular rank type assigned
 */
bool AbstractMatch::hasRank(
  SCORER_TYPE_T type ///< rank to test
  ) const {

  return (ranks_.find(type) != ranks_.end());
}

/**
 * \returns whether the match is a decoy or not (default false)
 */
bool AbstractMatch::isDecoy() {

  return false;
}
  
/**
 * \returns the beginning iterator for all set scores in the match
 */
ScoreMapIterator AbstractMatch::scoresBegin() {
  return scores_.begin();
}

/**
 * \returns the end iterator for all set scores in the match
 */
ScoreMapIterator AbstractMatch::scoresEnd() {
  return scores_.end();
}

/**
 * \returns the beginning iterator for all set ranks in the match
 */
RankMapIterator AbstractMatch::ranksBegin() {
  return ranks_.begin();
}

/**
 * \returns the end iterator for all set ranks in the match
 */
RankMapIterator AbstractMatch::ranksEnd() {
  return ranks_.end();
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
