#ifndef MATCHOBJECTS_H_
#define MATCHOBJECTS_H_

#include "objects.h"
#include <deque>
#include <map>
#include <string>

class ProteinMatchCollection;
class ProteinMatch;
class PeptideMatch;
class SpectrumMatch;
class AbstractMatch;

typedef std::deque<SpectrumMatch*> SpectrumMatchCollection;
typedef SpectrumMatchCollection::iterator SpectrumMatchIterator;

typedef std::deque<PeptideMatch*> PeptideMatchCollection;
typedef PeptideMatchCollection::iterator PeptideMatchIterator;

typedef std::deque<ProteinMatch*>::iterator ProteinMatchIterator;

typedef std::map<SCORER_TYPE_T, FLOAT_T> ScoreMap;
typedef std::map<SCORER_TYPE_T, int> RankMap;
typedef ScoreMap::iterator ScoreMapIterator;
typedef RankMap::iterator RankMapIterator;

#endif //MATCHOBJECTS_H


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
