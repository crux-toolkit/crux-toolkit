#ifndef MATCHOBJECTS_H_
#define MATCHOBJECTS_H_

#include "objects.h"
#include <vector>
#include <map>

class ProteinMatchCollection;
class ProteinMatch;
class PeptideMatch;
class SpectrumMatch;
class AbstractMatch;

typedef std::vector<SpectrumMatch*> SpectrumMatchCollection;
typedef SpectrumMatchCollection::iterator SpectrumMatchIterator;

typedef std::vector<PeptideMatch*> PeptideMatchCollection;
typedef PeptideMatchCollection::iterator PeptideMatchIterator;

typedef std::vector<ProteinMatch*>::iterator ProteinMatchIterator;

typedef std::map<SCORER_TYPE_T, FLOAT_T> ScoreMap;
typedef ScoreMap::iterator ScoreMapIterator;

#endif //MATCHOBJECTS_H


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
