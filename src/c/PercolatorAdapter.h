/**
 * \file PercolatorAdapter.h
 * $Revision$
 * \brief Converts Percolator result objects to Crux result objects.
 */

#ifndef PERCOLATORADAPTER_H_
#define PERCOLATORADAPTER_H_

#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "MatchCollection.h"
#include "PeptideMatch.h"
#include "PostProcessProtein.h"
#include "ProteinMatch.h"
#include "SpectrumMatch.h"
#include "src/external/percolator/src/Scores.h"

using namespace std;

class PercolatorAdapter {

public:
  PercolatorAdapter();
  virtual ~PercolatorAdapter();

  static MatchCollection* psmScoresToMatchCollection(Scores* scores);
  static void addPsmScores(ProteinMatchCollection* collection, Scores* scores);
  static void addProteinScores(ProteinMatchCollection* collection, Scores* scores);
  static void addPeptideScores(ProteinMatchCollection* collection, Scores* scores);

protected:
  static int parseChargeState(string psm_id);

};

#endif /* PERCOLATORADAPTER_H_ */


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
