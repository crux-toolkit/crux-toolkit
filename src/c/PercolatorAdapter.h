/**
 * \file PercolatorAdapter.h
 * $Revision$
 * \brief Converts Percolator result objects to Crux result objects.
 */

#ifndef PERCOLATORADAPTER_H_
#define PERCOLATORADAPTER_H_

#include <stdlib.h>
#include <stdio.h>

#include "MatchCollection.h"
#include "PostProcessProtein.h"
#include "src/external/percolator/src/Scores.h"

using namespace std;

class PercolatorAdapter {
public:
  PercolatorAdapter();
  PercolatorAdapter(Scores* scores);
  virtual ~PercolatorAdapter();

  void setScores(Scores* scores);

  MatchCollection* convertFromPsms();
  static MatchCollection* convertFromPsms(Scores* scores);
  MatchCollection* convertFromProteins();
  static MatchCollection* convertFromProteins(Scores* scores);
  MatchCollection* convertFromPeptides();
  static MatchCollection* convertFromPeptides(Scores* scores);

protected:
  Scores* scores_;

  int parseChargeState(string psm_id);
};

#endif /* PERCOLATORADAPTER_H_ */


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
