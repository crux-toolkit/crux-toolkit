#ifndef _CNPXPEPTIDEPROPHETRESULT_H
#define _CNPXPEPTIDEPROPHETRESULT_H

#include "CnpxSearchScoreSummary.h"
#include <string>
#include <vector>

class CnpxPeptideProphetResult {
public:

  CnpxPeptideProphetResult();
  CnpxPeptideProphetResult(bool b);

  bool present();
  void write(FILE* f);

  std::string all_ntt_prob;
  std::string analysis;
  double probability;
  double pep1_probability;
  double pep2_probability;
  
  CnpxSearchScoreSummary search_score_summary;

private:
  bool active;

};

#endif
