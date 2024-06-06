#ifndef _CNPXINTERPROPHETRESULT_H
#define _CNPXINTERPROPHETRESULT_H

#include "CnpxSearchScoreSummary.h"
#include <string>
#include <vector>

class CnpxInterprophetResult {
public:

  CnpxInterprophetResult();
  CnpxInterprophetResult(bool b);

  bool present();
  void write(FILE* f);

  std::string all_ntt_prob;
  double probability;

  CnpxSearchScoreSummary search_score_summary;

private:
  bool active;

};

#endif
