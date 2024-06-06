#ifndef _CNPXPEPXMLQUANTRESULT_H
#define _CNPXPEPXMLQUANTRESULT_H

#include "CnpxSearchScoreSummary.h"
#include <string>
#include <vector>

class CnpxPepXMLQuantResult {
public:

  CnpxPepXMLQuantResult();
  CnpxPepXMLQuantResult(bool b);

  bool present();
  void write(FILE* f);

  double area;
  double retention_time_sec;

  CnpxSearchScoreSummary search_score_summary;

private:
  bool active;

};

#endif
