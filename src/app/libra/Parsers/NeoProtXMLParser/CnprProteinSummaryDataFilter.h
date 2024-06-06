#ifndef _CNPRPROTEINSUMMARYDATAFILTER_H
#define _CNPRPROTEINSUMMARYDATAFILTER_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprProteinSummaryDataFilter {
public:

  CnprProteinSummaryDataFilter();

  void write(FILE* f, int tabs = -1);

  double min_probability;
  double sensitivity;
  double false_positive_error_rate;
  double predicted_num_correct;
  double predicted_num_incorrect;

private:

};

#endif