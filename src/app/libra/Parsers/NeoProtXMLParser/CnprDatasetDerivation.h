#ifndef _CNPRDATASETDERIVATION_H
#define _CNPRDATASETDERIVATION_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprDatasetDerivation {
public:

  void write(FILE* f, int tabs = -1);

  std::string generation_no;
  //std::vector<CnprDataFilter> data_filter;

private:

};

#endif