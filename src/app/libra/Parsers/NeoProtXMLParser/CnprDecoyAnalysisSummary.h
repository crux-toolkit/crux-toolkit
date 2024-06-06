#ifndef _CNPRDECOYANALYSISSUMMARY_H
#define _CNPRDECOYANALYSISSUMMARY_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprDecoyAnalysisSummary {
public:

  CnprDecoyAnalysisSummary();

  void write(FILE* f, int tabs = -1);

  std::string decoy_string;
  double decoy_ratio;
  std::string exclude_string;
  std::string use_confidence;

private:

};

#endif