#ifndef _CNPRSTPETERANALYSISSUMMARY_H
#define _CNPRSTPETERANALYSISSUMMARY_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprStPeterAnalysisSummary {
public:

  CnprStPeterAnalysisSummary();

  void write(FILE* f, int tabs = -1);

  std::string version;
  double probability;
  double FDR;
  double tolerance;
  std::string degenerate_peptides;
  double sampleLoad;

private:

};

#endif