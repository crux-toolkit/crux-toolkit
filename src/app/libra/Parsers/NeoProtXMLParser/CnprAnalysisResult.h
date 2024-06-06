#ifndef _CNPRANALYSISRESULT_H
#define _CNPRANALYSISRESULT_H

#include "NeoProtXMLStructs.h"
#include "CnprLibraResult.h"
#include "CnprStPeterQuant.h"
#include <string>
#include <vector>

class CnprAnalysisResult {
public:

  CnprAnalysisResult();

  void write(FILE* f, int tabs = -1);

  std::string analysis;
  int id;

  std::vector<CnprLibraResult> libra_result;
  std::vector<CnprStPeterQuant> StPeterQuant;

private:

};

#endif