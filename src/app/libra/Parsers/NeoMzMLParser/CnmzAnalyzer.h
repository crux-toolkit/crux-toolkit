#ifndef _CNMZANALYZER_H
#define _CNMZANALYZER_H

#include "NeoMzMLStructs.h"
#include "CnmzCvParam.h"
#include "CnmzReferenceableParamGroupRef.h"
#include "CnmzUserParam.h"
#include <vector>
#include <string>


class CnmzAnalyzer {
public:
  CnmzAnalyzer(){order=-1000;}
  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzReferenceableParamGroupRef> referenceableParamGroupRef;
  std::vector<CnmzCvParam> cvParam;
  std::vector<CnmzUserParam> userParam;

  int order;

private:

};

#endif