#ifndef _CNMZPROCESSINGMETHOD_H
#define _CNMZPROCESSINGMETHOD_H

#include "NeoMzMLStructs.h"
#include "CnmzCvParam.h"
#include "CnmzReferenceableParamGroupRef.h"
#include "CnmzUserParam.h"
#include <vector>
#include <string>


class CnmzProcessingMethod {
public:
  CnmzProcessingMethod(){ order = -1000; }
  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzReferenceableParamGroupRef> referenceableParamGroupRef;
  std::vector<CnmzCvParam> cvParam;
  std::vector<CnmzUserParam> userParam;

  int order;
  std::string softwareRef;

private:

};

#endif