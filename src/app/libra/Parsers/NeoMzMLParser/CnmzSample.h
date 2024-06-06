#ifndef _CNMZSAMPLE_H
#define _CNMZSAMPLE_H

#include "NeoMzMLStructs.h"
#include "CnmzCvParam.h"
#include "CnmzReferenceableParamGroupRef.h"
#include "CnmzUserParam.h"
#include <vector>
#include <string>


class CnmzSample {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzReferenceableParamGroupRef> referenceableParamGroupRef;
  std::vector<CnmzCvParam> cvParam;
  std::vector<CnmzUserParam> userParam;

  std::string id;
  std::string name;

private:

};

#endif