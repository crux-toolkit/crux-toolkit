#ifndef _CNMZREFERENCEABLEPARAMGROUP_H
#define _CNMZREFERENCEABLEPARAMGROUP_H

#include "NeoMzMLStructs.h"
#include "CnmzCvParam.h"
#include "CnmzReferenceableParamGroupRef.h"
#include "CnmzUserParam.h"
#include <vector>
#include <string>


class CnmzReferenceableParamGroup {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzCvParam> cvParam;
  std::vector<CnmzUserParam> userParam;

  std::string id;

private:

};

#endif