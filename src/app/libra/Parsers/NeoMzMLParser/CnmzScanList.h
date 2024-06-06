#ifndef _CNMZSCANLIST_H
#define _CNMZSCANLIST_H

#include "NeoMzMLStructs.h"
#include "CnmzCvParam.h"
#include "CnmzScan.h"

#include <string>
#include <vector>

class CnmzScanList {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  //std::vector<CnmzReferenceableParamGroupRef> referencableParamGroupRef;
  std::vector<CnmzCvParam> cvParam;
  //std::vector<CnmzUserParam> userParam;
  std::vector<CnmzScan> scan;

  int count;


private:

};

#endif