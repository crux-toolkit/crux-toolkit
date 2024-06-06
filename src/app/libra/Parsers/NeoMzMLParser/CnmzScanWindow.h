#ifndef _CNMZSCANWINDOW_H
#define _CNMZSCANWINDOW_H

#include "NeoMzMLStructs.h"
#include "CnmzCvParam.h"

#include <string>
#include <vector>

class CnmzScanWindow {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  //std::vector<CnmzReferenceableParamGroupRef> referencableParamGroupRef;
  std::vector<CnmzCvParam> cvParam;
  //std::vector<CnmzUserParam> userParam;


private:

};

#endif