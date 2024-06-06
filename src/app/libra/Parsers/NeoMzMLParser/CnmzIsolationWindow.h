#ifndef _CNMZISOLATIONWINDOW_H
#define _CNMZISOLATIONWINDOW_H

#include "NeoMzMLStructs.h"
#include "CnmzCvParam.h"
#include "CnmzReferenceableParamGroupRef.h"
#include "CnmzUserParam.h"
#include <vector>
#include <string>


class CnmzIsolationWindow {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzReferenceableParamGroupRef> referenceableParamGroupRef;
  std::vector<CnmzCvParam> cvParam;
  std::vector<CnmzUserParam> userParam;

private:

};

#endif