#ifndef _CNMZBINARYDATAARRAY_H
#define _CNMZBINARYDATAARRAY_H

#include "NeoMzMLStructs.h"
#include "CnmzBinary.h"
#include "CnmzCvParam.h"

#include <string>
#include <vector>

class CnmzBinaryDataArray {
public:

  void write(FILE* f, int tabs = -1, bool iterative=false);

  //std::vector<CnmzReferenceableParamGroupRef> referencableParamGroupRef;
  std::vector<CnmzCvParam> cvParam;
  //std::vector<CnmzUserParam> userParam;
  CnmzBinary binary;

  int arrayLength;
  std::string dataProcessingRef;
  int encodedLength;

private:

};

#endif