#ifndef _CNMZSPECTRUM_H
#define _CNMZSPECTRUM_H

#include "NeoMzMLStructs.h"
#include "CnmzBinaryDataArrayList.h"
#include "CnmzPrecursorList.h"
#include "CnmzReferenceableParamGroupRef.h"
#include "CnmzScanList.h"
#include "CnmzUserParam.h"

#include <string>
#include <vector>

class CnmzSpectrum {
public:

  void write(FILE* f, int tabs = -1, bool iterative=false);

  std::vector<CnmzReferenceableParamGroupRef> referencableParamGroupRef;
  std::vector<CnmzCvParam> cvParam;
  std::vector<CnmzUserParam> userParam;
  std::vector<CnmzScanList> scanList;
  std::vector<CnmzPrecursorList> precursorList;
  //std::vector<CnmzProductList> productList;
  std::vector<CnmzBinaryDataArrayList> binaryDataArrayList;

  std::string dataProcessingRef;
  int defaultArrayLength;
  std::string id;
  int index;
  std::string sourceFileRef;
  std::string spotID;


private:

};

#endif