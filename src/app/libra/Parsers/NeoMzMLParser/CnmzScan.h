#ifndef _CNMZSCAN_H
#define _CNMZSCAN_H

#include "NeoMzMLStructs.h"
#include "CnmzReferenceableParamGroupRef.h"
#include "CnmzCvParam.h"
#include "CnmzScanWindowList.h"
#include "CnmzUserParam.h"

#include <string>
#include <vector>

class CnmzScan {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzReferenceableParamGroupRef> referencableParamGroupRef;
  std::vector<CnmzCvParam> cvParam;
  std::vector<CnmzUserParam> userParam;
  std::vector<CnmzScanWindowList> scanWindowList;

  std::string externalSpectrumID;
  std::string instrumentConfigurationRef;
  std::string sourceFileRef;
  std::string spectrumRef;


private:

};

#endif