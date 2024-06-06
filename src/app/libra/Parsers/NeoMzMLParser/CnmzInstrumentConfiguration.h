#ifndef _CNMZINSTRUMENTCONFIGURATION_H
#define _CNMZINSTRUMENTCONFIGURATION_H

#include "NeoMzMLStructs.h"
#include "CnmzComponentList.h"
#include "CnmzCvParam.h"
#include "CnmzReferenceableParamGroupRef.h"
#include "CnmzSoftwareRef.h"
#include "CnmzUserParam.h"
#include <vector>
#include <string>


class CnmzInstrumentConfiguration {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzReferenceableParamGroupRef> referenceableParamGroupRef;
  std::vector<CnmzCvParam> cvParam;
  std::vector<CnmzUserParam> userParam;
  std::vector<CnmzComponentList> componentList;
  std::vector<CnmzSoftwareRef> softwareRef;

  std::string id;
  std::string scanSettingsRef;

private:

};

#endif