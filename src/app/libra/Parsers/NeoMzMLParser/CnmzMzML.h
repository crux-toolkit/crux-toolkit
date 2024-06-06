#ifndef _CNMZMZML_H
#define _CNMZMZML_H

#include "NeoMzMLStructs.h"
#include "CnmzCvList.h"
#include "CnmzDataProcessingList.h"
#include "CnmzFileDescription.h"
#include "CnmzInstrumentConfigurationList.h"
#include "CnmzReferenceableParamGroupList.h"
#include "CnmzRun.h"
#include "CnmzSampleList.h"
#include "CnmzSoftwareList.h"

#include <string>
#include <vector>

class CnmzMzML {
public:

  void write(FILE* f, int tabs = -1, bool iterative=false);

  CnmzCvList cvList;
  CnmzFileDescription fileDescription;
  std::vector<CnmzReferenceableParamGroupList> referencableParamGroupList;
  std::vector<CnmzSampleList> sampleList;
  CnmzSoftwareList softwareList;
  //std::vector<CnmzScanSettingsList> scanSettingsList;
  CnmzInstrumentConfigurationList instrumentConfigurationList;
  CnmzDataProcessingList dataProcessingList;
  CnmzRun run;

  std::string xmlns;
  std::string xmlns_xsi;
  std::string xsi_schemaLocation;
  std::string accession;
  std::string id;
  std::string version;


private:

};

#endif