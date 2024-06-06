#ifndef _CNMZRUN_H
#define _CNMZRUN_H

#include "NeoMzMLStructs.h"
#include "CnmzSpectrumList.h"


#include <string>
#include <vector>

class CnmzRun {
public:

  void write(FILE* f, int tabs = -1, bool iterative=false);

  //std::vector<CnmzReferenceableParamGroupRef> referencableParamGroupRef;
  //std::vector<CnmzCvParam> cvParam;
  //std::vector<CnmzUserParam> userParam;
  CnmzSpectrumList spectrumList;
  //CnmzChromatogramList chromatogramList;

  std::string defaultInstrumentConfigurationRef;
  std::string defaultSourceFileRef;
  std::string id;
  std::string sampleRef;
  nmzDateTime startTimeStamp;


private:

};

#endif