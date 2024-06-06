#ifndef _CNPRPDECOYANALYSIS_H
#define _CNPRPDECOYANALYSIS_H

#include "NeoProtXMLStructs.h"
#include "CnprPoint.h"
#include <string>
#include <vector>

class CnprDecoyAnalysis {
public:

  void write(FILE* f, int tabs = -1);

  std::vector<CnprPoint> point;

private:

};

#endif