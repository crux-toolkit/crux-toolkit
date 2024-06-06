#ifndef _CNPRLIBRARESULT_H
#define _CNPRLIBRARESULT_H

#include "NeoProtXMLStructs.h"
#include "CnprIntensity.h"
#include <string>
#include <vector>

class CnprLibraResult {
public:

  CnprLibraResult();

  void write(FILE* f, int tabs = -1);

  int number;

  std::vector<CnprIntensity> intensity;

private:

};

#endif