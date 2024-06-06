#ifndef _CNPRINTENSITY_H
#define _CNPRINTENSITY_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprIntensity {
public:

  CnprIntensity();

  void write(FILE* f, int tabs = -1);

  int channel;
  double mz;
  double ratio;
  double error;
  
private:

};

#endif