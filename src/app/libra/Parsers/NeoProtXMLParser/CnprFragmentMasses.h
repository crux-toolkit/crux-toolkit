#ifndef _CNPRFRAGMENTMASSES_H
#define _CNPRFRAGMENTMASSES_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprFragmentMasses {
public:

  CnprFragmentMasses();

  void write(FILE* f, int tabs = -1);

  int channel;
  double mz;

private:

};

#endif