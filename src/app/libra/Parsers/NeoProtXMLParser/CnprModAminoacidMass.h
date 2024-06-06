#ifndef _CNPRMODAMINOACIDMASS_H
#define _CNPRMODAMINOACIDMASS_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprModAminoacidMass {
public:
  CnprModAminoacidMass();

  void write(FILE* f, int tabs = -1);

  int position;
  double mass;

private:

};

#endif