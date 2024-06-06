#ifndef _CNPXINTENSITY_H
#define _CNPXINTENSITY_H

#include "NeoPepXMLStructs.h"
#include <string>
#include <vector>

class CnpxIntensity {
public:
  CnpxIntensity();

  void write(FILE* f, int tabs = -1);

  int channel;
  double target_mass;
  double absolute;
  double normalized;
  bool reject;

private:

};

#endif