#ifndef _CNPXLIBRARESULT_H
#define _CNPXLIBRARESULT_H

#include "NeoPepXMLStructs.h"
#include "CnpxIntensity.h"
#include <string>
#include <vector>

class CnpxLibraResult {
public:
  CnpxLibraResult();
  CnpxLibraResult(bool b);

  bool present();
  void write(FILE* f, int tabs = -1);

  bool is_rejected;

  std::vector<CnpxIntensity> intensity;

private:
  bool active;
};

#endif