#ifndef _CNPXSAMPLEENZYME_H
#define _CNPXSAMPLEENZYME_H

#include "CnpxSpecificity.h"
#include "NeoPepXMLStructs.h"
#include <string>
#include <vector>

class CnpxSampleEnzyme {
public:
  CnpxSampleEnzyme();

  void write(FILE* f, int tabs=-1);

  std::string description;
  std::string fidelity;
  bool independent;
  std::string name;

  std::vector<CnpxSpecificity> specificity;

private:

};

#endif
