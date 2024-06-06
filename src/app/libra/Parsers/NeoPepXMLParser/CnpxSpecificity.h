#ifndef _CNPXSPECIFICITY_H
#define _CNPXSPECIFICITY_H

#include "NeoPepXMLStructs.h"
#include <string>

class CnpxSpecificity {
public:
  CnpxSpecificity();

  void write(FILE* f, int tabs=-1);

  std::string cut;
  unsigned int min_spacing;
  std::string no_cut;
  std::string sense;

private:

};

#endif
