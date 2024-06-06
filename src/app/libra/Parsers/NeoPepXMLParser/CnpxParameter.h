#ifndef _CNPXPARAMETER_H
#define _CNPXPARAMETER_H

#include "NeoPepXMLStructs.h"
#include <string>

class CnpxParameter {
public:

  void write(FILE* f, int tabs=-1);

  std::string name;
  std::string type;
  std::string value;

private:

};

#endif
