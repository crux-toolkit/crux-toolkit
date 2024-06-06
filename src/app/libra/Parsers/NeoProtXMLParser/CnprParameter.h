#ifndef _CNPRPARAMETER_H
#define _CNPRPARAMETER_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprParameter {
public:

  void write(FILE* f, int tabs = -1);

  std::string name;
  std::string value;
  std::string type;

private:

};

#endif