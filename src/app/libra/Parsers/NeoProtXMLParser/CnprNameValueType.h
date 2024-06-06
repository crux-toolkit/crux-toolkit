#ifndef _CNPRNAMEVALUETYPE_H
#define _CNPRNAMEVALUETYPE_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprNameValueType {
public:

  void write(FILE* f, int tabs = -1);

  std::string name;
  std::string value;
  std::string type;

private:

};

#endif