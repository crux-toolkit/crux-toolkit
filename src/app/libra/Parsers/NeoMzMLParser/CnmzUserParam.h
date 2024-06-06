#ifndef _CNMZUSERPARAM_H
#define _CNMZUSERPARAM_H

#include "NeoMzMLStructs.h"
#include <string>


class CnmzUserParam {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::string name;
  std::string type;
  std::string unitAccession;
  std::string unitCvRef;
  std::string unitName;
  std::string value;

private:

};

#endif