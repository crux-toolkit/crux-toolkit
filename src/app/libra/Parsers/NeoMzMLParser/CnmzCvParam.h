#ifndef _CNMZCVPARAM_H
#define _CNMZCVPARAM_H

#include "NeoMzMLStructs.h"
#include <string>

class CnmzCvParam {
public:

  void write(FILE* f, int tabs = -1, bool iterative=false);

  std::string accession;
  std::string cvRef;
  std::string name;
  std::string unitAccession;
  std::string unitCvRef;
  std::string unitName;
  std::string value;

private:

};

#endif