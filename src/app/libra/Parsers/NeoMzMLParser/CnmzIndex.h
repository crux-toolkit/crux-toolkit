#ifndef _CNMZINDEX_H
#define _CNMZINDEX_H

#include "NeoMzMLStructs.h"
#include "CnmzOffset.h"

#include <string>
#include <vector>

class CnmzIndex {
public:

  void write(FILE* f, int tabs = -1, bool iterative=false);

  std::vector<CnmzOffset> offset;
  std::string name;


private:

};

#endif