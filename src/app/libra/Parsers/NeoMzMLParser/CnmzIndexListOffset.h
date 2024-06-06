#ifndef _CNMZINDEXLISTOFFSET_H
#define _CNMZINDEXLISTOFFSET_H

#include "NeoMzMLStructs.h"

#include <string>
#include <vector>

class CnmzIndexListOffset {
public:

  void write(FILE* f, int tabs = -1, bool iterative=false);

  std::string content;


private:

};

#endif