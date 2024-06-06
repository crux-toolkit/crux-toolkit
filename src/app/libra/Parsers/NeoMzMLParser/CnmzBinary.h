#ifndef _CNMZBINARY_H
#define _CNMZBINARY_H

#include "NeoMzMLStructs.h"

#include <string>
#include <vector>

class CnmzBinary {
public:

  void write(FILE* f, int tabs = -1, bool iterative=false);

  std::string content;


private:

};

#endif