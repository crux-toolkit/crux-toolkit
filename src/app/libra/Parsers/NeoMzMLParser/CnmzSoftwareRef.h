#ifndef _CNMZSOFTWAREREF_H
#define _CNMZSOFTWAREREF_H

#include "NeoMzMLStructs.h"
#include <string>


class CnmzSoftwareRef {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::string ref;

private:

};

#endif