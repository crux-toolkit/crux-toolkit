#ifndef _CNMZSOFTWARELIST_H
#define _CNMZSOFTWARELIST_H

#include "NeoMzMLStructs.h"
#include "CnmzSoftware.h"
#include <vector>


class CnmzSoftwareList {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzSoftware> software;

  int count;

private:

};

#endif