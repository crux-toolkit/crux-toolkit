#ifndef _CNMZSELECTEDIONLIST_H
#define _CNMZSELECTEDIONLIST_H

#include "NeoMzMLStructs.h"
#include "CnmzSelectedIon.h"
#include <vector>


class CnmzSelectedIonList {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzSelectedIon> selectedIon;

  int count;

private:

};

#endif