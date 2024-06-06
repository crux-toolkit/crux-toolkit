#ifndef _CNMZDATAPROCESSINGLIST_H
#define _CNMZDATAPROCESSINGLIST_H

#include "NeoMzMLStructs.h"
#include "CnmzDataProcessing.h"
#include <vector>


class CnmzDataProcessingList {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzDataProcessing> dataProcessing;

  int count;

private:

};

#endif