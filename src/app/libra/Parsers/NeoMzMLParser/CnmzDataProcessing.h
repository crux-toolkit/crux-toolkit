#ifndef _CNMZDATAPROCESSING_H
#define _CNMZDATAPROCESSING_H

#include "NeoMzMLStructs.h"
#include "CnmzProcessingMethod.h"
#include <vector>
#include <string>


class CnmzDataProcessing {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzProcessingMethod> processingMethod;

  std::string id;

private:

};

#endif