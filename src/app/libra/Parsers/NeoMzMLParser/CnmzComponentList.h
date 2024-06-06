#ifndef _CNMZCOMPONENTLIST_H
#define _CNMZCOMPONENTLIST_H

#include "NeoMzMLStructs.h"
#include "CnmzAnalyzer.h"
#include "CnmzDetector.h"
#include "CnmzSource.h"
#include <vector>


class CnmzComponentList {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzSource> source;
  std::vector<CnmzAnalyzer> analyzer;
  std::vector<CnmzDetector> detector;

  int count;

private:

};

#endif