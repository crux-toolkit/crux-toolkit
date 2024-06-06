#ifndef _CNMZSAMPLELIST_H
#define _CNMZSAMPLELIST_H

#include "NeoMzMLStructs.h"
#include "CnmzSample.h"
#include <vector>


class CnmzSampleList {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzSample> sample;

  int count;

private:

};

#endif