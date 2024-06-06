#ifndef _CNMZPRECURSORLIST_H
#define _CNMZPRECURSORLIST_H

#include "NeoMzMLStructs.h"
#include "CnmzPrecursor.h"
#include <vector>


class CnmzPrecursorList {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzPrecursor> precursor;

  int count;

private:

};

#endif