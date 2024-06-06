#ifndef _CNMZINDEXLIST_H
#define _CNMZINDEXLIST_H

#include "NeoMzMLStructs.h"
#include "CnmzIndex.h"

#include <string>
#include <vector>

class CnmzIndexList {
public:

  f_off_nmz getIndexListOffset();
  void write(FILE* f, int tabs = -1, bool iterative=false);

  std::vector<CnmzIndex> index;
  int count;


private:

  f_off_nmz fptr;

};

#endif