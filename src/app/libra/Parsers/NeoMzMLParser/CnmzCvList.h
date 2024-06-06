#ifndef _CNMZCVLIST_H
#define _CNMZCVLIST_H

#include "NeoMzMLStructs.h"
#include "CnmzCv.h"

#include <string>
#include <vector>

class CnmzCvList{
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzCv> cv;

  int count;


private:

};

#endif