#ifndef _CNMZOFFSET_H
#define _CNMZOFFSET_H

#include "NeoMzMLStructs.h"

#include <string>
#include <vector>

class CnmzOffset {
public:

  CnmzOffset(){ scanTime=0;}

  void write(FILE* f, int tabs = -1, bool iterative=false);

  std::string idRef;
  double scanTime;
  std::string spotID;

  std::string content;


private:

};

#endif