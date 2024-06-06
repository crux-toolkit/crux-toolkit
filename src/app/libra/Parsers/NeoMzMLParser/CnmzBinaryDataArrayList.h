#ifndef _CNMZBINARYDATAARRAYLIST_H
#define _CNMZBINARYDATAARRAYLIST_H

#include "NeoMzMLStructs.h"
#include "CnmzBinaryDataArray.h"

#include <string>
#include <vector>

class CnmzBinaryDataArrayList {
public:

  void write(FILE* f, int tabs = -1, bool interative=false);

  std::vector<CnmzBinaryDataArray> binaryDataArray;

  int count;


private:

};

#endif