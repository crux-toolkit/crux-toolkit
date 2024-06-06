#ifndef _CNMZREFERENCEABLEPARAMGROUPREF_H
#define _CNMZREFERENCEABLEPARAMGROUPREF_H

#include "NeoMzMLStructs.h"
#include <string>


class CnmzReferenceableParamGroupRef {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::string ref;

private:

};

#endif