#ifndef _CNMZREFERENCEABLEPARAMGROUPLIST_H
#define _CNMZREFERENCEABLEPARAMGROUPLIST_H

#include "NeoMzMLStructs.h"
#include "CnmzReferenceableParamGroup.h"
#include <vector>


class CnmzReferenceableParamGroupList {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzReferenceableParamGroup> referenceableParamGroup;

  int count;

private:

};

#endif