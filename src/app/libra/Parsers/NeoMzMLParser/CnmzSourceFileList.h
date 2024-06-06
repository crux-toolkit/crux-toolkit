#ifndef _CNMZSOURCEFILELIST_H
#define _CNMZSOURCEFILELIST_H

#include "NeoMzMLStructs.h"
#include "CnmzSourceFile.h"
#include <vector>


class CnmzSourceFileList {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzSourceFile> sourceFile;

  int count;

private:

};

#endif