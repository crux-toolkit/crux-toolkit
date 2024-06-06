#ifndef _CNMZFILEDESCRIPTION_H
#define _CNMZFILEDESCRIPTION_H

#include "NeoMzMLStructs.h"
#include "CnmzContact.h"
#include "CnmzFileContent.h"
#include "CnmzSourceFileList.h"

#include <string>
#include <vector>

class CnmzFileDescription{
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  CnmzFileContent fileContent;
  std::vector<CnmzSourceFileList> sourceFileList;
  std::vector<CnmzContact> contact;

private:

};

#endif