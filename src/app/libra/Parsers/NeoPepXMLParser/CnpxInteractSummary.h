#ifndef _CNPXINTERACTSUMMARY_H
#define _CNPXINTERACTSUMMARY_H

#include "CnpxInputFile.h"
#include <string>
#include <vector>

class CnpxInteractSummary {
public:

  void write(FILE* f);

  std::string directory;
  std::string filename;

  std::vector<CnpxInputFile> inputfile;

private:

};

#endif