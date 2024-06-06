#ifndef _CNPRPROTEINGROUP_H
#define _CNPRPROTEINGROUP_H

#include "NeoProtXMLStructs.h"
#include "CnprProtein.h"
#include <string>
#include <vector>

class CnprProteinGroup {
public:

  CnprProteinGroup();

  void write(FILE* f, int tabs = -1);

  std::string group_number;
  std::string pseudo_name;
  double probability;

  std::vector<CnprProtein> protein;

private:

};

#endif