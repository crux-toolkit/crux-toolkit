#ifndef _CNPRNSPINFORMATION_H
#define _CNPRNSPINFORMATION_H

#include "NeoProtXMLStructs.h"
#include "CnprNSPDistribution.h"
#include <string>
#include <vector>

class CnprNSPInformation {
public:

  void write(FILE* f, int tabs = -1);

  std::string neighboring_bin_smoothing;

  std::vector<CnprNSPDistribution> nsp_distribution;

private:

};

#endif