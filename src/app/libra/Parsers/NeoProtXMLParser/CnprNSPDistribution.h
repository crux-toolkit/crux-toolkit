#ifndef _CNPRNSPDISTRIBUTION_H
#define _CNPRNSPDISTRIBUTION_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprNSPDistribution {
public:

  CnprNSPDistribution();

  void write(FILE* f, int tabs = -1);

  int bin_no;
  double nsp_lower_bound_incl;
  std::string nsp_upper_bound_excl;
  double nsp_lower_bound_excl;
  std::string nsp_upper_bound_incl;
  double pos_freq;
  double neg_freq;
  double pos_to_neg_ratio;
  double alt_pos_to_neg_ratio;

private:

};

#endif