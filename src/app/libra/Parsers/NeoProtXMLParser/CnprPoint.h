#ifndef _CNPRPOINT_H
#define _CNPRPOINT_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprPoint {
public:

  CnprPoint();

  void write(FILE* f, int tabs = -1);

  double fdr_pp;
  double fdr_pp_decoy;
  double num_corr_pp;
  double num_corr_pp_decoy;
  double pp_decoy_uncert;
  double pp_uncert;
  double prob_cutoff;

private:

};

#endif