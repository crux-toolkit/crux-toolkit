#ifndef _CNPRERRORPOINT_H
#define _CNPRERRORPOINT_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprErrorPoint {
public:

  CnprErrorPoint();

  void write(FILE* f, int tabs = -1);

  double error;
  double min_prob;
  int num_corr;
  int num_incorr;

private:

};

#endif