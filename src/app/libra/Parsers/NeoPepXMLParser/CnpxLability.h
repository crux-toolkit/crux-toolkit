#ifndef _CNPXLABILITY_H
#define _CNPXLABILITY_H

#include <string>

class CnpxLability {
public:

  CnpxLability();

  void write(FILE* f);

  int numlosses;
  double pval;
  double probability;
  double oscore;
  double mscore;
  double cterm_score;
  double nterm_score;

private:

};

#endif
