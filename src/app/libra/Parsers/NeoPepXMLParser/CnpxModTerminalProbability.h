#ifndef _CNPXMODTERMINALPROBABILITY_H
#define _CNPXMODTERMINALPROBABILITY_H

#include <string>

class CnpxModTerminalProbability {
public:

  CnpxModTerminalProbability();

  void write(FILE* f);

  double probability;
  char terminus;
  double oscore;
  double mscore;
  double direct_oscore;
  double direct_mscore;
  double cterm_score;
  double nterm_score;
  char shift;

private:

};

#endif
