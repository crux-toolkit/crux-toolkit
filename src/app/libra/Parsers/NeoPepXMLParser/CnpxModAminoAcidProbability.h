#ifndef _CNPXMODAMINOACIDPROBABILITY_H
#define _CNPXMODAMINOACIDPROBABILITY_H

#include <string>

class CnpxModAminoAcidProbability {
public:

  CnpxModAminoAcidProbability();

  void write(FILE* f);

  int position;
  double probability;
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
