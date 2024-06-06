#ifndef _CNPXDECOYANALYSISSUMMARY_H
#define _CNPXDECOYANALYSISSUMMARY_H

#include <string>
#include <stdio.h>

class CnpxDecoyAnalysisSummary {
public:
  CnpxDecoyAnalysisSummary();

  void write(FILE* f);

  std::string decoy_string;
  double decoy_ratio;
  std::string exclude_string;
  std::string uniq_iproph_peps;
  std::string uniq_pproph_peps;
  std::string uniq_psm;
  std::string window_prob;

private:

};

#endif
