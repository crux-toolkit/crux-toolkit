#ifndef _CNPRINDISTINGUISHABLEPEPTIDE_H
#define _CNPRINDISTINGUISHABLEPEPTIDE_H

#include "NeoProtXMLStructs.h"
#include "CnprModificationInfo.h"
#include <string>
#include <vector>

class CnprIndistinguishablePeptide {
public:

  CnprIndistinguishablePeptide();

  void write(FILE* f, int tabs = -1);

  std::vector<CnprModificationInfo> modification_info;
  std::string peptide_sequence;
  int charge;
  double calc_neutral_pep_mass;

private:

};

#endif