#ifndef _CNPXPTMPROPHETRESULT_H
#define _CNPXPTMPROPHETRESULT_H

#include "CnpxLability.h"
#include "CnpxModAminoAcidProbability.h"
#include "CnpxModTerminalProbability.h"
#include "CnpxParameter.h"
#include <string>
#include <vector>

class CnpxPTMProphetResult {
public:

  CnpxPTMProphetResult();
  CnpxPTMProphetResult(bool b);

  bool present();
  void write(FILE* f);

  std::string ptm;
  std::string prior;
  std::string ptm_peptide;

  std::vector<CnpxParameter> parameter;
  std::vector<CnpxLability> lability;
  std::vector<CnpxModAminoAcidProbability> mod_amino_acid_probability;
  std::vector<CnpxModTerminalProbability> mod_terminal_probability;

private:
  bool active;

};

#endif
