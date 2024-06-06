#ifndef _CNPRMODIFICATIONINFO_H
#define _CNPRMODIFICATIONINFO_H

#include "NeoProtXMLStructs.h"
#include "CnprModAminoacidMass.h"
#include <string>
#include <vector>

class CnprModificationInfo {
public:

  CnprModificationInfo();

  void write(FILE* f, int tabs = -1);

  std::vector<CnprModAminoacidMass> mod_aminoacid_mass;
  double mod_nterm_mass;
  double mod_cterm_mass;
  std::string modified_peptide;

private:

};

#endif