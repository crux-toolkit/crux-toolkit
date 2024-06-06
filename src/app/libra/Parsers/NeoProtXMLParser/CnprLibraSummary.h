#ifndef _CNPRLIBRASUMMARY_H
#define _CNPRLIBRASUMMARY_H

#include "NeoProtXMLStructs.h"
#include "CnprFragmentMasses.h"
#include "CnprIsotopicContributions.h"
#include <string>
#include <vector>

class CnprLibraSummary {
public:

  CnprLibraSummary();

  void write(FILE* f, int tabs = -1);

  std::string version;
  double mass_tolerance;
  int centroiding_preference;
  int normalization;
  int output_type;
  std::string channel_code;
  double min_pep_prob;
  double min_pep_wt;
  double min_prot_prob;

  std::vector<CnprFragmentMasses> fragment_masses;
  std::vector<CnprIsotopicContributions> isotopic_contributions;

private:

};

#endif