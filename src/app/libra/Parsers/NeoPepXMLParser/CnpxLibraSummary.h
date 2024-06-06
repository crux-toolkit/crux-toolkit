#ifndef _CNPXLIBRASUMMARY_H
#define _CNPXLIBRASUMMARY_H

#include "NeoPepXMLStructs.h"
#include "CnpxFragmentMasses.h"
#include "CnpxIsotopicContributions.h"
#include <string>
#include <vector>

class CnpxLibraSummary {
public:
  CnpxLibraSummary();

  void write(FILE* f, int tabs = -1);

  double mass_tolerance;
  int centroiding_preference;
  int normalization;
  int output_type;
  std::string channel_code;

  std::vector<CnpxFragmentMasses> fragment_masses;
  std::vector<CnpxIsotopicContributions> isotopic_contributions;

private:

};

#endif