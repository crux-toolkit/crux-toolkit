#ifndef _CNPXSPECTRUMQUERY_H
#define _CNPXSPECTRUMQUERY_H

#include "CnpxSearchResult.h"
#include "NeoPepXMLStructs.h"
#include <iostream>
#include <string>
#include <vector>

class CnpxSpectrumQuery {
public:
  CnpxSpectrumQuery();
  
  CnpxSearchResult* addSearchResult();
  void write(FILE* f, int tabs=-1);

  std::string spectrum;
  std::string spectrumNativeID;
  int start_scan;
  int end_scan;
  double retention_time_sec;
  double collision_energy;
  double compensation_voltage;
  double precursor_intensity;
  std::string activation_method;
  double precursor_neutral_mass;
  int assumed_charge;
  std::string search_specification;
  int index;

  std::vector<CnpxSearchResult> search_result;

private:

};

#endif 
