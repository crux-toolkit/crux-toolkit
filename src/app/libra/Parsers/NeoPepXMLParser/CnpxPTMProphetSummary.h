#ifndef _CNPXPTMPROPHETSUMMARY_H
#define _CNPXPTMPROPHETSUMMARY_H

#include "CnpxInputFile.h"
#include "CnpxMixtureModel.h"
#include "CnpxROCErrorData.h"
#include <string>
#include <vector>

class CnpxPTMProphetSummary {
public:
  CnpxPTMProphetSummary();

  void write(FILE* f);

  std::string frag_ppm_tol;
  std::string min_o;
  std::string min_o_factors;
  std::string mod_string;
  std::string options;
  std::string version;

  std::vector<CnpxInputFile> inputfile;
  std::vector<CnpxMixtureModel> mixturemodel;
  std::vector<CnpxROCErrorData> roc_error_data;

private:

};

#endif