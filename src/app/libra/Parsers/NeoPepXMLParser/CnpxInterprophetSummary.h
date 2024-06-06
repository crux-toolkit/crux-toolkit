#ifndef _CNPXINTERPROPHETSUMMARY_H
#define _CNPXINTERPROPHETSUMMARY_H

#include "CnpxInputFile.h"
#include "CnpxMixtureModel.h"
#include "CnpxMixtureModelDistribution.h"
#include "CnpxROCErrorData.h"
#include <string>
#include <vector>

class CnpxInterprophetSummary {
public:
  CnpxInterprophetSummary();

  void write(FILE* f);

  double est_tot_num_correct_pep;
  double est_tot_num_correct_psm;
  std::string options;
  std::string version;  

  std::vector<CnpxInputFile> inputfile;
  std::vector<CnpxMixtureModel> mixturemodel;
  std::vector<CnpxROCErrorData> roc_error_data;
  std::vector<CnpxMixtureModelDistribution> mixturemodel_distribution;

private:

};

#endif