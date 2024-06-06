#ifndef _CNPXPEPTIDEPROPHETSUMMARY_H
#define _CNPXPEPTIDEPROPHETSUMMARY_H

#include "CnpxDistributionPoint.h"
#include "CnpxInputFile.h"
#include "CnpxMixtureModel.h"
#include "CnpxMixture_Model.h"
#include "CnpxROCErrorData.h"
#include <string>
#include <vector>

class CnpxPeptideprophetSummary {
public:
  CnpxPeptideprophetSummary();

  void write(FILE* f);

  std::string author;
  double est_tot_num_correct;
  double min_prob;
  std::string options;
  std::string type;
  std::string version;

  std::vector<CnpxDistributionPoint> distribution_point;
  std::vector<CnpxInputFile> inputfile;
  std::vector<CnpxMixture_Model> mixture_model;
  std::vector<CnpxROCErrorData> roc_error_data;

private:

};

#endif
