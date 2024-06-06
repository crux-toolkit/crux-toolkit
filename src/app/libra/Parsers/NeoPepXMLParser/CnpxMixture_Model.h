#ifndef _CNPXMIXTURE_MODEL_H
#define _CNPXMIXTURE_MODEL_H

#include "CnpxMixtureModel.h"
#include "CnpxMixtureModelDistribution.h"
#include <iostream>
#include <string>
#include <vector>

class CnpxMixture_Model {
public:
  CnpxMixture_Model();

  void write(FILE* f);

  std::string comments;
  double est_tot_correct;
  int num_iterations;
  int precursor_ion_charge;
  double prior_probability;
  int tot_num_spectra;
  
  std::vector<CnpxMixtureModelDistribution> mixturemodel_distribution;
  std::vector<CnpxMixtureModel> mixturemodel;

private:

};

#endif

