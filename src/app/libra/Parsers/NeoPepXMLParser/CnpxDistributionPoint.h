#ifndef _CNPXDISTRIBUTIONPOINT_H
#define _CNPXDISTRIBUTIONPOINT_H

#include <iostream>

class CnpxDistributionPoint {
public:

  CnpxDistributionPoint();

  void write(FILE* f);

  double fvalue;
  unsigned int obs_1_distr;
  double model_1_pos_distr;
  double model_1_neg_distr;
  unsigned int obs_2_distr;
  double model_2_pos_distr;
  double model_2_neg_distr;
  unsigned int obs_3_distr;
  double model_3_pos_distr;
  double model_3_neg_distr;
  unsigned int obs_4_distr;
  double model_4_pos_distr;
  double model_4_neg_distr;
  unsigned int obs_5_distr;
  double model_5_pos_distr;
  double model_5_neg_distr;
  unsigned int obs_6_distr;
  double model_6_pos_distr;
  double model_6_neg_distr;
  unsigned int obs_7_distr;
  double model_7_pos_distr;
  double model_7_neg_distr;


private:

};

#endif