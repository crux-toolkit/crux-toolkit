#include "CnpxDistributionPoint.h"

using namespace std;

CnpxDistributionPoint::CnpxDistributionPoint() {
  fvalue=0;
  obs_1_distr = 0;
  model_1_pos_distr = 0;
  model_1_neg_distr = 0;
  obs_2_distr = 0;
  model_2_pos_distr = 0;
  model_2_neg_distr = 0;
  obs_3_distr = 0;
  model_3_pos_distr = 0;
  model_3_neg_distr = 0;
  obs_4_distr = 0;
  model_4_pos_distr = 0;
  model_4_neg_distr = 0;
  obs_5_distr = 0;
  model_5_pos_distr = 0;
  model_5_neg_distr = 0;
  obs_6_distr = 0;
  model_6_pos_distr = 0;
  model_6_neg_distr = 0;
  obs_7_distr = 0;
  model_7_pos_distr = 0;
  model_7_neg_distr = 0;
}

void CnpxDistributionPoint::write(FILE* f) {
  fprintf(f, "<distribution_point fvalue=\"%.2lf\"", fvalue);
  fprintf(f, " obs_1_distr=\"%u\"", obs_1_distr);
  fprintf(f, " model_1_pos_distr=\"%.2lf\"", model_1_pos_distr);
  fprintf(f, " model_1_neg_distr=\"%.2lf\"", model_1_neg_distr);
  fprintf(f, " obs_2_distr=\"%u\"", obs_2_distr);
  fprintf(f, " model_2_pos_distr=\"%.2lf\"", model_2_pos_distr);
  fprintf(f, " model_2_neg_distr=\"%.2lf\"", model_2_neg_distr);
  fprintf(f, " obs_3_distr=\"%u\"", obs_3_distr);
  fprintf(f, " model_3_pos_distr=\"%.2lf\"", model_3_pos_distr);
  fprintf(f, " model_3_neg_distr=\"%.2lf\"", model_3_neg_distr);
  fprintf(f, " obs_4_distr=\"%u\"", obs_4_distr);
  fprintf(f, " model_4_pos_distr=\"%.2lf\"", model_4_pos_distr);
  fprintf(f, " model_4_neg_distr=\"%.2lf\"", model_4_neg_distr);
  fprintf(f, " obs_5_distr=\"%u\"", obs_5_distr);
  fprintf(f, " model_5_pos_distr=\"%.2lf\"", model_5_pos_distr);
  fprintf(f, " model_5_neg_distr=\"%.2lf\"", model_5_neg_distr);
  fprintf(f, " obs_6_distr=\"%u\"", obs_6_distr);
  fprintf(f, " model_6_pos_distr=\"%.2lf\"", model_6_pos_distr);
  fprintf(f, " model_6_neg_distr=\"%.2lf\"", model_6_neg_distr);
  fprintf(f, " obs_7_distr=\"%u\"", obs_7_distr);
  fprintf(f, " model_7_pos_distr=\"%.2lf\"", model_7_pos_distr);
  fprintf(f, " model_7_neg_distr=\"%.2lf\"", model_7_neg_distr);
  fprintf(f, "/>\n");
}