#ifndef _CNPXMIXTUREMODEL_H
#define _CNPXMIXTUREMODEL_H

#include <iostream>
#include <string>
#include <vector>

typedef struct npxPointM {
  float neg_dens;
  float neg_obs_dens;
  float pos_dens;
  float pos_obs_dens;
  float value;
  npxPointM(){
    neg_dens=0;
    neg_obs_dens = 0;
    pos_dens = 0;
    pos_obs_dens = 0;
    value = 0;
  }
  void write(FILE* f){
    fprintf(f, "<point value=\"%.6lf\"", value);
    fprintf(f, " pos_dens=\"%.6f\"", pos_dens);
    fprintf(f, " neg_dens=\"%.6f\"", neg_dens);
    if(neg_obs_dens>0) fprintf(f, " neg_obs_dens=\"%.6f\"", neg_obs_dens);
    if(pos_obs_dens>0) fprintf(f, " pos_obs_dens=\"%.6f\"", pos_obs_dens);
    fprintf(f, "/>\n");
  }
} npxPointM;

typedef struct npxBin{
  bool value;
  double pos_prob;
  double neg_prob;
  npxBin(){
    value=false;
    pos_prob=0;
    neg_prob=0;
  }
  void write(FILE* f){
    if(value) fprintf(f,"<bin value=\"true\"");
    else fprintf(f, "<bin value=\"false\"");
    fprintf(f, " pos_prob=\"%.6lf\"", pos_prob);
    fprintf(f, " neg_prob=\"%.6lf\"", neg_prob);
    fprintf(f, "/>\n");
  }
} npxBin;


class CnpxMixtureModel {
public:
  CnpxMixtureModel();

  void write(FILE* f);

  std::string name;
  float neg_bandwidth;
  float pos_bandwidth;

  std::vector<npxPointM> point;
  std::vector<npxBin> bin;

private:

};

#endif
