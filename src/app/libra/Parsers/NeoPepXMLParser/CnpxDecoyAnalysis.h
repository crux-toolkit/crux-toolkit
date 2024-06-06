#ifndef _CNPXDECOYANALYSIS_H
#define _CNPXDECOYANALYSIS_H

#include <iostream>
#include <string>
#include <vector>

typedef struct npxPointD{
  double fdr_pp;
  double fdr_pp_decoy;
  double fdr_ip;
  double fdr_ip_decoy;
  double num_corr_pp;
  double num_corr_pp_decoy;
  double num_corr_ip;
  double num_corr_ip_decoy;
  double pp_decoy_uncert;
  double pp_uncert;
  double ip_decoy_uncert;
  double ip_uncert;
  double prob_cutoff;
  npxPointD(){
    fdr_pp=0;
    fdr_pp_decoy = 0;
    fdr_ip = 0;
    fdr_ip_decoy = 0;
    num_corr_pp = 0;
    num_corr_pp_decoy = 0;
    num_corr_ip = 0;
    num_corr_ip_decoy = 0;
    pp_decoy_uncert = 0;
    pp_uncert = 0;
    ip_decoy_uncert = 0;
    ip_uncert = 0;
    prob_cutoff = 0;
  }
  void write(FILE* f){
    fprintf(f,"<point fdr_ip=\"%.16lf\"",fdr_ip);
    fprintf(f," fdr_ip_decoy=\"%.16lf\"",fdr_ip_decoy);
    fprintf(f, " fdr_pp=\"%.16lf\"", fdr_pp);
    fprintf(f, " fdr_pp_decoy=\"%.16lf\"", fdr_pp_decoy);
    fprintf(f, " ip_decoy_uncert=\"%.16lf\"", ip_decoy_uncert);
    fprintf(f, " ip_uncert=\"%.16lf\"", ip_uncert);
    fprintf(f, " num_corr_ip=\"%.16lf\"", num_corr_ip);
    fprintf(f, " num_corr_ip_decoy=\"%.16lf\"", num_corr_ip_decoy);
    fprintf(f, " num_corr_pp_decoy=\"%.16lf\"", num_corr_pp_decoy);
    fprintf(f, " pp_decoy_uncert=\"%.16lf\"", pp_decoy_uncert);
    fprintf(f, " pp_uncert=\"%.16lf\"", pp_uncert);
    fprintf(f, " prob_cutoff=\"%.16lf\"", prob_cutoff);
    fprintf(f, "/>\n");
  }
} npxPointD;

class CnpxDecoyAnalysis {
public:
  void write(FILE* f);

  std::vector<npxPointD> point;

private:

};

#endif
