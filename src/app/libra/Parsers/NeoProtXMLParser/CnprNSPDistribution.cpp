#include "CnprNSPDistribution.h"

using namespace std;

CnprNSPDistribution::CnprNSPDistribution(){
  bin_no=-1;
  nsp_lower_bound_incl=-1;
  nsp_lower_bound_excl=-1;
  pos_freq=-1;
  neg_freq=-1;
  pos_to_neg_ratio=-1;
  alt_pos_to_neg_ratio=-1;
}

void CnprNSPDistribution::write(FILE* f, int tabs){
  //required
  string el = "nsp_distribution";
  if (bin_no<0) NPRerrMsg(el, "bin_no");
  if (pos_freq<0) NPRerrMsg(el, "pos_freq");
  if (neg_freq<0) NPRerrMsg(el, "neg_freq");
  if (pos_to_neg_ratio<0) NPRerrMsg(el, "pos_to_neg_ratio");

  NPRprintTabs(f, tabs);
  fprintf(f, "<nsp_distribution");
  fprintf(f, " bin_no=\"%d\"", bin_no);
  if (nsp_lower_bound_incl>-1) fprintf(f, " nsp_lower_bound_incl=\"%.2lf\"", nsp_lower_bound_incl);
  if (!nsp_upper_bound_excl.empty()) fprintf(f, " nsp_upper_bound_excl=\"%s\"", nsp_upper_bound_excl.c_str());
  if (nsp_lower_bound_excl>-1) fprintf(f, " nsp_lower_bound_excl=\"%.2lf\"", nsp_lower_bound_excl);
  if (!nsp_upper_bound_incl.empty()) fprintf(f, " nsp_upper_bound_incl=\"%s\"", nsp_upper_bound_incl.c_str());
  fprintf(f, " pos_freq=\"%.2lf\"", pos_freq);
  fprintf(f, " neg_freq=\"%.2lf\"", neg_freq);
  fprintf(f, " pos_to_neg_ratio=\"%.2lf\"", pos_to_neg_ratio);
  if (alt_pos_to_neg_ratio>-1) fprintf(f, " alt_pos_to_neg_ratio=\"%.2lf\"", alt_pos_to_neg_ratio);
  fprintf(f, "/>\n");

}
