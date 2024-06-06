#include "CnprPoint.h"

using namespace std;

CnprPoint::CnprPoint(){
  fdr_pp=-1;
  fdr_pp_decoy = -1;
  num_corr_pp = -1;
  num_corr_pp_decoy = -1;
  pp_decoy_uncert = -1;
  pp_uncert = -1;
  prob_cutoff = -1;
}

void CnprPoint::write(FILE* f, int tabs){
  //no requirements

  NPRprintTabs(f, tabs);
  fprintf(f, "<point");
  if (fdr_pp>-1) fprintf(f, " fdr_pp=\"%.16lf\"", fdr_pp);
  if (fdr_pp_decoy>-1) fprintf(f, " fdr_pp_decoy=\"%.16lf\"", fdr_pp_decoy);
  if (num_corr_pp>-1) fprintf(f, " num_corr_pp=\"%.3lf\"", num_corr_pp);
  if (num_corr_pp_decoy>-1) fprintf(f, " num_corr_pp_decoy=\"%.16lf\"", num_corr_pp_decoy);
  if (pp_decoy_uncert>-1) fprintf(f, " pp_decoy_uncert=\"%.16lf\"", pp_decoy_uncert);
  if (pp_uncert>-1) fprintf(f, " pp_uncert=\"%.16lf\"", pp_uncert);
  if (prob_cutoff>-1) fprintf(f, " prob_cutoff=\"%.3lf\"", prob_cutoff);
  fprintf(f, "/>\n");
}

