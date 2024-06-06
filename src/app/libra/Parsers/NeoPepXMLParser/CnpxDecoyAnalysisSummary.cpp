#include "CnpxDecoyAnalysisSummary.h"

using namespace std;

CnpxDecoyAnalysisSummary::CnpxDecoyAnalysisSummary(){
  decoy_string.clear();
  decoy_ratio=0;
  exclude_string.clear();
  uniq_iproph_peps.clear();
  uniq_pproph_peps.clear();
  uniq_psm.clear();
  window_prob.clear();
}

void CnpxDecoyAnalysisSummary::write(FILE* f){
  fprintf(f, "<decoy_analysis_summary decoy_string=\"%s\"", decoy_string.c_str());
  if(decoy_ratio>0) fprintf(f, " decoy_ratio=\"%.8lf\"", decoy_ratio);
  if (exclude_string.size()>0) fprintf(f, " exclude_string=\"%s\"", exclude_string.c_str());
  if (uniq_iproph_peps.size()>0) fprintf(f, " uniq_iproph_peps=\"%s\"", uniq_iproph_peps.c_str());
  if (uniq_pproph_peps.size()>0) fprintf(f, " uniq_pproph_peps=\"%s\"", uniq_pproph_peps.c_str());
  if (uniq_psm.size()>0) fprintf(f, " uniq_psm=\"%s\"", uniq_psm.c_str());
  if (window_prob.size()>0) fprintf(f, " window_prob=\"%s\"", window_prob.c_str());
  fprintf(f, "/>\n");
}