#include "CnprDecoyAnalysisSummary.h"

using namespace std;

CnprDecoyAnalysisSummary::CnprDecoyAnalysisSummary(){
  decoy_ratio=-1;
}

void CnprDecoyAnalysisSummary::write(FILE* f, int tabs){
  //required
  string el = "decoy_analysis_summary";
  //if (decoy_string.empty()) NPRerrMsg(el, "decoy_string");  //required, but allowed to be empty
  if (decoy_ratio<0) NPRerrMsg(el, "decoy_ratio");

  NPRprintTabs(f, tabs);
  fprintf(f, "<decoy_analysis_summary");
  fprintf(f, " decoy_string=\"%s\"", decoy_string.c_str());
  fprintf(f, " decoy_ratio=\"%.16lf\"", decoy_ratio);
  if (!exclude_string.empty()) fprintf(f, " exclude_string=\"%s\"", exclude_string.c_str());
  if (!use_confidence.empty()) fprintf(f, " use_confidence=\"%s\"", use_confidence.c_str());
  fprintf(f, "/>\n");

}
