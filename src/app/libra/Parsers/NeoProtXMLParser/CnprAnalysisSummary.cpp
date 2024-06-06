#include "CnprAnalysisSummary.h"

using namespace std;

CnprAnalysisSummary::CnprAnalysisSummary(){
  time.date.year=0;
  id=-1;
}

void CnprAnalysisSummary::write(FILE* f, int tabs){
  //required
  string el = "analysis_summary";
  if (analysis.empty()) NPRerrMsg(el, "analysis");
  if (time.date.year==0) NPRerrMsg(el, "time");
  if (id<0) NPRerrMsg(el, "id");

  NPRprintTabs(f, tabs);
  fprintf(f, "<analysis_summary");
  fprintf(f, " analysis=\"%s\"", analysis.c_str());
  fprintf(f, " time=\"%s\"", time.write().c_str());
  fprintf(f, " id=\"%d\"", id);
  fprintf(f,">\n");

  int t = tabs;
  if (t>-1) t++;

  size_t j;
  for (j = 0; j<libra_summary.size(); j++) libra_summary[j].write(f, t);
  for (j = 0; j<StPeter_analysis_summary.size(); j++) StPeter_analysis_summary[j].write(f, t);
  for (j = 0; j<decoy_analysis_summary.size(); j++) decoy_analysis_summary[j].write(f, t);
  for (j = 0; j<decoy_analysis.size(); j++) decoy_analysis[j].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</analysis_summary>\n");
}