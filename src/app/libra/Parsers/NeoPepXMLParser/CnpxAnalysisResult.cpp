#include "CnpxAnalysisResult.h"

using namespace std;

CnpxAnalysisResult::CnpxAnalysisResult() {
  analysis.clear();
  id = 0;
}

void CnpxAnalysisResult::write(FILE* f, int tabs) {
  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<analysis_result analysis=\"%s\"", analysis.c_str());
  if (id > 0) fprintf(f, " id=\"%d\"", id);
  fprintf(f, ">\n");

  for(size_t i=0;i<parameter.size();i++) parameter[i].write(f);

  if (peptide_prophet_result.present()) peptide_prophet_result.write(f);
  if (interprophet_result.present()) interprophet_result.write(f);
  if (libra_result.present()) libra_result.write(f,t);
  if (pepxmlquant_result.present()) pepxmlquant_result.write(f);
  if (quantic_result.present()) quantic_result.write(f);
  if (expresslabelfree_result.present()) expresslabelfree_result.write(f);
  
  for (size_t i = 0; i < ptmprophet_result.size(); i++) ptmprophet_result[i].write(f);

  NPXprintTabs(f, tabs);
  fprintf(f, "</analysis_result>\n");
}