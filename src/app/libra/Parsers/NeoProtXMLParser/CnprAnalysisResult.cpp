#include "CnprAnalysisResult.h"

using namespace std;

CnprAnalysisResult::CnprAnalysisResult(){
  id=0;
}

void CnprAnalysisResult::write(FILE* f, int tabs){
  //required
  string el = "analysis_result";
  if (analysis.empty()) NPRerrMsg(el, "analysis");

  NPRprintTabs(f, tabs);
  fprintf(f, "<analysis_result");
  fprintf(f, " analysis=\"%s\"", analysis.c_str());
  if(id>0) fprintf(f, " id=\"%d\"", id);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;

  size_t j;
  for (j = 0; j<libra_result.size(); j++) libra_result[j].write(f, t);
  for(j=0;j<StPeterQuant.size();j++) StPeterQuant[j].write(f,t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</analysis_result>\n");

}
