#include "CnprStPeterAnalysisSummary.h"

using namespace std;

CnprStPeterAnalysisSummary::CnprStPeterAnalysisSummary(){
  probability = -1;
  FDR=-1;
  tolerance=-1;
  sampleLoad=-1;
}

void CnprStPeterAnalysisSummary::write(FILE* f, int tabs){
  //required
  string el = "StPeter_analysis_summary";
  if (version.empty()) NPRerrMsg(el, "version"); 
  if (probability<0) NPRerrMsg(el, "probability");
  if (FDR<0) NPRerrMsg(el, "FDR");
  if (tolerance<0) NPRerrMsg(el, "tolerance");

  NPRprintTabs(f, tabs);
  fprintf(f, "<StPeter_analysis_summary");
  fprintf(f, " version=\"%s\"", version.c_str());
  fprintf(f, " probability=\"%.4lf\"", probability);
  fprintf(f, " FDR=\"%.4lf\"", FDR);
  fprintf(f, " tolerance=\"%.4lf\"", tolerance);
  if (!degenerate_peptides.empty()) fprintf(f, " degenerate_peptides=\"%s\"", degenerate_peptides.c_str());
  if (sampleLoad>0) fprintf(f, " sampleLoad=\"%.4lf\"", sampleLoad);
  fprintf(f, "/>\n");

}
