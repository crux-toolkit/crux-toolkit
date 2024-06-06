#include "CnprDecoyAnalysis.h"

using namespace std;



void CnprDecoyAnalysis::write(FILE* f, int tabs){
  //no requirements

  NPRprintTabs(f, tabs);
  fprintf(f, "<decoy_analysis>\n");

  int t = tabs;
  if (t>-1) t++;

  size_t j;
  for (j = 0; j<point.size(); j++) point[j].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</decoy_analysis>\n");
}