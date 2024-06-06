#include "CnprDatasetDerivation.h"

using namespace std;

void CnprDatasetDerivation::write(FILE* f, int tabs){
  //required
  string el = "dataset_derivation";
  if (generation_no.empty()) NPRerrMsg(el, "generation_no");

  NPRprintTabs(f, tabs);
  fprintf(f, "<dataset_derivation");
  fprintf(f, " generation_no=\"%s\"", generation_no.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;

  //for (size_t j = 0; j<data_filter.size(); j++) data_filter[j].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</dataset_derivation>\n");

}
