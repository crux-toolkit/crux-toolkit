#include "CnprNSPInformation.h"

using namespace std;

void CnprNSPInformation::write(FILE* f, int tabs){
  //required
  string el = "nsp_information";
  if (neighboring_bin_smoothing.empty()) NPRerrMsg(el, "neighboring_bin_smoothing");

  NPRprintTabs(f, tabs);
  fprintf(f, "<nsp_information");
  fprintf(f, " neighboring_bin_smoothing=\"%s\"", neighboring_bin_smoothing.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;

  for(size_t j=0;j<nsp_distribution.size();j++) nsp_distribution[j].write(f,t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</nsp_information>\n");

}
