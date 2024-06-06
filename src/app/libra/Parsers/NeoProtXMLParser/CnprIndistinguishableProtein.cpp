#include "CnprIndistinguishableProtein.h"

using namespace std;

void CnprIndistinguishableProtein::write(FILE* f, int tabs){
  //required
  string el = "indistinguishable_protein";
  if (protein_name.empty()) NPRerrMsg(el, "protein_name");

  NPRprintTabs(f, tabs);
  fprintf(f, "<indistinguishable_protein");
  fprintf(f, " protein_name=\"%s\"", protein_name.c_str());
  if(parameter.empty() && annotation.empty()) {
    fprintf(f, "/>\n");
    return;
  } else fprintf(f,">\n");

  int t = tabs;
  if (t>-1) t++;

  size_t j;
  for (j = 0; j<parameter.size(); j++) parameter[j].write(f, t);
  for (j = 0; j<annotation.size(); j++) annotation[j].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</indistinguishable_protein>\n");

}
