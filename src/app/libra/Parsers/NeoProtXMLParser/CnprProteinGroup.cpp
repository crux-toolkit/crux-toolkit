#include "CnprProteinGroup.h"

using namespace std;

CnprProteinGroup::CnprProteinGroup(){
  probability=-1;
}

void CnprProteinGroup::write(FILE* f, int tabs){
  //required
  string el = "protein_group";
  if (group_number.empty()) NPRerrMsg(el, "group_number");
  if (probability<0)NPRerrMsg(el, "probability");
  if (protein.empty())NPRerrMsg(el, "protein");

  int t = tabs;
  if (t>-1) t++;

  NPRprintTabs(f, tabs);
  fprintf(f, "<protein_group");
  fprintf(f, " group_number=\"%s\"", group_number.c_str());
  if(!pseudo_name.empty()) fprintf(f, " pseudo_name=\"%s\"", pseudo_name.c_str());
  fprintf(f, " probability=\"%.4lf\"", probability);
  fprintf(f, ">\n");

  for(size_t j=0;j<protein.size();j++) protein[j].write(f,t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</protein_group>\n");

}

