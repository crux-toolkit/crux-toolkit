#include "CnmzReferenceableParamGroup.h"

using namespace std;

void CnmzReferenceableParamGroup::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "referenceableParamGroup";
  if (id.empty()) NMZerrMsg(el, "id");

  NMZprintTabs(f, tabs);
  fprintf(f, "<referenceableParamGroup");
  fprintf(f, " id=\"%s\"", id.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<cvParam.size(); a++) cvParam[a].write(f, t, iterative);
  for (size_t a = 0; a<userParam.size(); a++) userParam[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</referenceableParamGroup>\n");

}

