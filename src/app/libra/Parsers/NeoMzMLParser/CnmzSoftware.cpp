#include "CnmzSoftware.h"

using namespace std;

void CnmzSoftware::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "software";
  if (id.empty()) NMZerrMsg(el, "id");
  if (version.empty()) NMZerrMsg(el, "version");

  NMZprintTabs(f, tabs);
  fprintf(f, "<software");
  fprintf(f, " id=\"%s\"", id.c_str());
  fprintf(f, " version=\"%s\"", version.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<referenceableParamGroupRef.size(); a++) referenceableParamGroupRef[a].write(f, t, iterative);
  for (size_t a = 0; a<cvParam.size(); a++) cvParam[a].write(f, t, iterative);
  for (size_t a = 0; a<userParam.size(); a++) userParam[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</software>\n");

}

