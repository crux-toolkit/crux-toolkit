#include "CnmzSample.h"

using namespace std;

void CnmzSample::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "sample";
  if (id.empty()) NMZerrMsg(el, "id");

  NMZprintTabs(f, tabs);
  fprintf(f, "<sample");
  fprintf(f, " id=\"%s\"", id.c_str());
  if(!name.empty()) fprintf(f, " name=\"%s\"", name.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<referenceableParamGroupRef.size(); a++) referenceableParamGroupRef[a].write(f, t, iterative);
  for (size_t a = 0; a<cvParam.size(); a++) cvParam[a].write(f, t, iterative);
  for (size_t a = 0; a<userParam.size(); a++) userParam[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</sample>\n");

}

