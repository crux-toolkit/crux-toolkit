#include "CnmzSourceFile.h"

using namespace std;

void CnmzSourceFile::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "sourceFile";
  if (id.empty()) NMZerrMsg(el, "id");
  if (location.empty()) NMZerrMsg(el, "location");
  if (name.empty()) NMZerrMsg(el, "name");

  NMZprintTabs(f, tabs);
  fprintf(f, "<sourceFile");
  fprintf(f, " id=\"%s\"",id.c_str());
  fprintf(f, " name=\"%s\"", name.c_str());
  fprintf(f, " location=\"%s\"", location.c_str());
  fprintf(f,">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<referenceableParamGroupRef.size(); a++) referenceableParamGroupRef[a].write(f, t, iterative);
  for (size_t a = 0; a<cvParam.size(); a++) cvParam[a].write(f, t, iterative);
  for (size_t a = 0; a<userParam.size(); a++) userParam[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</sourceFile>\n");

}

