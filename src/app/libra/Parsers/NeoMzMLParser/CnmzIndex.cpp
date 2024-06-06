#include "CnmzIndex.h"

using namespace std;

void CnmzIndex::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "index";
  if (name.empty()) NMZerrMsg(el, "name");

  NMZprintTabs(f, tabs);
  fprintf(f, "<index");
  fprintf(f, " name=\"%s\"", name.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<offset.size(); a++) offset[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</index>\n");

}

