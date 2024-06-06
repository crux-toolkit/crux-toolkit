#include "CnmzIndexList.h"

using namespace std;

f_off_nmz CnmzIndexList::getIndexListOffset(){
  return fptr;
}

void CnmzIndexList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "indexList";
  if (index.empty()) NMZerrMsg(el, "index");

  NMZprintTabs(f, tabs);
  fptr = nmzftell(f);
  fprintf(f, "<indexList");
  fprintf(f, " count=\"%d\"", count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for(size_t a=0;a<index.size();a++) index[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</indexList>\n");

}

