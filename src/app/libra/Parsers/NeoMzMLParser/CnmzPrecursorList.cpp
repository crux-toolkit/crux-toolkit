#include "CnmzPrecursorList.h"

using namespace std;

void CnmzPrecursorList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "precursorList";
  if (precursor.empty()) NMZerrMsg(el, "precursor");
  if (count != precursor.size()) NMZerrMsg(el, "count match number of precursor");

  NMZprintTabs(f, tabs);
  fprintf(f, "<precursorList");
  fprintf(f, " count=\"%d\"", count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<precursor.size(); a++) precursor[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</precursorList>\n");

}

