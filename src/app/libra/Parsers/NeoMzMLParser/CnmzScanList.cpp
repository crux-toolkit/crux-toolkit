#include "CnmzScanList.h"

using namespace std;

void CnmzScanList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "scanList";
  if (scan.empty())NMZerrMsg(el, "scan");

  NMZprintTabs(f, tabs);
  fprintf(f, "<scanList");
  fprintf(f, " count=\"%d\"", count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<cvParam.size(); a++) cvParam[a].write(f, t, iterative);
  for (size_t a = 0; a<scan.size(); a++) scan[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</scanList>\n");

}
