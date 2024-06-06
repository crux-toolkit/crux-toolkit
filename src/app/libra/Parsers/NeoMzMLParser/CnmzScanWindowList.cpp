#include "CnmzScanWindowList.h"

using namespace std;

void CnmzScanWindowList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "scanWindowList";
  if (count!=(int)scanWindow.size()) NMZerrMsg(el, "count match number of scans");

  NMZprintTabs(f, tabs);
  fprintf(f, "<scanWindowList");
  fprintf(f, " count=\"%d\"", count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<scanWindow.size(); a++) scanWindow[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</scanWindowList>\n");

}
