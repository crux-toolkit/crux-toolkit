#include "CnmzScanWindow.h"

using namespace std;

void CnmzScanWindow::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "scanWindow";

  NMZprintTabs(f, tabs);
  fprintf(f, "<scanWindow>\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<cvParam.size(); a++) cvParam[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</scanWindow>\n");

}
