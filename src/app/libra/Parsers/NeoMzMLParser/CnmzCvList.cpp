#include "CnmzCvList.h"

using namespace std;

void CnmzCvList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "cvList";
  if (cv.empty()) NMZerrMsg(el, "cv");
  if (count != cv.size()) NMZerrMsg(el, "count match number of cv");

  NMZprintTabs(f, tabs);
  fprintf(f, "<cvList");
  fprintf(f, " count=\"%d\"", count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<cv.size(); a++) cv[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f,"</cvList>\n");

}

