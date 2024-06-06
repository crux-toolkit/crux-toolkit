#include "CnmzSourceFileList.h"

using namespace std;

void CnmzSourceFileList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "sourceFileList";
  if (sourceFile.empty()) NMZerrMsg(el, "sourceFile");
  if (count != sourceFile.size()) NMZerrMsg(el, "count match number of sourceFile");

  NMZprintTabs(f, tabs);
  fprintf(f, "<sourceFileList");
  fprintf(f, " count=\"%d\"", count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<sourceFile.size(); a++) sourceFile[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</sourceFileList>\n");

}

