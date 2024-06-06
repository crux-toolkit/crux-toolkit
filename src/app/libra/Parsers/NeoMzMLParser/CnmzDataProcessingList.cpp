#include "CnmzDataProcessingList.h"

using namespace std;

void CnmzDataProcessingList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "dataProcessingList";
  if (dataProcessing.empty()) NMZerrMsg(el, "dataProcessing");
  if (count != dataProcessing.size()) NMZerrMsg(el, "count match number of dataProcessing");

  NMZprintTabs(f, tabs);
  fprintf(f, "<dataProcessingList");
  fprintf(f, " count=\"%d\"", count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<dataProcessing.size(); a++) dataProcessing[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</dataProcessingList>\n");

}

