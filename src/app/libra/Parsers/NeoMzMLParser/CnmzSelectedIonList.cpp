#include "CnmzSelectedIonList.h"

using namespace std;

void CnmzSelectedIonList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "selectedIonList";
  if (selectedIon.empty()) NMZerrMsg(el, "selectedIon");
  if (count != selectedIon.size()) NMZerrMsg(el, "count match number of selectedIon");

  NMZprintTabs(f, tabs);
  fprintf(f, "<selectedIonList");
  fprintf(f, " count=\"%d\"", count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<selectedIon.size(); a++) selectedIon[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</selectedIonList>\n");

}

