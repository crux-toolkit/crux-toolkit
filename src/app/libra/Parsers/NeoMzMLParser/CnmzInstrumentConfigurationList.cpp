#include "CnmzInstrumentConfigurationList.h"

using namespace std;

void CnmzInstrumentConfigurationList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "instrumentConfigurationList";
  if (instrumentConfiguration.empty()) NMZerrMsg(el, "instrumentConfiguration");
  if (count != instrumentConfiguration.size()) NMZerrMsg(el, "count match number of instrumentConfiguration");

  NMZprintTabs(f, tabs);
  fprintf(f, "<instrumentConfigurationList");
  fprintf(f, " count=\"%d\"", count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<instrumentConfiguration.size(); a++) instrumentConfiguration[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</instrumentConfigurationList>\n");

}

