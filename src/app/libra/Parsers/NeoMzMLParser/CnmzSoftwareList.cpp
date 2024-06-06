#include "CnmzSoftwareList.h"

using namespace std;

void CnmzSoftwareList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "softwareList";
  if (software.empty()) NMZerrMsg(el, "software");
  if (count != software.size()) NMZerrMsg(el, "count match number of software");

  NMZprintTabs(f, tabs);
  fprintf(f, "<softwareList");
  fprintf(f, " count=\"%d\"", count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<software.size(); a++) software[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</softwareList>\n");

}

