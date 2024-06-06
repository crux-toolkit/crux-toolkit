#include "CnmzSampleList.h"

using namespace std;

void CnmzSampleList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "sampleList";
  if (sample.empty()) NMZerrMsg(el, "sample");
  if (count != sample.size()) NMZerrMsg(el, "count match number of sample");

  NMZprintTabs(f, tabs);
  fprintf(f, "<sampleList");
  fprintf(f, " count=\"%d\"", count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<sample.size(); a++) sample[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</sampleList>\n");

}

