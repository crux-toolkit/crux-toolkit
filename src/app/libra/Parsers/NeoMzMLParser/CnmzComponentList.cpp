#include "CnmzComponentList.h"

using namespace std;

void CnmzComponentList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "componentList";
  if (source.empty()) NMZerrMsg(el, "source");
  if (analyzer.empty()) NMZerrMsg(el, "analyzer");
  if (detector.empty()) NMZerrMsg(el, "detectir");
  if (count != source.size()+analyzer.size()+detector.size()) NMZerrMsg(el, "count match number of source+analyzer+detector");

  NMZprintTabs(f, tabs);
  fprintf(f, "<componentList");
  fprintf(f, " count=\"%d\"", count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<source.size(); a++) source[a].write(f, t, iterative);
  for (size_t a = 0; a<analyzer.size(); a++) analyzer[a].write(f, t, iterative);
  for (size_t a = 0; a<detector.size(); a++) detector[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</componentList>\n");

}

