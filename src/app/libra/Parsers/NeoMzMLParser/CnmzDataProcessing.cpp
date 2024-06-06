#include "CnmzDataProcessing.h"

using namespace std;

void CnmzDataProcessing::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "dataProcessing";
  if (processingMethod.empty()) NMZerrMsg(el, "processingMethod");
  if (id.empty()) NMZerrMsg(el, "id");

  NMZprintTabs(f, tabs);
  fprintf(f, "<dataProcessing");
  fprintf(f, " id=\"%s\"", id.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<processingMethod.size(); a++) processingMethod[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</dataProcessing>\n");

}

